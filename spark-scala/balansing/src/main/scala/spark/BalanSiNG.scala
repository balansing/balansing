package spark

import helper.SignedDirectedTriangleCount
import org.apache.log4j.Logger
import org.apache.spark.rdd.RDD
import org.apache.spark.{HashPartitioner, SparkConf, SparkContext}

import scala.collection.mutable
import scala.util.Random

import org.apache.log4j.Logger
import org.apache.log4j.Level


object BalanSiNG {

  var rand: Random = _

  val logger: Logger = Logger.getLogger(getClass)

  Logger.getLogger("org").setLevel(Level.OFF)
  Logger.getLogger("akka").setLevel(Level.OFF)

  /**
    * The main entry point
    * @param args [0]: level of recursion. [1]: the number of edges. [2]: degree of parallelism (the
    *             number of tasks to run). [3], [4], [5], [6]: the probabilities of quarters. [7]: noise
    *             for smoothing
    */
  def main(args: Array[String]): Unit = {

    val level = args(0).toInt
    val numEdges = args(1).toLong
    val numTasks = args(2).toInt
    val abcd = (args(3).toDouble, args(4).toDouble, args(5).toDouble, args(6).toDouble)
    val alpha = args(7).toDouble
    val noise = args(8).toDouble
    val output = args(9)
    val randSeed = args(10).toInt

    rand = new Random(randSeed)

    logger.info(f"level: $level, numEdges: $numEdges, numTasks: $numTasks, " +
      f"(p11,m12,m21,p22): $abcd")

    val conf = new SparkConf().setAppName(f"[BalanSiNG] lv: $level, E: $numEdges, abcd: $abcd") //.setMaster("local")

    val sc = new SparkContext(conf)

    val split_ratio = 0
    val res = run(level, numEdges, numTasks, abcd, split_ratio, alpha, noise, sc, randSeed)

    res.map { case (u, v, p) => u + "\t" + v + "\t" + (if (p) "+1" else "-1") }
       .saveAsTextFile(output)

    val data = res.collect()

    val pos = data.filter(_._3 == true).map{case (x, y, _) => (x, y)}
    val neg = data.filter(_._3 == false).map{case (x, y, _) => (x, y)}
    val numTrianglesBySign = SignedDirectedTriangleCount.countTriangles(data.map{case (u, v, plus) => (u.toInt, v.toInt, plus)}.toSeq)

    val balanced = numTrianglesBySign(0) + numTrianglesBySign(3) + numTrianglesBySign(5) + numTrianglesBySign(6) +
      numTrianglesBySign(8)+ numTrianglesBySign(10)
    val unbalanced = numTrianglesBySign(1) + numTrianglesBySign(2) + numTrianglesBySign(4) + numTrianglesBySign(9) +
      numTrianglesBySign(7) + numTrianglesBySign(11)


    val sign_total = (pos.size + neg.size).toDouble
    val pos_ratio = pos.size / sign_total
    val neg_ratio = neg.size / sign_total

    val tri_total = (balanced + unbalanced).toDouble
    val bal_ratio = balanced / tri_total
    val unbal_ratio =  unbalanced / tri_total

    sc.stop()

    println("Ratio of positive edges: " + pos_ratio)
    println("Ratio of negative edges: " + neg_ratio)
    println(s"Ratio of balanced triangles: $bal_ratio")
    println(s"Ratio of unbalanced triangles: $unbal_ratio")
  }


  /**
    * submit the spark job.
    *
    * @param level of recursion
    * @param numEdges the number of edges
    * @param numTasks degree of parallelism (the number of tasks to run)
    * @param abcd the probabilities of quarters.
    * @param gamma
    * @param alpha
    * @param noise for smoothing
    * @param sc spark context
    * @return an RDD containing edges of the generated graph
    */
  def run(level: Int, numEdges: Long, numTasks: Int,
          abcd: (Double, Double, Double, Double), split_ratio: Double, alpha: Double, noise: Double, sc: SparkContext, randSeed: Int): RDD[(Long, Long, Boolean)] = {

    // rand = new Random(randSeed)
    rand = new Random()

    val regionQueue = new mutable.PriorityQueue[Region]()(Ordering.by(regionOrder))
    val regionSlotQueue = new mutable.PriorityQueue[RegionSlot]()(Ordering.by(regionSlotOrder))

    for (i <- 0 until numTasks) {
      val numEdgeAssigned = (numEdges * (i + 1) / numTasks) - (numEdges * i / numTasks)
      regionSlotQueue.enqueue(RegionSlot(numEdgeAssigned, mutable.Buffer[Region]()))
    }

    val numNodes = 1L << (level + 1)

    //val noiseList = sc.broadcast[Array[Double]]((0L to level).map(_ => Math.random() * 2 * noise - noise).toArray)
    val noiseList = sc.broadcast[Array[Double]]((0L to level).map(_ => rand.nextDouble * 2 * noise - noise).toArray)

    regionQueue.enqueue(Region(0, 0, numNodes - 1, numNodes - 1, numEdges, level, 1.0))

    val rmat = new BalanSiNGSingle(abcd, split_ratio, alpha, noiseList.value, randSeed)

    while (regionQueue.nonEmpty) {
      val region = regionQueue.dequeue
      val slot = regionSlotQueue.dequeue

      if(region.numEdges <= BalanSiNGSingle.THRESHOLD_REGION_EDGE_NUM || region.numEdges <= slot.edgeCapacity){
        slot.add(region)
      }
      else{
        rmat.splitRegion(region).foreach{r => regionQueue.enqueue(r)}
      }

      regionSlotQueue.enqueue(slot)
    }

    val regionSlotRDD = sc.parallelize(regionSlotQueue.toSeq.zipWithIndex.map{case (x,y) => (y,x)})
      .partitionBy(new HashPartitioner(numTasks))

    regionSlotRDD.flatMap{ case (_, x) =>
      val l_rmat = new BalanSiNGSingle(abcd, split_ratio, alpha, noiseList.value, randSeed)
      x.regions.toStream.flatMap(l_rmat.generate)
    }

  }

  case class Region(srcMin: Long, dstMin: Long, srcMax: Long, dstMax: Long, numEdges: Long, level: Int, plus_joint: Double)

  case class RegionSlot(var edgeCapacity: Long, regions: mutable.Buffer[Region]){

    /**
      * add a region to this region slot
      * @param r region to add
      * @return true if r is added or false
      */
    def add(r: Region): Unit ={
      regions.append(r)
      edgeCapacity -= r.numEdges
    }
  }

  def regionOrder(r: Region): Long = r.numEdges

  def regionSlotOrder(s: RegionSlot): Long = s.edgeCapacity
}
