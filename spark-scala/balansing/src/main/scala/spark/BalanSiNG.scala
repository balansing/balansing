package spark

import org.apache.log4j.Logger
import org.apache.spark.rdd.RDD
import org.apache.spark.{HashPartitioner, SparkConf, SparkContext}

import scala.collection.mutable
import scala.util.Random

object BalanSiNG {

  var rand: Random = _

  val logger: Logger = Logger.getLogger(getClass)

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
    val split_ratio = args(7).toDouble
    val alpha = args(8).toDouble
    val noise = args(9).toDouble
    val output = args(10)
    val randSeed = args(11).toInt

    rand = new Random(randSeed)

    logger.info(f"level: $level, numEdges: $numEdges, numTasks: $numTasks, " +
      f"(a,b,c,d): $abcd")

    val conf = new SparkConf().setAppName(f"[RMAT] lv: $level, E: $numEdges, abcd: $abcd") //.setMaster("local")

    val sc = new SparkContext(conf)

    run(level, numEdges, numTasks, abcd, split_ratio, alpha, noise, sc, randSeed)
      .map { case (u, v, p) => u + "\t" + v + "\t" + (if (p) "+1" else "-1") }
      .saveAsTextFile(output)

    sc.stop()
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

    rand = new Random(randSeed)

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

    val rmat = new BalanSiNGSingle(abcd, split_ratio, alpha, noiseList.value)

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
      val l_rmat = new BalanSiNGSingle(abcd, split_ratio, alpha, noiseList.value)
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
