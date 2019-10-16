package spark

import helper.SignedDirectedTriangleCount
import org.apache.spark.{SparkConf, SparkContext}
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.mutable.ListBuffer

class BalanSiNGSpec extends FlatSpec with Matchers {

  def mean(list:List[Double]):Double =
    if(list.isEmpty) 0 else list.sum/list.size

  def stddev(list:List[Double]):Double =
    Math.sqrt(mean(list.map(x => Math.pow(x, 2))) - Math.pow(mean(list), 2))

  it should "run without an error" in {

    val conf = new SparkConf().setAppName("BalanSiNG-test").setMaster("local[1]")
    val sc = new SparkContext(conf)

    val level = 17
    val numEdges = 841372
    val numTasks = 1
    val A = (0.57,0.19,0.19,0.05)
    val alpha = 0.15
    val gamma = 0.1
    val randSeed = 0
    var pos_ratios = ListBuffer[Double]()
    var bal_ratios = ListBuffer[Double]()

    //for(t <- 0 to 100) {
    //println("Start trial " + t)
      val data = BalanSiNG.run(level, numEdges, numTasks, A, 0.0, alpha, gamma, sc, randSeed).collect()

      val pos = data.filter(_._3 == true).map { case (x, y, _) => (x, y) }
      val neg = data.filter(_._3 == false).map { case (x, y, _) => (x, y) }

      val t1 = System.nanoTime
      val numTrianglesBySign = SignedDirectedTriangleCount.countTriangles(data.map { case (u, v, plus) => (u.toInt, v.toInt, plus) }.toSeq)
      val duration = (System.nanoTime - t1) / 1e9d
      println("count_time: " + duration)


      val balanced = numTrianglesBySign(0) + numTrianglesBySign(3) + numTrianglesBySign(5) + numTrianglesBySign(6) +
        numTrianglesBySign(8) + numTrianglesBySign(10)
      val unbalanced = numTrianglesBySign(1) + numTrianglesBySign(2) + numTrianglesBySign(4) + numTrianglesBySign(9) +
        numTrianglesBySign(7) + numTrianglesBySign(11)

      val sign_total = (pos.size + neg.size).toDouble
      val pos_ratio = pos.size / sign_total
      val neg_ratio = neg.size / sign_total

      val tri_total = (balanced + unbalanced).toDouble
      val bal_ratio = balanced / tri_total
      val unbal_ratio = unbalanced / tri_total

      pos_ratios += pos_ratio
      bal_ratios += bal_ratio
    //}

    sc.stop()

    val avg_pos_ratio = mean(pos_ratios.toList)
    val std_pos_ratio = stddev(pos_ratios.toList)
    val avg_bal_ratio = mean(bal_ratios.toList)
    val std_bal_ratio = stddev(bal_ratios.toList)
    println("[pos_ratio] mean: " + avg_pos_ratio + ", stddev: " + std_pos_ratio)
    println("[bal_ratio] mean: " + avg_bal_ratio + ", stddev: " + std_bal_ratio)

    /*println("Ratio of positive edges: " + pos_ratio)
    println("Ratio of negative edges: " + neg_ratio)
    println(s"Ratio of balanced triangles: $bal_ratio")
    println(s"Ratio of unbalanced triangles: $unbal_ratio")*/
  }
}