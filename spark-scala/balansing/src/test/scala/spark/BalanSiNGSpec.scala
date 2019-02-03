package spark

import helper.SignedDirectedTriangleCount
import org.apache.spark.{SparkConf, SparkContext}
import org.scalatest.{FlatSpec, Matchers}

class BalanSiNGSpec extends FlatSpec with Matchers {

  it should "run without an error" in {

    val conf = new SparkConf().setAppName("BalanSiNG-test").setMaster("local[1]")
    val sc = new SparkContext(conf)

    val level = 11
    val numEdges = 36000
    val numTasks = 1
    val A = (0.57,0.19,0.19,0.05)
    val alpha = 0.2
    val gamma = 0.1

    val data = BalanSiNG.run(level,numEdges,numTasks, A, 0.0, alpha, gamma, sc).collect()
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
    println(s"Ratio of negative edges: $unbal_ratio")

  }
}