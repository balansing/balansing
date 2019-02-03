package spark

import org.apache.commons.math3.distribution.MultivariateNormalDistribution
import spark.BalanSiNG.Region

import scala.collection.mutable

/**
  * rmat generator
  *
  * @param abcd the probabilities
  * @param noiseList a list of noises applied to each level
  */
class BalanSiNGSingle(abcd: (Double, Double, Double, Double), gamma: Double, alpha: Double, noiseList: Array[Double]) {


  private val probs = noiseList.map{ mu =>
    (abcd._1 - 2 * mu * abcd._1 / (abcd._1 + abcd._4),
      abcd._2 + mu,
      abcd._3 + mu,
      abcd._4 - 2 * mu * abcd._4 / (abcd._1 + abcd._4))
  }


  /**
    * generate edges within a given region
    * @param region to process
    * @return list of edges in the region
    */
  def generate(region: Region): Stream[(Long, Long, Boolean)] = {
    if(region.numEdges < BalanSiNGSingle.THRESHOLD_REGION_EDGE_NUM)
      (0L until region.numEdges).map { _ => makeOneEdge(region) }.distinct.toStream
    else
      splitRegion(region).toStream.flatMap(generate)
  }

  /**
    * generate one edge
    * @param region
    * @return
    */
  def makeOneEdge(region: Region): (Long, Long, Boolean) = {

    val A = probs(region.level)

    val cumul: (Double, Double, Double) = (A._1, A._1 + A._2, A._1 + A._2 + A._3)

    var x = region.srcMin
    var y = region.dstMin

    var plusProb = region.plus_joint

    //var p = Math.random
    var p = BalanSiNG.rand.nextDouble

    for (l <- region.level to 0 by -1) {

      plusProb =
        if (p > cumul._3) {
          x += (1L << l)
          y += (1L << l)
          p = (p - cumul._3) / A._4
          //(plusProb * (1 - gamma) + (1 - plusProb) * gamma) * (1 - alpha) + alpha
          (plusProb * (1 - gamma) + (1 - plusProb) * gamma) * (alpha) + (1 - alpha)
        }
        else if (p > cumul._2){
          x += (1L << l)
          p = (p - cumul._2) / A._3
          (plusProb * gamma + (1 - plusProb) * (1 - gamma)) * (alpha) + (1 - alpha)
        }
        else if (p > cumul._1){
          y += (1L << l)
          p = (p - cumul._1) / A._2
          (plusProb * gamma + (1 - plusProb) * (1 - gamma)) * (alpha) + (1 - alpha)
        }
        else{
          p = p / A._1
          (plusProb * (1 - gamma) + (1 - plusProb) * gamma) * (alpha) + (1 - alpha)
        }

    }

    // stochastic sign decision
    // if (Math.random < plusProb) (x, y, true)
    if (BalanSiNG.rand.nextDouble < plusProb) (x, y, true)
    else (x, y, false)

  }

  /**
    * split a region
    * @param region to split
    * @return list of splitted regions
    */
  def splitRegion(region: Region): Array[Region] = {


    val t = probs(region.level)

    val splitEdgeNums = multinomial(region.numEdges, Array(t._1, t._2, t._3))
    val lastSplitEdgeNum = region.numEdges - splitEdgeNums.sum

    val srcMid = (region.srcMin + region.srcMax) >> 1
    val dstMid = (region.dstMin + region.dstMax) >> 1

    var maxSplitNumNodes = (region.srcMax - srcMid) * (region.srcMax - srcMid)

    if(region.srcMax - srcMid > 0 && maxSplitNumNodes == 0){
      maxSplitNumNodes = Long.MaxValue
    }

    val res = new mutable.ArrayBuffer[Region]()

    if(splitEdgeNums(0) > 0){
      val plus_joint = region.plus_joint * (1 - gamma) + (1 - region.plus_joint) * gamma
      res += Region(region.srcMin, region.dstMin, srcMid, dstMid,
        Math.min(splitEdgeNums(0), maxSplitNumNodes), region.level - 1, plus_joint * alpha + 1 - alpha)
    }

    if(splitEdgeNums(1) > 0) {
      val plus_joint = region.plus_joint * gamma + (1 - region.plus_joint) * (1 - gamma)
      res += Region(region.srcMin, dstMid + 1, srcMid, region.dstMax,
        Math.min(splitEdgeNums(1), maxSplitNumNodes), region.level - 1, plus_joint * alpha + 1 - alpha)
    }

    if(splitEdgeNums(2) > 0) {
      val plus_joint = region.plus_joint * gamma + (1 - region.plus_joint) * (1 - gamma)
      res += Region(srcMid + 1, region.dstMin, region.srcMax, dstMid,
        Math.min(splitEdgeNums(2), maxSplitNumNodes), region.level - 1, plus_joint * alpha + 1 - alpha)
    }

    if(lastSplitEdgeNum > 0) {
      val plus_joint = region.plus_joint * (1 - gamma) + (1 - region.plus_joint) * gamma
      res += Region(srcMid + 1, dstMid + 1, region.srcMax, region.dstMax,
        Math.min(lastSplitEdgeNum, maxSplitNumNodes), region.level - 1, plus_joint * alpha + 1 - alpha)
    }


    res.toArray
  }





  /**
    * get a sample following a multinomial distribution
    *
    * @param n the number of trials
    * @param p list of probability
    * @return a sample following a multinomial distribution
    */
  def multinomial(n: Long, p: Array[Double]): Array[Long] ={

    val nM = Array.ofDim[Double](p.length, p.length)
    val np = p.map(_*n)

    for (i <- p.indices; j <- p.indices) {
      nM(i)(j) = if (i == j) p(i) else 0
      nM(i)(j) -= p(i) * p(j)
      nM(i)(j) *= n
    }

    val distribution = new MultivariateNormalDistribution(np, nM)
    var sample: Array[Long] = null

    do{
      sample = distribution.sample().map(Math.round)
    }while(sample.exists(_<0))


    sample
  }

}

object BalanSiNGSingle{
  // if a region contains less than thresholdRegionEdgeNum, this region will not be splitted.
  val THRESHOLD_REGION_EDGE_NUM = 10000
}
