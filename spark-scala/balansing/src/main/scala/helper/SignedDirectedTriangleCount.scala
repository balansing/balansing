package helper

import java.util.StringTokenizer

import scala.io.Source

object SignedDirectedTriangleCount {

  def countTriangles(inputPath: String): Array[Long] ={
    val edges = Source.fromFile(inputPath).getLines().filter(_.nonEmpty).map{line =>

      val st = new StringTokenizer(line)
      val u = st.nextToken().toInt
      val v = st.nextToken().toInt
      val plus = st.nextToken().toInt > 0

      (u, v, plus)
    }.toSeq

    countTriangles(edges)
  }

  def countTriangles(edges : Seq[(Int, Int, Boolean)]): Array[Long] ={

    val outNeighbors = edges.groupBy(_._1) // group by u
                            .mapValues { seq =>
                              seq.map { case (u, v, plus) => (v, plus) }.sorted.toArray
                            }
                            .withDefaultValue(Array.empty[(Int, Boolean)])

    val inNeighbors = edges.groupBy(_._2) // group by v
                           .mapValues { seq =>
                             seq.map { case (u, v, plus) => (u, plus) }.sorted.toArray
                           }
                           .withDefaultValue(Array.empty[(Int, Boolean)])


    /*
     * type 0 ~ type 7: non-cycle triangle types
     * The type of a non-cycle triangle is defined by the sign of pivot, left, and right edges
     * * pivot edge: the edge whose both end nodes have out-going edges.
     * * left edge: the edge having the same source node with the pivot edge
     * * right edge: the edge whose source node is the target node of the pivot edge
     * type 0: +++
     * type 1: ++-
     * type 2: +-+
     * type 3: +--
     * type 4: -++
     * type 5: -+-
     * type 6: --+
     * type 7: ---
     *
     * type 8 ~ type 11: cycle triangle types
     * The type of a cycle triangle is defined by the number of plus signs.
     * type 8: +++
     * type 9: ++-
     * type 10: +--
     * type 11: ---
     *
     */
    val numTrianglesByType = new Array[Long](12)


    // count non-cycle triangles
    edges.foreach { case (u, v, plus) =>

      val uN = outNeighbors(u)
      val vN = outNeighbors(v)

      //p: plus, m: minus
      val (pp, pm, mp, mm) = countIntersectBySign(uN, vN)

      if(plus){
        numTrianglesByType(0) += pp
        numTrianglesByType(1) += pm
        numTrianglesByType(2) += mp
        numTrianglesByType(3) += mm
      }
      else{
        numTrianglesByType(4) += pp
        numTrianglesByType(5) += pm
        numTrianglesByType(6) += mp
        numTrianglesByType(7) += mm
      }

    }

    // count cycle triangles
    // each triangle is counted 3 times.
    edges.foreach{ case (u, v, plus) =>

      val uN = inNeighbors(u)
      val vN = outNeighbors(v)

      //p: plus, m: minus
      val (pp, pm, mp, mm) = countIntersectBySign(uN, vN)

      if(plus){
        numTrianglesByType(8) += pp
        numTrianglesByType(9) += pm
        numTrianglesByType(9) += mp
        numTrianglesByType(10) += mm
      }
      else{
        numTrianglesByType(9) += pp
        numTrianglesByType(10) += pm
        numTrianglesByType(10) += mp
        numTrianglesByType(11) += mm
      }

    }

    //correct the cycle triangle numbers
    numTrianglesByType(8) /= 3
    numTrianglesByType(9) /= 3
    numTrianglesByType(10) /= 3
    numTrianglesByType(11) /= 3

    numTrianglesByType
  }


  def countIntersectBySign(uN: Array[(Int, Boolean)], vN: Array[(Int, Boolean)]): (Long, Long, Long, Long) = {

    var (pp, pm, mp, mm) = (0L,0L,0L,0L)
    var uCur = 0
    var vCur = 0
    val uD = uN.length
    val vD = vN.length

    while (uCur < uD && vCur < vD) {
      if (uN(uCur)._1 < vN(vCur)._1) {
        uCur += 1
      } else if (vN(vCur)._1 < uN(uCur)._1)
        vCur += 1
      else {

        if(uN(uCur)._2) {
          if (vN(vCur)._2) pp += 1
          else pm += 1
        }
        else{
          if (vN(vCur)._2) mp += 1
          else mm += 1
        }

        uCur += 1
        vCur += 1
      }
    }

    return (pp, pm, mp, mm)
  }
}
