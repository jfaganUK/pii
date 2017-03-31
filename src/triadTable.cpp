#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::depends("RcppArmadillo")]]
//' @export
// [[Rcpp::export]]

DataFrame triadTable(NumericVector edgeDistance, IntegerMatrix shortPaths, IntegerMatrix triads, CharacterVector vertices, NumericVector edgevalence){
  IntegerVector edgeDistanceDims = edgeDistance.attr("dim");
  arma::cube cubeEdgeDistance(edgeDistance.begin(), edgeDistanceDims[0], edgeDistanceDims[1], edgeDistanceDims[2], false);
  DataFrame triadtable;
  CharacterVector triadID(0);
  NumericVector nodeID(0);
  CharacterVector direction(0);
  NumericVector valence(0);
  NumericVector distance(0);

  int node1, node2, node3;
  String triname, edge1_2, edge1_3, edge2_3;
  int ed, dis1, dis2, dis3;
  String ref;
  int truev, diff;

  for(int i = 0; i < triads.nrow(); i++) {
    node1 = triads(i,0);
    node2 = triads(i,1);
    node3 = triads(i,2);
    std::string snode1;
    std::string snode2;
    std::string snode3;
    std::stringstream out;
    out << node1;
    snode1 = out.str();
    out.str(std::string());
    out << node2;
    snode2 = out.str();
    out.str(std::string());
    out << node3;
    snode3 = out.str();
    triname = snode1 + "-" + snode2 + "-" + snode3;
    edge1_2 = snode1 + "-" + snode2;
    edge1_3 = snode1 + "-" + snode3;
    edge2_3 = snode2 + "-" + snode3;
    for(int j = 0; j < vertices.length(); j++){
      dis1 = cubeEdgeDistance((node1 - 1), (node2 - 1), j);
      dis2 = cubeEdgeDistance((node1 - 1), (node3 - 1), j);
      dis3 = cubeEdgeDistance((node2 - 1), (node3 - 1), j);
      if (dis1 == dis2 & dis1 == dis3) {
        int nodedis1_2 = shortPaths[node1, node2];
        int nodedis1_3 = shortPaths[node1, node3];
        int nodedis2_3 = shortPaths[node2, node3];
        if (nodedis1_2 !=  nodedis1_3) {
          if (nodedis1_2 != nodedis2_3) {
            ref = edge1_2;
          }else{
            ref = edge1_3;
          }
        }else{
          ref = edge2_3;
        }
        //NumericVector valence_set = edgevalence[triads(i, 0)];
        if (ref == edge1_2) {
          truev = edgevalence[triads(i, 0)];
        }
        if (ref == edge1_3) {
          truev = edgevalence[triads(i, 2)];
        }
        if (ref == edge2_3) {
          truev = edgevalence[triads(i, 1)];
        }
        triadID.push_back(triname);
        nodeID.push_back(j+1);
        direction.push_back("IN ");
        valence.push_back(truev);
        distance.push_back(dis1);
      }
      else{
        if (dis1 != dis2) {
          if (dis1 != dis3) {
            diff = dis1;
          }else{
            diff = dis2;
          }
        }else{
          diff = dis3;
        }

        //NumericVector valence_set = edgevalence[triads(i)];
        if (diff == dis1) {
          truev = edgevalence[triads(i, 0)];
        }
        if (diff == dis2) {
          truev = edgevalence[triads(i, 2)];
        }
        if (diff == dis3) {
          truev = edgevalence[triads(i, 1)];
        }
        triadID.push_back(triname);
        nodeID.push_back(j+1);
        direction.push_back("OUT");
        valence.push_back(truev);
        distance.push_back(diff);
      }
    }
  }

  triadtable = DataFrame::create(_["triadID"]= triadID, _["nodeID"]= nodeID, _["direction"]= direction, _["valence"]= valence, _["distance"]= distance);
  return triadtable;
}
