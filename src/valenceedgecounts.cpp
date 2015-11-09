#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
List valenceEdgeCounts(IntegerMatrix edgeDistance, NumericVector edgevalence, int maxDistance) {
  int nEdges = edgeDistance.nrow(), nNodes = edgeDistance.ncol();
  NumericVector negCount(maxDistance+1), posCount(maxDistance+1);
  NumericMatrix pos(maxDistance+1, nNodes), neg(maxDistance+1, nNodes);
  //List crosslist;

  for(int i = 0; i < nNodes; i++) {
    for(int j = 0; j < nEdges; j++) {
      if(edgevalence[j] < 0) {
        negCount[edgeDistance(j,i)]++;
      }
      if(edgevalence[j] > 0) {
        posCount[edgeDistance(j,i)]++;
      }
    }
    for(int k = 0; k < maxDistance+1; k++){
      pos(k, i) = posCount[k];
      neg(k, i) = negCount[k];
      posCount[k] = 0;
      negCount[k] = 0;
    }
  }
  return List::create(Named("pos") = pos,
                           Named("neg") = neg);
}
