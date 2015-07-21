#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;


//' Calculate the PII scores after the edge distances are calculated
//'
//' @param edgeDistance An integer matrix of distances from the nodes to the edges
//' @param valence The postive / negative valence of each of the edges.
//' @param piiBeta The beta attenuation value (ideally -1.0 to -0.01)
//' @param piiX The normalization factor based on the maximum degree.
//' @param maxDistance The maximum distance in the edgeDistance matrix (or some user-provided max distance for very large networks)
//' @export
// [[Rcpp::export]]
NumericVector piiCalc(IntegerMatrix edgeDistance, NumericVector edgevalence, double piiBeta, double piiX, int maxDistance) {
	int nEdges = edgeDistance.nrow(), nNodes = edgeDistance.ncol();
	IntegerVector negCount(maxDistance+1), posCount(maxDistance+1);
	NumericVector piiBetaVector(maxDistance + 1);
	NumericVector piIndex(nNodes);


	// Initialize piiBetaVector
	for(int k = 0; k <= maxDistance; k++) {
		piiBetaVector[k] = pow(piiBeta, k);
	}

	for(int i = 0; i < nNodes; i++) {
		for(int j = 0; j < nEdges; j++) {
			if(edgevalence[j] < 0) {
				negCount[edgeDistance(j,i)]++;
			}
			if(edgevalence[j] > 0) {
				posCount[edgeDistance(j,i)]++;
			}
		}
		for(int k = 0; k <= maxDistance; k++) {
			piIndex[i] += piiBetaVector[k] * (pow(posCount[k], piiX) - pow(negCount[k], piiX));
			negCount[k] = 0;
			posCount[k] = 0;
		}
	}

	return piIndex;
}

//' @export
// [[Rcpp::export]]
NumericVector piiTriadicCalc(IntegerMatrix edgeDistance, NumericVector edgevalence, double piiBeta, double piiX, int maxDistance, DataFrame triadTable, double piiDelta) {
	int nEdges = edgeDistance.nrow(), nNodes = edgeDistance.ncol();
	IntegerVector negCount(maxDistance+1), posCount(maxDistance+1);
  IntegerVector negOutCount(maxDistance+1), negInCount(maxDistance+1),
                posOutCount(maxDistance+1), posInCount(maxDistance+1);
	NumericVector piiBetaVector(maxDistance + 1);
	NumericVector piIndex(nNodes);
	double triadicPart;
	int ed; // a holder for the current edge distance

  CharacterVector triadID = triadTable["triadID"];
  CharacterVector direction = triadTable["direction"];
  IntegerVector valence = triadTable["valence"];
  IntegerVector distance = triadTable["distance"];
  IntegerVector nodeNum = triadTable["nodeID"];

	// Initialize piiBetaVector
	for(int k = 0; k <= maxDistance; k++) {
		piiBetaVector[k] = pow(piiBeta, k);
	}

	for(int i = 0; i < nNodes; i++) {
		for(int j = 0; j < nEdges; j++) {
			ed = edgeDistance(j,i);
			// if this edge is negative
			if(edgevalence[j] < 0) {
				// Add it to the negative count at the edge distance
				negCount[ed]++;
			}
			if(edgevalence[j] > 0) {
				posCount[ed]++;
			}
		}
    for(int r = 0; r < triadID.size(); r++){
        String currDir = direction[r];
        int currVal = valence[r];
        int currDist = distance[r];
        int currNode = nodeNum[r]-1;
        if(i == currNode){
          if(currDir == "OUT"){
            if(currVal < 0){
              negOutCount[currDist]++;
            }
            else{
              posOutCount[currDist]++;
            }
          }
          else{
            if(currVal < 0){
              negInCount[currDist]++;
            }
            else{
              posInCount[currDist]++;
            }
          }
        }
    }
		for(int k = 0; k <= maxDistance; k++) {
			triadicPart = piiDelta * (pow(negOutCount[k], piiX) - pow(posOutCount[k], piiX) +
                                pow(negInCount[k], piiX) - pow(posInCount[k], piiX));
			piIndex[i] += piiBetaVector[k] * (pow(posCount[k], piiX) - pow(negCount[k], piiX) + triadicPart);
			posOutCount[k] = 0;
			negOutCount[k] = 0;
      posInCount[k] = 0;
			negInCount[k] = 0;
			negCount[k] = 0;
			posCount[k] = 0;
      //Rcout << "Node: " << i << " Dist: " << k << " triadpart: " << triadicPart << std::endl;
		}
	}

	return piIndex;
}
