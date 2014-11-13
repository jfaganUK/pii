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
NumericVector piiCalc(IntegerMatrix edgeDistance, NumericVector valence, double piiBeta, double piiX, int maxDistance) {
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
			if(valence[j] < 0) {
				negCount[edgeDistance(j,i)]++;
			}
			if(valence[j] > 0) {
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

// [[Rcpp::export]]
NumericVector piiTriadicCalc(IntegerMatrix edgeDistance, NumericVector valence, double piiBeta, double piiX, int maxDistance, LogicalVector edgeTriadic, double piiDelta) {
	int nEdges = edgeDistance.nrow(), nNodes = edgeDistance.ncol();
	IntegerVector negCount(maxDistance+1), posCount(maxDistance+1);
	IntegerVector negTriadCount(maxDistance+1), posTriadCount(maxDistance+1);
	NumericVector piiBetaVector(maxDistance + 1);
	NumericVector piIndex(nNodes);
	double triadicPart;
	int ed; // a holder for the current edge distance

	// Initialize piiBetaVector
	for(int k = 0; k <= maxDistance; k++) {
		piiBetaVector[k] = pow(piiBeta, k);
	}

	for(int i = 0; i < nNodes; i++) {
		for(int j = 0; j < nEdges; j++) {
			ed = edgeDistance(j,i);

			// if this edge is negative
			if(valence[j] < 0) {
				// Add it to the negative count at the edge distance
				negCount[ed]++;

				// If the edge is a member of a triad
				if(edgeTriadic[j]) {
					negTriadCount[ed]++;
				}
			}
			if(valence[j] > 0) {
				posCount[ed]++;
				if(edgeTriadic[j]) {
					posTriadCount[ed]++;
				}
			}
		}
		for(int k = 0; k <= maxDistance; k++) {
			triadicPart = piiDelta * (pow(posTriadCount[k], piiX) - pow(negTriadCount[k], piiX));
			piIndex[i] += piiBetaVector[k] * (pow(posCount[k], piiX) - pow(negCount[k], piiX) + triadicPart);
			posTriadCount[k] = 0;
			negTriadCount[k] = 0;
			negCount[k] = 0;
			posCount[k] = 0;
		}
	}

	return piIndex;
}