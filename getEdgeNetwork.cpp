#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix getEdgeNetworkCalc(IntegerMatrix inputEdgeList) {
	int nr = inputEdgeList.nrow();
	int maxNodeID = 0;
	IntegerVector edgeNames(nr);
	IntegerMatrix outputEdgeList(nr * 2, 2);

	// Get the maximum node id
	for(int i = 0; i < nr; i++) {
		if(inputEdgeList(i,1) > maxNodeID) {
			maxNodeID = inputEdgeList(i,0);
		} 
		if(inputEdgeList(i,2) > maxNodeID) {
			maxNodeID = inputEdgeList(i,1);
		}
	}

	for(int i = 1; i <= nr; i++) {
		edgeNames(i-1) = maxNodeID + i;
	}


	for(int i = 0; i < nr; i++) {
		// from v1 to the edge
		outputEdgeList(i * 2, 0) = inputEdgeList(i, 0);
		outputEdgeList(i * 2, 1) = edgeNames(i);

		// from the edge to v2
		outputEdgeList(i * 2 + 1, 0) = edgeNames(i);
		outputEdgeList(i * 2 + 1, 1) = inputEdgeList(i, 1);
	}

	return outputEdgeList;
}


