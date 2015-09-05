#include <stdio.h>
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
CharacterVector triadTable(IntegerMatrix edgeDistance, IntegerMatrix shortPaths, IntegerMatrix triads, CharacterVector vertices){
  CharacterVector triadID(0);
  CharacterVector direction(0);
  IntegerVector valence(0);
  NumericVector distance(0);
  IntegerVector nodeNum(0);

  //   triadID = character(), direction = character(), valence = integer(),
  //     distance = numeric(), nodeNum = integer()

  for(int i = 0; i < triads.nrow(); i++) {
    triadID.push_back(vertices[i]);
  }

  return triadID;
//   DataFrame triad_table = new DataFrame();
//   int ed;
//   for(int i = 0; i < triads.length(); i++){
//     string node1;
//     string node2;
//     string node3;
//     String triname = node1+"-"+node2+"-"+node3;
//     String edge1_2 = node1+"-"+node2;
//     String edge1_3 = node1+"-"+node3;
//     String edge2_3 = node2+"-"+node3;
//     for(int j = 0; j<vertices.length(); j++){
//       //String refnode = vertices[j];
//       //ed = edgeDistance(j,i);
//       int dis1;
//       int dis2;
//       int dis3;
//       String ref;
//       String truev;
//       String diff;
//       if (dis1 == dis2 & dis1 == dis3) {
//         int nodedis1_2;
//         int nodedis1_3;
//         int nodedis2_3;
//         if (nodedis1_2 !=  nodedis1_3) {
//           if (nodedis1_2 != nodedis2_3) {
//             ref = edge1_2;
//           }else{
//             ref = edge1_3;
//           }
//         }else{
//           ref = edge2_3;
//         }
//         if (ref == edge1_2) {
//           //truev = valence_set[1]
//         }
//         if (ref == edge1_3) {
//           //truev = valence_set[3]
//         }
//         if (ref == edge2_3) {
//           //truev = valence_set[2]
//         }
//       }
//       else{
//         if (dis1 != dis2) {
//           if (dis1 != dis3) {
//             diff = edge1_2;
//           }else{
//             diff = edge1_3;
//           }
//         }else{
//           diff = edge2_3;
//         }
//
//         if (diff == edge1_2) {
//           //truev = valence_set[1]
//         }
//         if (diff == edge1_3) {
//           //truev = valence_set[3]
//         }
//         if (diff == edge2_3) {
//           //truev = valence_set[2]
//         }
//       }
//     }
//   }
//   return;
}
