library("Rcpp")
library("RcppArmadillo")

cppFunction(code='
NumericVector arrayC(NumericVector input, IntegerVector dim) {
  input.attr("dim") = dim;
  Rcpp::Rcout << input[0,0,0];
  return input;
}
')

m <- edge.distance(g, lookup.mat = T)

cppFunction(code = '
void getCell(NumericVector inputEdgeDistance, IntegerVector Vcell) {
  IntegerVector inputEdgeDistanceDims = inputEdgeDistance.attr("dim");
  arma::cube cubeEdgeDistance(inputEdgeDistance.begin(), inputEdgeDistanceDims[0], inputEdgeDistanceDims[1], inputEdgeDistanceDims[2], false);
  Rcout << "Hello world!: " << cubeEdgeDistance(Vcell[0]-1, Vcell[1]-1, Vcell[2]-1);
}', depends = 'RcppArmadillo')

getCell(m, c(1,2,3))
#zed <- arrayC(e.dist.l, c(dimNum, dimNum, dimNum))
#zed[1,1,1]
