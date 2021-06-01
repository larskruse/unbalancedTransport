#ifndef SHIFTS_H
#define SHIFTS_H



//' Calculating u(k) and v(k) for the south-west-shifts
//'
//'
//' @param costMatrix A numeric matrix that fulfills the Monge property
//' @param iList A numeric vector holing the row index of the transport path
//' @param jList A numeric vector holing the row index of the transport path.
//' @param ukSW Reference to return the u(k) vector
//' @param vkSW Reference to return the v(k) vector
//' @return The vectors u(k) and v(k) for the south-west shifts
//' @noRd
//'
void CalcSWShift(Rcpp::NumericVector costMatrix, Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::NumericVector &ukSW, Rcpp::NumericVector &vkSW);



//' Calculating u(k) and v(k) for the north-east shifts
//'
//'
//' @param costMatrix A numeric matrix that fulfills the Monge property
//' @param iList A numeric vector holing the row index of the transport path
//' @param jList A numeric vector holing the row index of the transport path.
//' @param ukSW Reference to return the u(k) vector
//' @param vkSW Reference to return the v(k) vector
//' @return The vectors u(k) and v(k) for the north-east shifts
//' @noRd
//'
void CalcNOShift(Rcpp::NumericVector costMatrix, Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::NumericVector &ukNO, Rcpp::NumericVector &vkNO);


#endif // SHIFTS_H
