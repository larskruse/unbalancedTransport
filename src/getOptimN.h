#ifndef GET_OPTIM_N_H
#define GET_OPTIM_N_H


//' Calculating the best south-west shift and the cost savings
//'
//' This function calculates the optimal south-west shift and the cost
//' change by applying it.
//'
//' @param ukSW A numeric vector holding the u(k) values for the south-west shifts
//' @param vkSW A numeric vector holding the v(k) values for the south-west shifts
//' @param indexGroupsSWi A list of numeric vectors
//' @param indexGroupsSWj A list of numeric vectors
//' @param iList A numeric vector holing the row index of the transport path
//' @param jList A numeric vector holing the row index of the transport path
//' @param nt Reference to pass result: Index of the optimal "shift group"
//' @param k0 Reference to pass results: Position in the supply vector where mass is removed
//' @param k1 Reference to pass results: Position in the demand vector where mass is removed
//' @param shiftTransportCost Reference to pass results: Cost savings by applying the shift
//' @return nt, k0, k1 and shiftTransportCost by reference
//'
void getNtSW(Rcpp::NumericVector ukSW, Rcpp::NumericVector vkSW, Rcpp::List indexGroupsSWi, Rcpp::List indexGroupsSWj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost);





//' Calculating the best north-east shift and the cost savings
//'
//' This function calculates the optimal north-east shift and the cost
//' change by applying it.
//'
//'
//' @param ukNO A numeric vector holding the u(k) values for the north-east shifts
//' @param vkNO A numeric vector holding the v(k) values for the north-east shifts
//' @param indexGroupsNOi A list of numeric vectors
//' @param indexGroupsNOj A list of numeric vectors
//' @param iList A numeric vector holing the row index of the transport path
//' @param jList A numeric vector holing the row index of the transport path
//' @param nt Reference to pass result: Index of the optimal "shift group"
//' @param k0 Reference to pass results: Position in the supply vector where mass is removed
//' @param k1 Reference to pass results: Position in the demand vector where mass is removed
//' @param shiftTransportCost Reference to pass results: Cost savings by applying the shift
//' @return nt, k0, k1 and shiftTransportCost by reference
//'
void getNtNO(Rcpp::NumericVector ukNO, Rcpp::NumericVector vkNO, Rcpp::List indexGroupsNOi, Rcpp::List indexGroupsNOj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost);





//' Calculating the best north-east shift and the cost savings
//'
//' This function calculates the optimal north-east shift and the cost
//' change by applying it.
//'
//'
//' @param ukSW A numeric vector holding the u(k) values for the south-west shifts
//' @param vkSW A numeric vector holding the v(k) values for the south-west shifts
//' @param indexGroupsSWi A list of numeric vectors
//' @param indexGroupsSWj A list of numeric vectors
//' @param ukNO A numeric vector holding the u(k) values for the north-east shifts
//' @param vkNO A numeric vector holding the v(k) values for the north-east shifts
//' @param indexGroupsNOi A list of numeric vectors
//' @param indexGroupsNOj A list of numeric vectors
//' @param iList A numeric vector holing the row index of the transport path
//' @param jList A numeric vector holing the row index of the transport path
//' @param nt Reference to pass result: Index of the optimal "shift group"
//' @param k0 Reference to pass results: Position in the supply vector where mass is removed
//' @param k1 Reference to pass results: Position in the demand vector where mass is removed
//' @param shiftTransportCost Reference to pass results: Cost savings by applying the shift
//' @return nt, k0, k1 and shiftTransportCost by reference
//'
void getNtSWNO(Rcpp::NumericVector ukSW, Rcpp::NumericVector vkSW, Rcpp::List indexGroupsSWi, Rcpp::List indexGroupsSWj, Rcpp::NumericVector ukNO, Rcpp::NumericVector vkNO, Rcpp::List indexGroupsNOi, Rcpp::List indexGroupsNOj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost);


#endif // GET_OPTIM_N_H
