#ifndef UTILMONGE_H
#define UTILMONGE_H

//' C++ implantation of the North-West Corner Rule
//'
//' This function calculates the optimal transport cost for a balanced optimal
//' transport problem if its cost matrix fulfills the Monge property. Otherwise,
//' it calculates a feasible solution for the transport problme.
//'
//' @param costMatrix A numeric matrix.
//' @param supply A numeric supply vector.
//' @param demand A numeric demand vector.
//' @param iList A numeric vector holding the row index of the transport path
//' @param jList A numeric vector holding the row index of the transport path
//' @param weightList A numeric vector holding the transport path weights
//' @return The transport cost.
//' @noRd

double Nw_Corner_Rule(Rcpp::NumericMatrix costMatrix, Rcpp::NumericVector supply, Rcpp::NumericVector demand, Rcpp::NumericVector &iList, Rcpp::NumericVector &jList, Rcpp::NumericVector &weightList);


//' Expanding an unbalanced to a balanced optimal transport problem
//'
//' This function expands a given unbalanced optimal transport problem to a balance
//' transport problem by adding mass to the supply or demand vector and expanding
//' the cost matrix. The Monge property of the cost matrix is maintained.
//' The resulting transport problem can be solved using the North-West Corner Rule.
//'
//' @param costMatrix A numeric matrix.
//' @param supply A numeric supply vector.
//' @param demand A numeric demand vector.
//' @return None.
//' @noRd
void UnbalancedToBalanced(Rcpp::NumericMatrix &costMatrix, Rcpp::NumericVector &p, Rcpp::NumericVector &q);





//' Calculating which shifts are possible and if they are in the north-east
//' or south-west direction
//'
//' @param iList A numeric vector holding the row index of the transport path
//' @param jList A numeric vector holding the row index of the transport path
//' @param indexGroupsSWi A numeric vector.
//' @param indexGroupsSWj A numeric vector.
//' @param indexGroupsSWi A numeric vector.
//' @param indexGroupsSWj A numeric vector.
//' @noRd
void GetIndexGroups(Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::List &indexGroupsSWi, Rcpp::List &indexGroupsSWj, Rcpp::List &indexGroupsNOi, Rcpp::List &indexGroupsNOj);


#endif // UTILMONGE_H
