// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// StabilizedScaling_Rcpp
Rcpp::List StabilizedScaling_Rcpp(Eigen::Map<Eigen::MatrixXd> costMatrix, Eigen::Map<Eigen::VectorXd> supply, Eigen::Map<Eigen::VectorXd> demand, double lambdaSupply, double alphaSupply, double betaSupply, double lambdaDemand, double alphaDemand, double betaDemand, int DivSupply, int DivDemand, int iterMax, Eigen::Map<Eigen::VectorXd> epsvec, double tol);
RcppExport SEXP _unbalancedTransport_StabilizedScaling_Rcpp(SEXP costMatrixSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP lambdaSupplySEXP, SEXP alphaSupplySEXP, SEXP betaSupplySEXP, SEXP lambdaDemandSEXP, SEXP alphaDemandSEXP, SEXP betaDemandSEXP, SEXP DivSupplySEXP, SEXP DivDemandSEXP, SEXP iterMaxSEXP, SEXP epsvecSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaSupply(lambdaSupplySEXP);
    Rcpp::traits::input_parameter< double >::type alphaSupply(alphaSupplySEXP);
    Rcpp::traits::input_parameter< double >::type betaSupply(betaSupplySEXP);
    Rcpp::traits::input_parameter< double >::type lambdaDemand(lambdaDemandSEXP);
    Rcpp::traits::input_parameter< double >::type alphaDemand(alphaDemandSEXP);
    Rcpp::traits::input_parameter< double >::type betaDemand(betaDemandSEXP);
    Rcpp::traits::input_parameter< int >::type DivSupply(DivSupplySEXP);
    Rcpp::traits::input_parameter< int >::type DivDemand(DivDemandSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type epsvec(epsvecSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(StabilizedScaling_Rcpp(costMatrix, supply, demand, lambdaSupply, alphaSupply, betaSupply, lambdaDemand, alphaDemand, betaDemand, DivSupply, DivDemand, iterMax, epsvec, tol));
    return rcpp_result_gen;
END_RCPP
}
// lambertInit
Rcpp::NumericVector lambertInit(Rcpp::NumericVector& x);
RcppExport SEXP _unbalancedTransport_lambertInit(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(lambertInit(x));
    return rcpp_result_gen;
END_RCPP
}
// Sinkhorn_Rcpp
Rcpp::List Sinkhorn_Rcpp(Rcpp::NumericMatrix costMatrix, Rcpp::NumericVector& supply, Rcpp::NumericVector& demand, double lambdaSupply, double param1Supply, double param2Supply, double lambdaDemand, double param1Demand, double param2Demand, int DivSupply, int DivDemand, int iterMax, Rcpp::NumericVector& epsvec, double tol);
RcppExport SEXP _unbalancedTransport_Sinkhorn_Rcpp(SEXP costMatrixSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP lambdaSupplySEXP, SEXP param1SupplySEXP, SEXP param2SupplySEXP, SEXP lambdaDemandSEXP, SEXP param1DemandSEXP, SEXP param2DemandSEXP, SEXP DivSupplySEXP, SEXP DivDemandSEXP, SEXP iterMaxSEXP, SEXP epsvecSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaSupply(lambdaSupplySEXP);
    Rcpp::traits::input_parameter< double >::type param1Supply(param1SupplySEXP);
    Rcpp::traits::input_parameter< double >::type param2Supply(param2SupplySEXP);
    Rcpp::traits::input_parameter< double >::type lambdaDemand(lambdaDemandSEXP);
    Rcpp::traits::input_parameter< double >::type param1Demand(param1DemandSEXP);
    Rcpp::traits::input_parameter< double >::type param2Demand(param2DemandSEXP);
    Rcpp::traits::input_parameter< int >::type DivSupply(DivSupplySEXP);
    Rcpp::traits::input_parameter< int >::type DivDemand(DivDemandSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type epsvec(epsvecSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Sinkhorn_Rcpp(costMatrix, supply, demand, lambdaSupply, param1Supply, param2Supply, lambdaDemand, param1Demand, param2Demand, DivSupply, DivDemand, iterMax, epsvec, tol));
    return rcpp_result_gen;
END_RCPP
}
// Hausdorff_Vec_Rcpp
Rcpp::NumericVector Hausdorff_Vec_Rcpp(Rcpp::NumericMatrix costMatrix, Rcpp::NumericVector& distribution, Rcpp::NumericVector& f, double lambda, double param1, double param2, int Div, double eps);
RcppExport SEXP _unbalancedTransport_Hausdorff_Vec_Rcpp(SEXP costMatrixSEXP, SEXP distributionSEXP, SEXP fSEXP, SEXP lambdaSEXP, SEXP param1SEXP, SEXP param2SEXP, SEXP DivSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type param1(param1SEXP);
    Rcpp::traits::input_parameter< double >::type param2(param2SEXP);
    Rcpp::traits::input_parameter< int >::type Div(DivSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(Hausdorff_Vec_Rcpp(costMatrix, distribution, f, lambda, param1, param2, Div, eps));
    return rcpp_result_gen;
END_RCPP
}
// treegkr_Rcpp
Rcpp::List treegkr_Rcpp(Rcpp::List& tree, Rcpp::NumericVector& supply, Rcpp::NumericVector& demand, Rcpp::NumericVector& creation, Rcpp::NumericVector& destruction);
RcppExport SEXP _unbalancedTransport_treegkr_Rcpp(SEXP treeSEXP, SEXP supplySEXP, SEXP demandSEXP, SEXP creationSEXP, SEXP destructionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type supply(supplySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type demand(demandSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type creation(creationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type destruction(destructionSEXP);
    rcpp_result_gen = Rcpp::wrap(treegkr_Rcpp(tree, supply, demand, creation, destruction));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_unbalancedTransport_StabilizedScaling_Rcpp", (DL_FUNC) &_unbalancedTransport_StabilizedScaling_Rcpp, 14},
    {"_unbalancedTransport_lambertInit", (DL_FUNC) &_unbalancedTransport_lambertInit, 1},
    {"_unbalancedTransport_Sinkhorn_Rcpp", (DL_FUNC) &_unbalancedTransport_Sinkhorn_Rcpp, 14},
    {"_unbalancedTransport_Hausdorff_Vec_Rcpp", (DL_FUNC) &_unbalancedTransport_Hausdorff_Vec_Rcpp, 8},
    {"_unbalancedTransport_treegkr_Rcpp", (DL_FUNC) &_unbalancedTransport_treegkr_Rcpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_unbalancedTransport(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
