#ifndef DIVERGENCES_H
#define DIVERGENCES_H

#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::VectorXd div0(Eigen::VectorXd a, Eigen::VectorXd b);
    
Eigen::VectorXd axb0(Eigen::VectorXd a, Eigen::VectorXd b);


double vectorDivergence (Eigen::VectorXd r, Eigen::VectorXd s, int DivFun,
                         double param1, double param2 = 0, double param3 = 0);
    
        
double fVectorDivergence (Eigen::VectorXd p, Eigen::VectorXd u, int DivFun,
                          double param1, double param2 = 0, double param3 = 0);
            
double dualSolSummand(Eigen::VectorXd u, Eigen::VectorXd v, double eps,
                      Eigen::MatrixXd Kernel);
                
            
            
Eigen::MatrixXd updateK(Eigen::VectorXd u, Eigen::VectorXd v, double eps,
                        Eigen::MatrixXd costMatrix);
                
double dualSolSummandSink(Eigen::VectorXd u, Eigen::VectorXd v, double eps,
                        Eigen::MatrixXd costMatrix, Eigen::MatrixXd supdem);
            
#endif
