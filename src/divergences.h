#ifndef DIVERGENCES_H
#define DIVERGENCES_H

#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::VectorXd div0(const Eigen::VectorXd a, const Eigen::VectorXd b);
    
Eigen::VectorXd axb0(const Eigen::VectorXd a, const Eigen::VectorXd b);


double vectorDivergence (const Eigen::VectorXd r, const Eigen::VectorXd s, 
                         const int DivFun, const double param1, 
                         const double param2 = 0, const double param3 = 0);
        
double fVectorDivergence (const Eigen::VectorXd p, const Eigen::VectorXd u,
                          const int DivFun, const double param1, 
                          const double param2 = 0, const double param3 = 0);
            
double dualSolSummand(const Eigen::VectorXd u, const Eigen::VectorXd v, 
                      const double eps, const Eigen::MatrixXd Kernel);

Eigen::MatrixXd updateK(const Eigen::VectorXd u, const Eigen::VectorXd v, 
                        const double eps, const Eigen::MatrixXd costMatrix);
                
double dualSolSummandSink(const Eigen::VectorXd u, const Eigen::VectorXd v, 
                          const double eps, const Eigen::MatrixXd costMatrix,
                          const Eigen::MatrixXd supdem);
            
#endif
