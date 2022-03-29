#include <algorithm>
#include <RcppEigen.h>
#include <math.h>
#include <Rcpp.h>
#include "divergences.h"


//' The proxdiv operator
//'
//' The proxdiv operators for different divergences in the form of 'lambda * DivFun(.|p)'
//' Implemented are the operators for the Kullback-Leibler divergence and total variation.
//'
//' @param lambda Regularization parameter
//' @param p A numeric vector
//' @param s A numeric vector
//' @param u A numeric vector
//' @param eps The epsilon value
//' @param DivFun A numeric value indicating the function to be used.'1' gives
//'   the proxdiv operator for the Kullback-Leibler divergence and '2' the opterator
//'   for the total variation.
//' @param alpha num value
//' @param beta num value
//' @return A vector holding the proxdiv evaluation
//' @noRd
Eigen::VectorXd proxdiv(const double lambda,
                        const Eigen::VectorXd p,
                        const Eigen::VectorXd s,
                        const Eigen::VectorXd u,
                        const double eps,
                        const int DivFun,
                        const double alpha,
                        const double beta){
  
  Eigen::VectorXd temp;
    
    // Kullback-Leibler
  if (DivFun == 1){
      
    temp = s.array()*exp(u.array()/lambda);
    temp = div0(p,temp);
    temp = temp.array().pow(lambda/(lambda+eps));
    return temp;

    //Total variation
  }else if(DivFun == 2){
      
      temp = ((lambda-u.array())/eps).array().exp().array().min(div0(p,s).array().max((-(lambda+u.array())/eps).array().exp()));
      
      for (int i =0 ; i < u.size(); i++){
          
          if(isinf(u(i))){
              
                temp(i) = 0;
              }
          
          }
      
    return (temp);
      // Range Constraint 
  }else{
      
      temp= (beta*div0(p,s)).array().min((alpha*div0(p,s)).array().max((-u.array()/eps).array().exp()));
      
      for (int i =0 ; i < u.size(); i++){
          
          if(isinf(u(i))){
              
              temp(i) = 0;
          }
          
      }
      
    return(temp);
  }

}



//' The stabilized Scaling Algorithm
//'
//' C++ implementation of the log-domain stabilized Version of the Scaling
//' Algorithm.
//'
//' @param costMatrix A numeric matrix.
//' @param supply A numeric vector
//' @param demand A numeric vector
//' @param lambdaSupply Parameter for the supply proxdiv function
//' @param lambdaDemand Parameter for the demand proxdiv function
//' @param DivSupply Parameter indicating the divergence function to use for the supply proxdiv function
//' @param DivDemand Parameter indicating the divergence function to use for the demand proxdiv function
//' @param iterMax Maximum number of iterations
//' @param epsvec A numeric vector of decreasing epsilon values.
//' @param alphaSupply numeric Value
//' @param betaSupply numeric value
//' @param alphaDemand numeric Value
//' @param betaDemand numeric value
//' @param tol num vale
//'
//' @return The optimal transport plan
//' @noRd
//[[Rcpp::export]]
Rcpp::List StabilizedScaling_Rcpp(const Eigen::Map<Eigen::MatrixXd> costMatrix,
                                  const Eigen::Map<Eigen::VectorXd> supply,
                                  const Eigen::Map<Eigen::VectorXd> demand,
                                  const double lambdaSupply,
                                  const double alphaSupply,
                                  const double betaSupply,
                                  const double lambdaDemand,
                                  const double alphaDemand,
                                  const double betaDemand,
                                  const int DivSupply,
                                  const int DivDemand,
                                  const size_t iterMax,
                                  const Eigen::Map<Eigen::VectorXd> epsvec,
                                  const double tol = 1e-7){
   
   
   
   

    // number of absorptions
    size_t numAbs {0};

    // number of points in the reference measures
    size_t Nx = supply.size();
    size_t Ny = demand.size();
    
    
    
    // initializing vectors
    Eigen::VectorXd a = Eigen::VectorXd::Ones(Nx);
    Eigen::VectorXd b = Eigen::VectorXd::Ones(Ny);
    
    // stabilization vectors
    Eigen::VectorXd u = Eigen::VectorXd::Zero(Nx);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(Ny);
    
    Eigen::VectorXd uPrev = Eigen::VectorXd::Zero(Nx);
    

    Eigen::VectorXd u0 = Eigen::VectorXd::Zero(Nx);
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(Ny);
    
    
    // main loop iteration counter
    size_t i {1};
    
    // epsilon value index counter
    size_t epsind {0};
    
    // setting first epsilon value
    double eps = epsvec(0);
    
    
    // computing the initial kernel
    // since u and v are 0, the updateK function returns the Gibbs kernel
    Eigen::MatrixXd Kernel = updateK(u, v, eps, costMatrix);
    Eigen::MatrixXd gaussKernel = Kernel;

    // primal and dual cost variables
    double pCost {0};
    double dCost {0};
    
    // if true: decrease epsilon
    bool incEps = false;
    
    // Vector representatin of the kernel and gauss kernel
    Eigen::VectorXd KVec;
    Eigen::VectorXd gKVec(Eigen::Map<Eigen::VectorXd>(gaussKernel.data(), gaussKernel.cols()*gaussKernel.rows()));  

    while(i < iterMax ){
        // update iteration counter
        i++;
        
        // calculate scaling iterates
        a = Kernel * b;
        a = proxdiv(lambdaSupply, supply, a, u, eps, DivSupply, alphaSupply, betaSupply);
    
    
        b = Kernel.transpose() * a;
        b = proxdiv(lambdaDemand, demand, b, v, eps, DivDemand, alphaDemand, betaDemand);
    
        //Stabilizing step and changing epsilon
        // called when:
        //  1. a or b are too large,
        //  2. a new value for epsilon has to be assigned
        //  3. in the last iteration to calculate the transport map
        //  4. to calculate the difference between ´
        if ((abs(a.array()) > 1e+100).any() ||
            (abs(b.array()) > 1e+100).any() ||
            (static_cast<double>(i)/static_cast<double>(iterMax)) >
            static_cast<double>(epsind + 1)/static_cast<double>(epsvec.size()) ||
            i == iterMax - 1  || i % 50 == 0 || incEps){
                
            uPrev = u;
                
            // absorbing a/b in u/v
            u = u.array() + eps*(a.array().log());
            v = v.array() + eps*(b.array().log());


            //updating epsilon
            if((static_cast<double>(i)/static_cast<double>(iterMax)) >
               static_cast<double>(epsind + 1)/static_cast<double>(epsvec.size())||
               incEps){
                epsind++;
                eps = epsvec(epsind);

                gaussKernel = updateK(Eigen::VectorXd::Zero(Nx),
                                      Eigen::VectorXd::Zero(Ny),
                                      eps,
                                      costMatrix);
                gKVec = Eigen::Map<Eigen::VectorXd>(gaussKernel.data(),
                                                    gaussKernel.cols()*gaussKernel.rows());

                incEps = false;
            }
              // update Kernel according to u and v
            Kernel = updateK(u, v, eps, costMatrix);
            KVec = Eigen::Map<Eigen::VectorXd>(Kernel.data(),
                                               Kernel.cols()*Kernel.rows());
            
            //reset a and b
            a = Eigen::VectorXd::Ones(Nx);
            b = Eigen::VectorXd::Ones(Ny);

            // check if stopping criterion is reached
            if(abs((u.array()-uPrev.array()).maxCoeff()) < tol){
            
                if(epsind == epsvec.size()-1){

                    return Rcpp::List::create(Rcpp::Named("TransportPlan") = Kernel,
                                                  Rcpp::Named("dual_f") = u,
                                                  Rcpp::Named("dual_g") = v,
                                                  Rcpp::Named("pCost") = pCost,
                                                  Rcpp::Named("dCost") = dCost,
                                                  Rcpp::Named("iterations") = i);
                }else{
                        incEps = true;
                }
            
            } 
            
        }
    
    }
  
        
    for(size_t i = 0; i < Nx; i ++){
            
        for(size_t j= 0; j < Ny; j++){
        
            Kernel(i,j) = Kernel(i,j)*a(i)*b(j);  

        }
      
    }
  
  
    // calculate primal cost
    u0 = u;
    v0 = v;
    KVec = Eigen::Map<Eigen::VectorXd>(Kernel.data(),
                                       Kernel.cols()*Kernel.rows());
    
    pCost =  vectorDivergence(KVec,
                              gKVec,
                              1,
                              eps);
    
    pCost += vectorDivergence(Kernel.rowwise().sum(),
                              supply,
                              DivSupply,
                              lambdaSupply,
                              alphaSupply,
                              betaSupply);
  
    pCost += vectorDivergence(Kernel.colwise().sum(),
                              demand,
                              DivDemand,
                              lambdaDemand,
                              alphaDemand,
                              betaDemand);
  

    // calculate dual cost
    dCost = -dualSolSummand(u,v,eps,gaussKernel);
    for(size_t i = 0; i < u.size(); i++){
          if(supply(i) == 0){
              u0(i) = 0;
        }
    }
    for(size_t i = 0; i < v.size(); i++){
        if(demand(i) == 0){
            v0(i) = 0;
        }
    }

  
    dCost -= fVectorDivergence(supply,
                               -u0,
                               DivSupply,
                               lambdaSupply,
                               alphaSupply,
                               betaSupply);
  
  

    dCost -= fVectorDivergence(demand,
                               -v0,
                               DivDemand,
                               lambdaDemand,
                               alphaDemand,
                               betaDemand);
  

 
  // returning the transport plan
  // since the absorption is called in the last iteration of the loop,
  // the transport plan is equal to the kernel.
    return Rcpp::List::create(Rcpp::Named("TransportPlan") = Kernel,
                          Rcpp::Named("dual_f") = u,
                          Rcpp::Named("dual_g") = v,
                          Rcpp::Named("pCost") = pCost,
                          Rcpp::Named("dCost") = dCost,
                          Rcpp::Named("iterations") = i);

}
