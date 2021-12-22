// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include "divergences.h"


// Worker for RcppParallel 
struct Lse : public RcppParallel::Worker
{
    
    const RcppParallel::RVector<double> logMeasure;
    const RcppParallel::RVector<double> dualPot;
    const RcppParallel::RMatrix<double> costMatrix;
    const double eps;
    
        
    RcppParallel::RVector<double> lseVec;
    
    tthread::mutex mutexDualPot_;
    tthread::mutex mutexDualMeasure_;
    
    // constructors
    Lse(
        const Rcpp::NumericVector logMeasure,
        const Rcpp::NumericVector dualPot,
        const Rcpp::NumericMatrix costMatrix,
        const double eps,
        Rcpp::NumericVector lseVec) :
        logMeasure(logMeasure),
        dualPot(dualPot), costMatrix(costMatrix), eps(eps),
        lseVec(lseVec) {}
    
    
    void operator()(std::size_t begin, std::size_t end) {
        
        for(std::size_t i = begin; i < end; i++){
            
            // the lse is calcualted for each colum of the cost matrix
            RcppParallel::RMatrix<double>::Column col = costMatrix.column(i);
            std::vector<double> resVec(col.length());
            
            
            mutexDualPot_.lock();
            std::transform(col.begin(), col.end(), dualPot.begin(), resVec.begin(),
                           [](double i, double j){return (j - i);});
            mutexDualPot_.unlock();

            
            mutexDualMeasure_.lock();
            std::transform(resVec.begin(), resVec.end(), logMeasure.begin(), resVec.begin(),
                           [&](double i, double j){return (j + i/eps) ;});
            mutexDualMeasure_.unlock();
            

            double maxVal = *std::max_element(resVec.begin(), resVec.end());
            std::for_each(resVec.begin(), resVec.end(),
                          [&](double &i) {i = std::exp(i - maxVal);});
            

            
            lseVec[i] =
                std::log(std::accumulate(resVec.begin(), resVec.end(),0.0)) +  maxVal;
        
        }
        
    }
    
};






//' initial value for the lambert W function
//' @param x numeric vector
//' @return initial value 
//' @noRd
Rcpp::NumericVector lambertInit(Rcpp::NumericVector& x){
    Rcpp::NumericVector temp;
    Rcpp::NumericVector z(x.length());
    z[x>1] = x[x>1];
    
    temp = x[x<-2];
    temp = Rcpp::exp(temp)*(1-Rcpp::exp(temp));
    z[x<-2] = temp;
    
    temp = Rcpp::exp(x[x <= 1 & x>=-2]);
    temp = temp*(3+6*temp+temp*temp)/(3+9*temp+5*temp*temp);
    z[x<=1 & x>= -2] = temp;
    z[z == 0] = 1e-6;
    
    return(z);
}



//' computing the lambert w function
//' @param x numeric vector
//' @return the function value
//' @noRd 
// 
Rcpp::NumericVector lambertWFunctionLog(Rcpp::NumericVector& x){
    
    
    Rcpp::NumericVector z = lambertInit(x);
    
    Rcpp::NumericVector c;
    Rcpp::NumericVector b;
    
    double eps = std::numeric_limits<double>::epsilon();
    
    for(int i = 0; i < 4; i++){
        
        c = z*(Rcpp::log(z+eps)+z-x)/(1+z);
        b = -1/(z*(1+z));
        z = Rcpp::pmax(z-c/(1-0.5*c*b),eps);
    }

    return(z);
}



//' computing the inital values for the iterators
//' @param x numeric vector
//' @return the function value
//' @noRd
Rcpp::NumericVector inital(Rcpp::NumericVector& x){
    
    Rcpp::NumericVector z = lambertInit(x);
    Rcpp::NumericVector c;
    Rcpp::NumericVector b;
    
    double eps = std::numeric_limits<double>::epsilon();
    
    for(int i = 0; i < 4; i++){
        
        c = z*(Rcpp::log(z+eps)+z-x)/(1+z);
        b = -1/(z*(1+z));
        z = Rcpp::pmax(z-c/(1-0.5*c*b),eps);
        
    }
    
    return(z);
}




//' The aprox operator
//'
//' The aprix operators for different divergences in the form of 'lambda * DivFun(.|p)'
//' Implemented are the operators for the Kullback-Leibler divergence and total variation.
//'
//' @param lambda Regularization parameter
//' @param p A numeric vector
//' @param eps The epsilon value
//' @param DivFun A numeric value indicating the function to be used.'1' gives
//'   the proxdiv operator for the Kullback-Leibler divergence and '2' the opterator
//'   for the total variation.
//' @param param1 num value
//' @param param2 num value
//' @return A vector holding the proxdiv evaluation
//' @noRd
Rcpp::NumericVector aprox(double lambda,
                          Rcpp::NumericVector& p,
                          double eps,
                          int DivFun,
                          double param1,
                          double param2){
    
    Rcpp::NumericVector temp(p.length());
    
    if (DivFun == 1){
        
        temp = (lambda/(lambda+eps)) * p;
        return temp;
    
    }else if(DivFun == 2){
        
        temp = Rcpp::pmax(p,-lambda);
        temp = Rcpp::pmin(temp,lambda);
        return(temp);
        
    }else if(DivFun == 3){
        
        temp = Rcpp::pmax(0, p-(eps*std::log(param2)));
        temp = Rcpp::pmin(temp, p-(eps*std::log(param1)));
        
        return(temp);
        
    // Power
    }else if(DivFun == 4){
        
        temp = -p/(eps*(1-param1));
        temp = temp + lambda/eps;
        temp = temp + std::log(lambda/eps);
        temp = (1-param1)*(lambda-eps*lambertWFunctionLog(temp));
        
        return(temp);
        
    }else{
        
        return(p);
        
    }
}

//' matrix mulitplication
//' @param mat num mat
//' @param vec num vec
//' @param Nx num val
//' @param Ny num val
//' @noRd
Rcpp::NumericVector matMul(Rcpp::NumericMatrix& mat,
                           Rcpp::NumericVector& vec,
                           int Nx, int Ny){
    
    Rcpp::NumericVector res(Nx);
    
    for(size_t i = 0; i < Nx; i++){
        for(size_t j = 0; j < Ny; j++){
       
            res(i) = res(i) + vec(j)*mat(i,j);
            
        }    
        
    }
    
    return(res);
    
}








//' Inital values
//'
//' The aprix operators for different divergences in the form of 'lambda * DivFun(.|p)'
//' Implemented are the operators for the Kullback-Leibler divergence and total variation.
//'
//' @param lambda Regularization parameter
//' @param costMatrix A numeric matrix
//' @param distribution A numeric vector
//' @param secDistribution A numeric vector
//' @param DivFun A numeric value indicating the function to be used.'1' gives
//'   the proxdiv operator for the Kullback-Leibler divergence and '2' the opterator
//'   for the total variation.
//' @param param1 num value
//' @param param2 num value
//' @param Nx num value
//' @param Ny num value
//' @param eps eps value
//' @return A vector holding the proxdiv evaluation
//' @noRd
Rcpp::NumericVector init_vectors(double lambda,
                                 Rcpp::NumericMatrix costMatrix,
                                 Rcpp::NumericVector& distribution, 
                                 Rcpp::NumericVector& secDistribution,
                                 int DivFun,
                                 double param1,
                                 double param2,
                                 int Nx,
                                 int Ny,
                                 double eps){
    
    
    Rcpp::NumericVector temp(Nx);

    double val;
    
    
   
    
    //KL
    if (DivFun == 1){
        
        // val = Rcpp::sum(distribution);
        
        temp.fill(-lambda*std::log(Rcpp::sum(distribution)));

        return temp;

    //TV
    }else if(DivFun == 2){
        
        if(std::log(Rcpp::sum(distribution)) < 0){
            temp.fill(lambda);
        }else if(std::log(Rcpp::sum(distribution)) > 0){
            temp.fill(-lambda);
        }else{
            
            temp = matMul(costMatrix, secDistribution, Nx, Ny);
 
            temp = temp + 0.5*distribution*temp;
            
            temp = aprox(lambda, temp, eps, DivFun, param1, param2);
            
        }
        
        return(temp);
        
    //Range
    }else if(DivFun == 3){
        return(temp);
        
        
    //Power
    }else if(DivFun == 4){
        temp.fill(lambda*(1-param1)*(std::pow(Rcpp::sum(distribution),(1/(param1-1)))-1));
        
        return(temp);
    }else{
        return(temp);
    }
    
    
    
}

//' The LogSumExp Function
//'
//' @param vec a numeric vector
//' @return The function LogSumExp applied to vec
//' @noRd
double logsumexp(Rcpp::NumericVector vec){
    
    if(vec.length() >  0){
        
        
        double maxVal = Rcpp::max(vec);
        size_t len = vec.length();
        
        double lse = 0;
        
        for(int i = 0; i < len; i++){
        
            lse += exp(vec(i) - maxVal);    
            
        }
        
        return(log(lse) + maxVal);
            
    
    }else{
    
        return(0);    
        
    }
    
}




//' Log-Sum-Exp Operator
//'
//' A parallel implementation of the LSE operator used in the iterators of the
//' Sinkhorn algorithm.
//'
//' @param logSup The log values of the input measures
//' @param f the current value of the dual iterator
//' @param costMatrix the cost matrix
//' @param eps Regularization Parameter
//' @param Ny Number of elements in the result vector
//' @noRd
Rcpp::NumericVector parallelVectorLse(Rcpp::NumericVector logSup,
                                      Rcpp::NumericVector f,
                                      Rcpp::NumericMatrix costMatrix,
                                      double eps,
                                      int Ny) {
    
    Rcpp::NumericVector ret(Ny);
    
    // declare the Lse instance 
    Lse lse(logSup, f, costMatrix, eps, ret);
    
    // call parallel_For to start the parallel computation.
    parallelFor(0, costMatrix.ncol(), lse);
    
    ret = lse.lseVec;
    
    // return the computed lse vector
    return eps * ret;
}





//' The Sinkhorn Algorithm
//'
//' C++ implementation of the Sinkhorn Algorithm.
//'
//' @param costMatrix A numeric matrix.
//' @param supply A numeric vector
//' @param demand A numeric vector
//' @param lambdaSupply Parameter for the supply proxdiv function
//' @param lambdaDemand Parameter for the demand proxdiv function
//' @param DivSupply Parameter indicating the divergence function to use for the supply proxdiv function
//' @param DivDemand Parameter indicating the divergence function to use for the demand proxdiv function
//' @param iterMax Maximum number of iterations
//' @param eps A numeric vector of decreasing epsilon values.
//' @param param1Supply numeric Value
//' @param param2Supply numeric value
//' @param param1Demand numeric Value
//' @param param2Demand numeric value
//' @param tol numeric value
//'
//' @return The optimal transport plan
//' @noRd
//[[Rcpp::export]]
Rcpp::List Sinkhorn_Rcpp(Rcpp::NumericMatrix &costMatrix,
                                  Rcpp::NumericVector& supply,
                                  Rcpp::NumericVector& demand,
                                  double lambdaSupply,
                                  double param1Supply,
                                  double param2Supply,
                                  double lambdaDemand,
                                  double param1Demand,
                                  double param2Demand,
                                  int DivSupply,
                                  int DivDemand,
                                  int iterMax,
                                  Rcpp::NumericVector& epsvec,
                                  double tol,
                                  Eigen::Map<Eigen::MatrixXd> supdem){
                                  
                                 
    // set first epsilon
    int epsind = 0;
    double eps = epsvec(0);
    
    
   
    // number of points in the reference measures
    int Nx = supply.length();
    int Ny = demand.length();
    
   
    
    Rcpp::NumericVector logSup = Rcpp::log(supply);
    Rcpp::NumericVector logDem = Rcpp::log(demand);

    // initializing vectors
    Rcpp::NumericVector f = init_vectors(lambdaDemand,
                                         Rcpp::transpose(costMatrix),
                                         demand,
                                         supply,
                                         DivDemand,
                                         param1Demand,
                                         param2Demand,
                                         Ny, Nx, eps);

    Rcpp::NumericVector g = init_vectors(lambdaSupply,
                                             costMatrix,
                                             supply,
                                             demand,
                                             DivSupply,
                                             param1Supply,
                                             param2Supply,
                                             Nx, Ny, eps); 
    
    
    
    Rcpp::NumericVector f_prev(Nx);
    
    //Transport map
    Rcpp::NumericMatrix Kernel(Nx,Ny);
    
    // indicates of eps has to be decreased
    bool incEps = false;
    
    
    size_t k = 0;
    
    while(k < iterMax){
        
        f_prev = Rcpp::clone(f);
        
        // calculating the iterators
        g = parallelVectorLse(logSup, f, costMatrix, eps, Ny);
        g = -aprox(lambdaDemand, g, eps, DivDemand, param1Demand, param2Demand);
        
        f = parallelVectorLse(logDem, g, Rcpp::transpose(costMatrix), eps , Nx);
        f = -aprox(lambdaSupply, f, eps, DivSupply, param1Supply, param2Supply);

        
        // check if difference of iterators is smaller than the given tolerance
        if(Rcpp::max(Rcpp::abs(f-f_prev)) < tol ){
            
            if(epsind == epsvec.length()-1){
                
                break;
                
            }else{
                
                incEps = true;
                
            }
            
        }

        // decrease epsilon 
        if((static_cast<double>(k)/static_cast<double>(iterMax)) >
               static_cast<double>(epsind + 1)/static_cast<double>(epsvec.length()) || 
           incEps){
            
            epsind = epsind + 1;
            eps = epsvec(epsind);
            
            incEps = false;
            
        }
        
      k++;
        
    }


    
    // compute transport plan
    for(int i = 0; i < Nx; i++){
        
        for(int j = 0; j < Ny ; j++){
            
            Kernel(i,j) = exp((f(i) + g(j) - costMatrix(i,j))/eps)*supdem(i,j);
        }
        
    }

    
    double pCost = 0;
    double dCost = 0;
    
    
    Eigen::Map<Eigen::VectorXd> f0(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(f));
    Eigen::Map<Eigen::VectorXd> g0(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(g));
    
    // vectorized versions of the kernels
    Eigen::VectorXd gKVec;
    Eigen::VectorXd KVec;


    Eigen::Map<Eigen::MatrixXd> EKernel(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(Kernel));
    KVec = Eigen::Map<Eigen::VectorXd>(EKernel.data(), EKernel.cols()*EKernel.rows());


 
    Eigen::MatrixXd gaussKernel;

    gaussKernel = updateK(Eigen::VectorXd::Zero(Nx),
                          Eigen::VectorXd::Zero(Ny),
                          eps,
                          Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(costMatrix));
    
    gKVec = Eigen::Map<Eigen::VectorXd>(gaussKernel.data(),
                                        gaussKernel.cols()*gaussKernel.rows());

    
    
    // calcualte primal cost
    pCost =  vectorDivergence(KVec, gKVec, 1, eps);
    
    pCost += vectorDivergence(EKernel.rowwise().sum(),
                              Rcpp::as<Eigen::Map<Eigen::VectorXd> >(supply),
                              DivSupply,
                              lambdaSupply,
                              param1Supply,
                              param2Supply);

    
    
    pCost += vectorDivergence(EKernel.colwise().sum(),
                              Rcpp::as<Eigen::Map<Eigen::VectorXd> >(demand),
                              DivDemand,
                              lambdaDemand,
                              param1Demand,
                              param2Demand);
    // calcualte dual cost
    dCost  = - dualSolSummandSink(f0,
                                  g0,
                                  eps,
                                  Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(costMatrix),
                                  supdem);



    dCost -= fVectorDivergence(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(supply),
                               -f0,
                               DivSupply,
                               lambdaSupply,
                               param1Supply,
                               param2Supply);

    dCost -= fVectorDivergence(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(demand),
    -g0, DivDemand, lambdaDemand,  param1Demand, param2Demand);
    

    
    double conv = Rcpp::max(Rcpp::abs(f-f_prev));
    
    // returnING the transport plan
    // since the absorbtion is called in the last iteration of the loop,
    // the transport plan is equal to the kernel.
    return Rcpp::List::create(Rcpp::Named("TransportPlan") = Kernel,
                              Rcpp::Named("cost") = pCost,
                              Rcpp::Named("dualCost") = dCost,
                              Rcpp::Named("dual_f") = f,
                              Rcpp::Named("dual_g") = g,
                              Rcpp::Named("converge") = conv,
                              Rcpp::Named("Iterations") = k);

}



//' Helper Function Hausdorff divergence
//'
//' Extending the symmetrical dual solution to the support of the other input measure
//'
//' @param costMatrix A numeric matrix.
//' @param f A numeric vector
//' @param lambda Parameter for the supply proxdiv function
//' @param Div Parameter indicating the divergence function to use for the supply proxdiv function
//' @param eps A numeric vector of decreasing epsilon values.
//' @param param1 numeric Value
//' @param param2 numeric value
//' @param distribution num distrie
//'
//' @return The optimal transport plan
//' @noRd
//[[Rcpp::export]]
Rcpp::NumericVector Hausdorff_Vec_Rcpp(Rcpp::NumericMatrix costMatrix,
                                       Rcpp::NumericVector& distribution,
                                       Rcpp::NumericVector& f,
                                       double lambda,
                                       double param1,
                                       double param2,
                                       int Div,
                                       double eps){
    
    int Nx = costMatrix.nrow();
    int Ny = costMatrix.ncol();
    
    Rcpp::NumericMatrix temp(Nx,Ny);
    Rcpp::NumericVector g(Ny);
 
    Rcpp::NumericVector logDistribution = Rcpp::log(distribution);
    g = parallelVectorLse(logDistribution, f, costMatrix, eps, Ny);
    g = - aprox(lambda, g, eps, Div, param1, param2);
    
    return g;
    
}
