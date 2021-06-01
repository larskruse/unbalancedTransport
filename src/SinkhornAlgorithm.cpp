#include <algorithm>
#include <chrono>
#include <Rcpp.h>




//' init the lambert w function
//' @param x numeric vector
//' @return initial value 
//' @noRd
//[[Rcpp::export]]
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
    
    Rcpp::Rcout << "init z: " << z << "\n";
    
    Rcpp::NumericVector c;
    Rcpp::NumericVector b;
    
    // TODO: change eps
    double eps = std::numeric_limits<double>::epsilon();
    
    for(int i = 0; i < 4; i++){
        
        c = z*(Rcpp::log(z+eps)+z-x)/(1+z);
        b = -1/(z*(1+z));
        z = Rcpp::pmax(z-c/(1-0.5*c*b),eps);
    }
    
    Rcpp::Rcout << "end z: " << z << "\n";
    
    return(z);
}



//' computing the inital values for f and g
//' @param x numeric vector
//' @return the function value
//' @noRd
Rcpp::NumericVector inital(Rcpp::NumericVector& x){
    
    
    Rcpp::NumericVector z = lambertInit(x);
    Rcpp::NumericVector c;
    Rcpp::NumericVector b;
    
    // TODO: change eps
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
Rcpp::NumericVector aprox(double lambda, Rcpp::NumericVector& p, double eps,
                          int DivFun, double param1, double param2){
    Rcpp::NumericVector temp(p.length());
    
    if (DivFun == 1){
        temp = lambda/(lambda+eps) * p;
        
        return temp;
    
    }else if(DivFun == 2){
        
        temp = Rcpp::pmax(p,-lambda);
        temp = Rcpp::pmin(temp,lambda);
        
        return(temp);
    }else if(DivFun == 3){
        
        temp = Rcpp::pmin(0, p-(eps*std::log(param1)));
        temp = Rcpp::pmax(temp, p-(eps*std::log(param2)));
        
        return(temp);
        
    // Power
    }else if(DivFun == 4){
        
        
        
        temp = p/(eps*(1-param1));
        Rcpp::Rcout << "temp : " << temp << "\n";
        temp = temp + lambda/eps;
        Rcpp::Rcout << temp << "\n";
        temp = temp + std::log(lambda/eps);
        Rcpp::Rcout << temp << "\n";
        
        Rcpp::Rcout << "param1 : " << param1 << "\n";
        Rcpp::Rcout << "eps : " << eps << "\n";
        
        temp = (1-param1)*(lambda-eps*lambertWFunctionLog(temp));
        
        
        return(-temp);
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
Rcpp::NumericVector matMul(Rcpp::NumericMatrix& mat, Rcpp::NumericVector& vec, int Nx, int Ny){
    
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
Rcpp::NumericVector init_vectors(double lambda, Rcpp::NumericMatrix costMatrix,
                                 Rcpp::NumericVector& distribution, 
                                 Rcpp::NumericVector& secDistribution, int DivFun,
                                 double param1, double param2, int Nx, int Ny,
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
Rcpp::List Sinkhorn_Rcpp(Rcpp::NumericMatrix costMatrix, Rcpp::NumericVector& supply,
                                  Rcpp::NumericVector& demand, double lambdaSupply, double param1Supply,
                                  double param2Supply, double lambdaDemand, double param1Demand,
                                  double param2Demand,int DivSupply, int DivDemand,
                                  int iterMax, double eps, double tol ){
                                  
                                 
   
   
    // number of points in the reference measures
    int Nx = supply.length();
    int Ny = demand.length();
    
    Rcpp::NumericVector logSup = Rcpp::log(supply);
    Rcpp::NumericVector logDem = Rcpp::log(demand);

    // initializing vectors
    Rcpp::NumericVector f = init_vectors(lambdaSupply, costMatrix, supply, demand, DivSupply, param1Supply, param2Supply, Nx, Ny, eps);
    Rcpp::NumericVector g = init_vectors(lambdaDemand, Rcpp::transpose(costMatrix), demand, supply, DivDemand, param1Demand, param2Demand, Ny, Nx, eps);
    
    
    // Rcpp::Rcout << "init: \n" << "f: " << f << "\n" << "g: " << g << "\n\n";
    
    Rcpp::NumericVector f_prev;
    //Rcpp::NumericVector temp;
    Rcpp::NumericMatrix temp(Nx,Ny);
    
    
    Rcpp::Environment LogSumExp("package:logSumExp");
    Rcpp::Function lse = LogSumExp["colLogSumExps"]; 
    
    Rcpp::NumericMatrix transportPlan(Nx,Ny);

    
    // ProfilerStart("sink.log");
    
    
    
    for(size_t k=0; k < iterMax; k++){
        
        f_prev = Rcpp::clone(f);
        //Rcpp::Rcout << "k:" <<k << "\n";

       // Rcpp::Rcout << "f: " << f << "\n";
        
        for(int j = 0; j < Ny; j++){
            temp(Rcpp::_,j) = logSup + (f-costMatrix(Rcpp::_,j))/eps;
        }
        
        
        g = lse(temp);
        g = -eps*g;
        
        // Rcpp::Rcout << "g_j\n";
        // Rcpp::Rcout << g << "\n";
        
        g = aprox(lambdaSupply, g, eps, DivSupply, param1Supply, param2Supply);
        // Rcpp::Rcout << g << "\n";
        
        
        for(int i = 0; i < Nx; i++){
            temp(Rcpp::_,i) = log(demand)+(g-costMatrix(i,Rcpp::_))/eps;
        }
        
        f = lse(temp);
        f = -eps*f;
        
        // Rcpp::Rcout << "f_i\n";
        // Rcpp::Rcout << f << "\n";
        
        f = aprox(lambdaDemand, f, eps, DivDemand, param1Demand, param2Demand);
        // Rcpp::Rcout << f << "\n\n";
        
        
        if(Rcpp::max(Rcpp::abs(f-f_prev)) < tol){
            Rcpp::Rcout << "converged at: " << k << " \n";

            break;
        }
        
    }
    // ProfilerStop();
    
    
    
    
    
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny ; j++){
            transportPlan(i,j) = exp((f(i) + g(j) - costMatrix(i,j))/eps);
        }
    }
    

    // returnING the transport plan
    // since the absorbtion is called in the last iteration of the loop,
    // the transport plan is equal to the kernel.
    return Rcpp::List::create(Rcpp::Named("TransportPlan") = transportPlan,
                              Rcpp::Named("dual_f") = f,
                              Rcpp::Named("dual_g") = g);

}












//' The symmetric stabilized Scaling Algorithm
//'
//' C++ implementation of the log-domain stabilized Version of the Scaling
//' Algorithm.
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
Rcpp::NumericVector Hausdorff_Vec_Rcpp(Rcpp::NumericMatrix costMatrix,Rcpp::NumericVector& distribution, Rcpp::NumericVector& f,
                         double lambda, double param1,
                         double param2, int Div, double eps){
    int Nx = costMatrix.nrow();
    int Ny = costMatrix.ncol();
    
    Rcpp::NumericMatrix temp(Nx,Ny);
    Rcpp::NumericVector g(Ny);
 
    Rcpp::Environment LogSumExp("package:logSumExp");
    Rcpp::Function lse = LogSumExp["colLogSumExps"]; 
 
    Rcpp::NumericVector logDistribution = Rcpp::log(distribution);
 
    Rcpp::Rcout << "f: " << f << "\n" << "g: " << g << "\n\n";
    Rcpp::Rcout << "temp: " << temp << "\n\n";
 
    for(int j = 0; j < Ny; j++){
        temp(Rcpp::_,j) = logDistribution + (f-costMatrix(Rcpp::_,j))/eps;
    }
    Rcpp::Rcout << "temp: " << temp << "\n\n";
    
 
    g = lse(temp);
    g = -eps*g;
    
    
    g = aprox(lambda, g, eps, Div, param1, param2);
    
    return g;
    
}