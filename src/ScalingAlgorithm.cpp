#include <algorithm>
#include <chrono>
#include <RcppEigen.h>


// Elementwise division of two vectors with 'x/0 = 0'
//
// @param a The dividend vector
// @param b The divisor vector
// @return solution of elementwise division of a and b

Eigen::VectorXd div0(Eigen::VectorXd a, Eigen::VectorXd b){
  // if all values in b are unequal 0, the result is the standard elementwise division
  if((b.array() > 0).all()){
    return a.array()/b.array();
  }else{
    // compute the standard division and check for elements that violate the rule x/0 = 0
    Eigen::VectorXd res = a.array()/b.array();
    for(int i = 0; i < a.size(); i++){
      if(b(i) == 0){
        res(i) = 0;
      }

    }
    return res;
  }

}

//
//function fdiv(φ,x,p,dx)
//  if φ[:type] == "OT" || (in(φ[:type],["KL","TV"]) && φ[:param] == Inf)
//    return sum( axb0(dx,exp(abs(x-p))-1) ) # barrier approx
//  elseif φ[:type] == "KL"
//    λ = φ[:param]
//    return  λ*sum( axb0(dx,axb0(x,log(div0(x,p))) - x + p ))
//  elseif φ[:type] == "TV"
//    λ = φ[:param]
//    return λ*sum( dx.*abs(x-p) )
//  elseif φ[:type] == "RG"
//    β1,β2 = φ[:param]
//    @assert 0 <= β1 <= β2 < Inf
//    return sum( axb0(dx, exp(max(0,β1*p-x))+ exp(max(0,x-β2*p))-2)) #barrier
//  else
//    error("Type of φ not recognized")
//end
//end


//' Computing the divergence functions values
//'
//' @param r input vector
//' @param s comparision vector
//' @param DivFun kind of function to use
//' @param param1 lambda or alpha
//' @param param2 beta or 0
//'
//' @export
//'
//[[Rcpp::export]]
double vectorDivergence (Eigen::VectorXd r, Eigen::VectorXd s, int DivFun, double param1, double param2 = 0){
  if (DivFun == 1){
    // return  λ*sum(axb0(x,log(div0(x,p))) - x + p )
    Eigen::VectorXd temp = (r.array()*(div0(r,s).array().log()) - r.array() + s.array());
    return(param1 * temp.sum());

   }else if(DivFun == 2){
     Eigen::VectorXd temp = (r.array()-s.array()).array().abs();

     return(param1 * temp.sum());
  // }else if(DivFun == 3 && param1 <= param2 && 0 <= param1){
  //   Eigen::VectorXd temp =  ;
  //   return(temp.sum())
  //   //    return sum( exp(max(0,β1*p-x))+ exp(max(0,x-β2*p))-2)


  }

  return(0);
}




// //
// //function fdivstar(φ,u,p,dx)
// //  if φ[:type] == "OT" || (in(φ[:type],["KL","TV"]) && φ[:param] == Inf)
// //    return sum( axb0(p.*dx,u) )
// //  elseif φ[:type] == "KL"
// //    λ = φ[:param]
// //    return λ*sum( axb0(p.*dx,exp(u/λ)-1) )
// //  elseif φ[:type] == "TV"
// //    λ = φ[:param]
// //    return λ*sum( dx.*min(p, max(-p,axb0(p,u/λ))) )
// //  elseif φ[:type] == "RG"
// //    β1,β2 = φ[:param]
// //    @assert 0 <= β1 <= β2 < Inf
// //    return sum( axb0(dx, max(β1*axb0(p,u), β2*axb0(p,u)) ))
// //  else
// //    error("Type of φ not recognized")
// //  end
// //end
//
// //' Computing the divergence functions values
// //'
// //' @param p input vector
// //' @param u comparision vector
// //' @param DivFun kind of function to use
// //' @param param1 lambda or alpha
// //' @param param2 beta or 0
// //'
// //' @export
// //'
// //[[Rcpp::export]]
// Eigen::VectorXd fVectorDivergence (Eigen::VectorXd p, Eigen::VectorXd u, int DivFun, double param1, double param2 = 0){
//   if (DivFun == 1){
//     return(param1 * (p.array() * (u.array()/param1).array().exp()-1).sum());
//
//   }else if(DivFun == 2){
//     return(param1*(p.array().min(-p.array().max(p.array()*u.array()/param1)).sum()));
//
//   }else if(DivFun == 3 && param1 <= param2 && 0 <= param1){
//     return(((param1*p.array()*u.array())).max(param2*p.array()*u.array()).sum());
//
//   }
//




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
//' @export
//[[Rcpp::export]]
Eigen::VectorXd proxdiv(double lambda, Eigen::VectorXd p, Eigen::VectorXd s, Eigen::VectorXd u, double eps, int DivFun, double alpha, double beta){
  if (DivFun == 1){
    Eigen::VectorXd temp = s.array()*exp(u.array()/lambda);
    Eigen::VectorXd temp1 = div0(p,temp);
    Eigen::VectorXd temp2 = temp1.array().pow(lambda/(lambda+eps));

    return temp2;

  }else if(DivFun == 2){
    return ((lambda-u.array())/eps).array().exp().array().min(div0(p,s).array().max((-(lambda+u.array())/eps).array().exp()));
  }else{
    return((beta*div0(p,s)).array().min(alpha*div0(p,s).array().max((-u.array()/eps).array().exp())));
  }

}

// Updating the Kernel
//
// Calculating and updating the log-domain stabilized kernel. For 0 vectors u and v
// it calculates the Gibbs kernel.
//
// @param u A numeric vector
// @param v A numeric vector
// @param eps The epsilon value
// @param costMatrix A numeric matrix
// @return The updated kernel

Eigen::MatrixXd updateK(Eigen::VectorXd u, Eigen::VectorXd v, double eps, Eigen::MatrixXd costMatrix){

    int size = u.size();

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(size, size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size ; j++){
          K(i,j) = exp((u(i) + v(j) - costMatrix(i,j))/eps);
        }
    }
    return K;

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
//'
//' @return The optimal transport plan
//' @export
// [[Rcpp::export]]

Rcpp::List StabilizedScaling_Rcpp(Eigen::Map<Eigen::MatrixXd> costMatrix,Eigen::Map<Eigen::VectorXd> supply,
                                  Eigen::Map<Eigen::VectorXd> demand, double lambdaSupply, double alphaSupply,
                                  double betaSupply, double lambdaDemand, double alphaDemand, double betaDemand,
                                  int DivSupply, int DivDemand, int iterMax, Eigen::Map<Eigen::VectorXd> epsvec){
  // number of absorptions
  int numAbs = 0;

  // number of points in the reference measures
  int Nx = supply.size();
  int Ny = demand.size();

  // initializing vectors
  Eigen::VectorXd a = Eigen::VectorXd::Ones(Nx);
  Eigen::VectorXd b = Eigen::VectorXd::Ones(Ny);

  // stabilization vectors
  Eigen::VectorXd u = Eigen::VectorXd::Zero(Nx);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(Ny);


  // main loop iteration counter
  int i = 1;

  // epsilon value index counter
  int epsind = 0;

  // setting first epsilon value
  double eps = epsvec(epsind);


  // computing the initial kernel
  // since u and v are 0, the updateK function returns the Gibbs kernel
  Eigen::MatrixXd Kernel = updateK(u, v, eps, costMatrix);

  while(i < iterMax){


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
    if ((abs(a.array()) > 1e+100).any() || (abs(b.array()) > 1e+100).any() ||
        (static_cast<double>(i)/static_cast<double>(iterMax)) > static_cast<double>(epsind + 1)/static_cast<double>(epsvec.size()) ||
        i == iterMax - 1){


      if(round(100*static_cast<double>(i)/static_cast<double>(iterMax)) < 100){
        Rcpp::Rcout << round(100*static_cast<double>(i)/static_cast<double>(iterMax)) << "% done. \n";
      }

      // update number of absorptions
      numAbs = numAbs +1;

      // absorbing a/b in u/v
      u = u.array() + eps*(a.array().log());
      v = v.array() + eps*(b.array().log());



      // updating epsilon
      if((static_cast<double>(i)/static_cast<double>(iterMax)) > static_cast<double>(epsind + 1)/static_cast<double>(epsvec.size())){
        epsind = epsind + 1;
        eps = epsvec(epsind);
      }
      // update Kernel according to u and v
      Kernel = updateK(u, v, eps, costMatrix);
      //reset a and b
      a = Eigen::VectorXd::Ones(Nx);
      b = Eigen::VectorXd::Ones(Ny);

    }

    // update iteration counter
    i = i + 1;

  }

  Rcpp::Rcout << 100 << "% done. \n";
  // returnING the transport plan
  // since the absorbtion is called in the last iteration of the loop,
  // the transport plan is equal to the kernel.
return Rcpp::List::create(Rcpp::Named("TransportMap") = Kernel);

}


