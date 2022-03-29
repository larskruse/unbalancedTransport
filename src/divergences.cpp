#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>

//' Elementwise division of two vectors with 'x/0 = 0'
//'
//' @param a The dividend vector
//' @param b The divisor vector
//' @return solution of elementwise division of a and b
//' @noRd
Eigen::VectorXd div0(Eigen::VectorXd a, Eigen::VectorXd b){
    
    // if all values in b are unequal 0,
    // the result is the standard elementwise division
    if((b.array() > 0).all()){
        
        return a.array()/b.array();
        
    }else{
        
        // compute the standard division and check for elements that
        // violate the rule x/0 = 0
        Eigen::VectorXd res = a.array()/b.array();
        
        for(int i = 0; i < a.size(); i++){
            
            if(b(i) == 0){
                
                res(i) = 0;
            }

        }
        return res;
    }

}
//' Elementwise multiplication of two vectors with infinity * 0 = 0
//'
//' @param a The dividend vector
//' @param b The divisor vector
//' @return solution of elementwise division of a and b
//' @noRd
Eigen::VectorXd axb0(Eigen::VectorXd a, Eigen::VectorXd b){
    
    // if all values in b are unequal 0, the result is the
    // standard elementwise division
    if((a.array() > 0).all()){
        
        return a.array()*b.array();
        
    }else{
        
        // check for elements that are zero and set the corresponding result to 0
        Eigen::VectorXd res = a.array()*b.array();
        
        for(int i = 0; i < a.size(); i++){
            
            if(a(i) == 0){
                
                res(i) = 0;
                
            }
            
        }
        return res;
    }
    
}






//' Computing the divergence functions values
//'
//' @param r input vector
//' @param s comparision vector
//' @param DivFun kind of function to use
//' @param param1 lambda or alpha
//' @param param2 beta or 0
//' @param param3 beta or 0
//' @noRd
double vectorDivergence (Eigen::VectorXd r,
                         Eigen::VectorXd s,
                         const int DivFun,
                         const double param1,
                         const double param2 = 0,
                         const double param3 = 0){
    
    // Kullback-Leibler
    Eigen::VectorXd temp;
    
    if (DivFun == 1){
        
        temp = div0(r,s).array().log();
        temp = axb0(r,temp);
        temp = (temp.array() - r.array() + s.array());
        return(param1 * temp.sum());
        
    // Total variation
    }else if(DivFun == 2){
        
        Eigen::VectorXd temp = (r.array()-s.array()).array().abs();
        return(param1 * temp.sum());
        
    // Range Constraint
    }else if(DivFun == 3){
        
        Eigen::VectorXd temp1 = (param2*s.array()).array()-r.array();
        temp1 = temp1.array().max(Eigen::VectorXd::Zero(r.size()).array());
        
        Eigen::VectorXd temp2 = (r.array()-(param3*s.array()).array());
        temp2 = temp2.array().max(Eigen::VectorXd::Zero(r.size()).array());
        
        temp1 = (temp1.array().exp() +temp2.array().exp()).array()-2;
        
        return(temp1.sum());
        
    // Power divergence
    }else if(DivFun == 4){
        
        if(param2 == 0){
            
            temp = div0(r,s).array().log();
            temp = axb0(s,temp);
            temp = r.array() - s.array() - temp.array();
            
            for(int i = 0; i < temp.size(); i++){
                
                if(s(i) == 0){
                    
                    temp(i) = r(i);
                    
                }
                
            }
            
            return(param1 * temp.sum());
            
        }else{
            
            double paramS = param2/(param2 - 1);
            temp = div0(r,s);
            temp = temp.array().pow(paramS) - paramS*(temp.array() - 1 ) - 1;
            temp = temp.array()*s.array();
            
            for(int i = 0; i < temp.size(); i++){
                
                if(s(i) == 0){
                    
                    temp(i) = -paramS * r(i);
                    
                }
                
            }
            
            return(param1/(paramS * (param2-1)) * temp.sum());
            
        }
        
    }else{
    
        return(0);
        
    }
    
}





//' Computing the legendre transformations of the divergence functions 
//'
//' @param p input vector
//' @param u comparision vector
//' @param DivFun kind of function to use
//' @param param1 lambda or alpha
//' @param param2 beta or 0
//' @param param3 beta or 0
//' @noRd
double fVectorDivergence (const Eigen::VectorXd p,
                          const Eigen::VectorXd u,
                          const int DivFun,
                          const double param1,
                          const double param2 = 0,
                          const double param3 = 0){
    
    Eigen::VectorXd temp;
    // Kullback-Leibler
    if (DivFun == 1){
        
        temp = (u.array()/param1).array().exp().array()-1;
        temp = axb0(p,temp);
        return(param1  * temp.sum());
    // Total variation
    }else if(DivFun == 2){
        
        temp = p.array().min((-p).array().max(p.array()*u.array()/param1));
        return(param1 * temp.sum());
        
    // Range Constraint
    }else if(DivFun == 3 && param2 <= param3 && 0 <= param2){
        
        return((((param2*p.array()*u.array())).max(param3*p.array()*u.array())).sum());
        
    // Power divergence
    }else if (DivFun == 4){
        
        if(param2 == 0.5){
            
            temp = 1 - u.array()/param1;
            temp = temp.array().log();
            temp = p.array() * temp.array();
            return(-param1*temp.sum());
            
        }else{
        
            temp = 1+u.array()/(param1*(param2 - 1));
            temp = temp.array().pow(param2)-1;
            temp = param1*(param2-1)/param2 *temp.array();
            temp = p.array() * temp.array();
                
            return(temp.sum());
            
        }
        
    }else{
        
        return(0);
        
    }
    
    
}

//' Computing the entropy term of the dual solution
//'
//' @param u dual potential 
//' @param v dual potential
//' @param eps regularization parameter
//' @param Kernel the Kernel
//' @noRd
double dualSolSummand(const Eigen::VectorXd u,
                      const Eigen::VectorXd v,
                      const double eps,
                      const Eigen::MatrixXd Kernel){
    
    size_t Nx = u.size();
    size_t Ny = v.size();
    
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(Nx, Ny);
    
    for(size_t i = 0; i < Nx; i++){
        
        for(size_t j = 0; j < Ny ; j++){
            
            K(i,j) = (u(i) + v(j));
            
        }
        
    }
   

    K = K.array()/eps;
    K = K.array().exp();
    K = K.array() - 1;
    K = K.array() * Kernel.array();
    
    return(eps * K.sum());
    
}

//' Computing the entropy term of the dual solution for the Sinkhorn algorithm
//'
//' @param u dual potential 
//' @param v dual potential
//' @param eps regularization parameter
//' @param Kernel the Kernel
//' @noRd
double dualSolSummandSink(const Eigen::VectorXd u,
                          const Eigen::VectorXd v,
                          const double eps,
                          const Eigen::MatrixXd costMatrix,
                          const Eigen::MatrixXd supdem){
    
    size_t Nx = u.size();
    size_t Ny = v.size();
    
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(Nx, Ny);
    
    for(size_t i = 0; i < Nx; i++){
        
        for(size_t j = 0; j < Ny ; j++){
            
            K(i,j) = (u(i) + v(j) - costMatrix(i,j));
            
        }
        
    }

    K = K.array()/eps;
    K = K.array().exp();
    K = K.array()-1;
    K = K.array() *supdem.array();
    
    return(eps*K.sum());
    
}




//' Updating the Kernel
//'
//' Calculating and updating the log-domain stabilized kernel.
//' For 0 vectors u and v it calculates the Gibbs kernel.
//'
//' @param u A numeric vector
//' @param v A numeric vector
//' @param eps The epsilon value
//' @param costMatrix A numeric matrix
//' @return The updated kernel
//' @noRd
Eigen::MatrixXd updateK(const Eigen::VectorXd u,
                        const Eigen::VectorXd v,
                        const double eps,
                        const Eigen::MatrixXd costMatrix){
    
    size_t Nx = u.size();
    size_t Ny = v.size();
    
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(Nx, Ny);
    
    for(size_t i = 0; i < Nx; i++){
        
        for(size_t j = 0; j < Ny ; j++){
            
            K(i,j) = exp((u(i) + v(j) - costMatrix(i,j))/eps);
            
        }
        
    }
    
    return K;
    
}
