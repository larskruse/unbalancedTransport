

#' @title legendre Entropy
#' @param lambda num value
#' @param x num value
#' @param DivFun div fun value
#' @param param1 num val
#' @param param2 num val
#' @noRd
#'
legendre_entropy <- function(lambda, x, DivFun, param1 = 0, param2 = 0){
    
    if(DivFun == "KL"){
        
        return(lambda*(exp(x/lambda) - 1))
        
    }else if(DivFun == "Power"){
        print("le power")
        print(x)
        
        if(param1 == 0){
            return( -lambda *log(1-(x/lambda)))
        }else{
            return(lambda*(1- 1/param1) * ((1+x/(lambda *(param1 - 1)))**param1 - 1) )
        }
        
    }else if(DivFun == "RG"){
        return(pmax(param1 * x, param2*x))
        
    }else if(DivFun == "TV"){
        return(x)
    }
    
    
}

#' @title legendre grad
#' @param lambda num value
#' @param x num value
#' @param DivFun div fun value
#' @param param1 num val
#' @param param2 num val
#'
#' @noRd
grad_legrende <- function(lambda, x, DivFun, param1 = 0, param2 = 0){
    
    if(DivFun == "KL"){
        
        return(exp(x/(lambda)))
        
    }else if(DivFun == "Power"){
        
        return((1-(x/(lambda * (1-param1))))**(param1-1))
        
    }else if(DivFun == "RG"){
        return(pmax(-param1*sign(x), param1 * sign(x)))
        
    }else if(DivFun == "TV"){
        return(rep(1, length(x)))
        
        
    }
    
}

#' @title exp exponent
#' @param f num vec
#' @param g num  vec 
#' @param C cost mat
#'
#' @noRd
expC <- function(f, g, C){
    
    res = matrix(0, length(f), length(g))
    
    for(i in 1:length(f)){
        for(j in 1:length(g)){
            
            res[i,j] <- f[i]+g[j]-C[i,j] 
            
        }
        
    }
    
    return(res)
    
}

#' @title Sinkhorn Divergence
#' 
#' A function to calculate the sinkhorn divergence for regularized unbalanced optimal transport.
#' 
#' For supply and demand measures \eqn{\alpha}{a} and \eqn{\beta}{b} the Sinkhorn divergence is defined as 
#' \eqn{S_{\varepsilon} = OT_\varepsilon(\alpha,\beta) - \frac{1}{2}OT_\varepsilon(\alpha,\alpha) -
#' \frac{1}{2}OT_\varepsilon(\beta,\beta) + \frac{\varepsilon}{2}(m(\alpha) - m(\beta))^2 }{S_eps(a,b) 
#' = OT_eps(a,b) - 0.5 OT_eps(a,a) - 0.5 OT_eps(b,b) + 0.5*eps (m(a) - m(b))²}
#' 
#' with transport cost \eqn{OT_\varepsilon(\alpha,\beta)}{OT_eps(a,b)}, regularization parameter \eqn{\varepsilon}{eps}
#' and supply and demand mass \eqn{m(\alpha)}{m(a)} and \eqn{m(\beta)}{m(b)}.
#' 
#' The transport costs are calculated using the Sinkhorn algorithm \code{\link[unbalancedTransport]{sinkhornAlgorithm}}.
#' 
#' @param supplyList A list containing the information about the supply measure, 
#' divergence function and parameter. The first element is the supply measure itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation and "RG" for range constraint) followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL divergence \eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()} and the 
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The same structure is used for the TV divergence: \eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()} and the 
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The divergence function associated with the range constraint needs two parameters, that define the upper and 
#' lower bound. 
#' 
#' If the cost matrix is not provided, the support of the supply measure has to be provied as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' 
#' [REF]
#' 
#' @param demandList A list containing the information about the demand measure. It has to have the same structure as the supplyList.
#' @param eps num val
#' @param iterMax The maximum number of algorithm iterations
#' @param tol  (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-5.
#' @param method (optional) Determines the method that is used to compute the cost matrix.
#' \itemize{
#' \item "euclidean"
#' \item "minkowski"
#' \item "maximum" 
#' }
#' The default value is "euclidean".
#' @param exp (optional) The exponent that is applied to the cost matrix. Can be used to compute quadratic cost. The default value is 1.
#' @param p (optional) Parameter for the minkwski cost function. Can be omitted if either "euclidean" or "maximum" is used. The default value is 2.
#' @param wfr (optional) Computes the cost matrix needed for the Wasserstein-Fisher-Rao distance \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}.
#' The default value is "false". 
#' @param Cxy (optional) 
#' @param Cxx (optional)
#' @param Cyy (optional)
#'
#' @export
sinkhorn_divergence <- function(supplyList, demandList, eps, iterMax = 1000, tol = 1e-5, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){

    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy) | is.null(Cyy) | is.null(Cxx)){
        
        X <- supplyList[[lenSup]]
        Y <- demandList[[lenDem]]
        
        
        if(is.null(Cxy)){
            Cxy <- costMatrix(X, Y , method, exp, wfr, p)
        }
        
        if(is.null(Cxx)){
            Cxx <- costMatrix(X, X, method, exp, wfr, p)
        }
        
        if(is.null(Cyy)){
            Cyy <- costMatrix(Y, Y , method, exp, wfr, p)
        }
        
    }
    
    
    if(supplyList[[2]] != demandList[[2]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithm(supplyList, demandList, iterMax, eps, tol, costMatrix = Cxy)
    res_x = sinkhornAlgorithm(supplyList, supplyList, iterMax, eps, tol, costMatrix = Cxx)
    res_y = sinkhornAlgorithm(demandList, demandList, iterMax, eps, tol, costMatrix = Cyy)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
    f_x1 = res_x$dual_f
    f_x2 = res_x$dual_g
    
    g_y1 = res_y$dual_f
    g_y2 = res_y$dual_g
    
    
    if(supplyList[[2]] == "TV"){
        
        
        func = sum(supplyList[[1]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) + sum(demandList[[1]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
        
        
        supxdem = supplyList[[1]] %*% t(demandList[[1]])
        supxsup = supplyList[[1]] %*% t(supplyList[[1]])
        demxdem = demandList[[1]] %*% t(demandList[[1]])
        

        func = func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)
        
        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
      
        supxdem = supplyList[[1]] %*% t(demandList[[1]])
        supxsup = supplyList[[1]] %*% t(supplyList[[1]])
        demxdem = demandList[[1]] %*% t(demandList[[1]])
        
        
       
        func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)-
                                          0.5 * (-legendre_entropy(0, -f_x1, supplyList[[2]], param1, param2))-
                                          0.5 * (-legendre_entropy(0, -f_x2, supplyList[[2]], param1, param2)))) +
                                        sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2) -
                                            0.5 * (-legendre_entropy(0, -g_y1, supplyList[[2]], param1, param2)) -
                                            0.5 * (-legendre_entropy(0, -g_y2, supplyList[[2]], param1, param2))))
                  
        func = func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        
        return(func)

        
    }else{
        
        param1 <- 0
        
        if(supplyList[[2]] == "Berg"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- 0
            
            
        }else if(supplyList[[2]] == "Hellinger"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- -1
            
        }else if(supplyList[[2]] == "Power"){
            
            param1 <- supplyList[[4]]
        }
        
        
        
        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[2]], param1)
        outf_xx <- -legendre_entropy(supplyList[[3]], -f_x1, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_x1, supplyList[[2]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1)
        outg_yy <- -legendre_entropy(supplyList[[3]], -g_y1, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_y1, supplyList[[2]], param1)
      
        
        func <- sum(supplyList[[1]]*(outf_xy-outf_xx)) + sum(demandList[[1]] * (outg_xy - outg_yy))
        
        return(func)
    }
    
    
    
}





#' @title Regularized transport cost
#' 
#' A function to calculate the regularized transport cost using the dual optimal optimal potentials.
#' 
#' The regularized transport cost is defined as \eqn{OT_\varepsilon = sup_{f,g} -<\alpha, \phi^*(-f)> - 
#' <\beta, \phi^*(-g)> - \varepilon <\alpha \times \beta , \exp{(f+g-C)/\varepsilon}- 1>}
#' 
#' with supply and demand distributions \eqn{\alpha}{a} and \eqn{\beta}{b}, optimal dual potentials \eqn{f} and \eqn{g}, 
#' Legendre transforms of the divergence function \eqn{\phi^*}, cost matrix \eqn{C} and regularization parameter \eqn{\varepsilon}{eps}.
#' 
#' The function uses the Sinkhorn algorithm to compute the optimal dual potentials.
#' 
#' 
#' 
#' [REF]
#' 
#' @param supplyList A list containing the information about the supply measure, 
#' divergence function and parameter. The first element is the supply measure itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation and "RG" for range constraint) followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL divergence \eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()} and the 
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The same structure is used for the TV divergence: \eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()} and the 
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The divergence function associated with the range constraint needs two parameters, that define the upper and 
#' lower bound. 
#' 
#' If the cost matrix is not provided, the support of the supply measure has to be provied as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' @param demandList A list containing the information about the demand measure. It has to have the same structure as the supplyList.
#' 
#' @param maxIteration The maximum number of iterations.
#' @param epsVector A vector containing a decreasing sequence of epsilon values. If no epsilon scaling is needed a vector with a single 
#' value can be used.
#' @param tol (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-10.
#' @param algorithm (optional) Determines which algorithm to use. Possible values are "sinkhorn"
#'  for \code{\link[unbalancedTransport]{sinkhornAlgorithm}} and "scaling" for \code{\link[unbalancedTransport]{scalingAlgorithm}}
#' @param method (optional) Determines the method that is used to compute the cost matrix.
#' \itemize{
#' \item "euclidean"
#' \item "minkowski"
#' \item "maximum" 
#' }
#' The default value is "euclidean".
#' @param exp (optional) The exponent that is applied to the cost matrix. Can be used to compute quadratic cost. The default value is 1.
#' @param p (optional) Parameter for the minkwski cost function. Can be omitted if either "euclidean" or "maximum" is used. The default value is 2.
#' @param wfr (optional) Computes the cost matrix needed for the Wasserstein-Fisher-Rao distance \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}.
#' The default value is "false". 
#' @param Cxy (optional) Instead of having the algorithm compute the cost matrix, a custom cost matrix can be passed to the algorithm. 
#' @export
#' 
regularized_ot <- function(supplyList, demandList, eps, iterMax = 100, tol = 1e-3, algorithm = "sinkhorn", method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL){
    
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
    }
    
    if(supplyList[[2]] != demandList[[2]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithm(supplyList, demandList,
                                       iterMax, eps, tol, costMatrix = Cxy)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
  
    if(supplyList[[2]] == "TV"){
  
        func = sum(supplyList[[1]] * f_xy) + sum(demandList[[1]] * g_xy)
      
        expFun <- supplyList[[1]] %*% t(demandList[[1]]) * (1-exp(expC(f_xy,g_xy,Cxy)/eps))
        
        
        func = func + sum(eps * expFun)
        return(func)
        
       
        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        
        
        func = sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)))+ 
            sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2)))

        func = func + sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        
        return(func)
        

    
    }else{
        
        param1 <- 0
        
        if(supplyList[[2]] == "Berg"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- 0
            
            
        }else if(supplyList[[2]] == "Hellinger"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- -1
            
        }else if(supplyList[[2]] == "Power"){
            
            param1 <- supplyList[[4]]
        }
        

        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[2]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1)
        
        res = sum(supplyList[[1]] * outf_xy) + sum(demandList[[1]] * outg_xy) + 
            eps*(sum(supplyList[[1]]) * sum(demandList[[1]]))
        
        return(res)
        
    }

}





#' @title regularized output
#' @param supplyList sup list
#' @param demandList dem list
#' @param eps num val
#' @param costMatrix cost matrix
#' @param f_xy dual va
#' @param g_xy dual va
#' @noRd
#' 
regularized_ot_intern <- function(supplyList, demandList, f_xy, g_xy, eps, costMatrix){



    if(supplyList[[2]] == "TV"){
        
        func = sum(supplyList[[1]] * f_xy) + sum(demandList[[1]] * g_xy)
        
        expFun <- supplyList[[1]] %*% t(demandList[[1]]) * (1-exp(expC(f_xy,g_xy,costMatrix)/eps))
        
        
        func = func + sum(eps * expFun)
        return(func)
        
        
        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        
        
        func = sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)))+ 
            sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2)))
        
        func = func + sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,costMatrix)/eps))))
        
        return(func)
        
        
        
    }else{
        
        param1 <- 0
        
        if(supplyList[[2]] == "Berg"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- 0
            
            
        }else if(supplyList[[2]] == "Hellinger"){
            supplyList[[2]] <- "Power"
            demandList[[2]] <- "Power"
            
            param1 <- -1
            
        }else if(supplyList[[2]] == "Power"){
            
            param1 <- supplyList[[4]]
        }
        
        
        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[2]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1)
        
        res = sum(supplyList[[1]] * outf_xy) + sum(demandList[[1]] * outg_xy) + 
            eps*(sum(supplyList[[1]]) * sum(demandList[[1]]))
        
        
        
        return(res)
        
    }
    
}







#' 
#' #' @title hausdorff div
#' #' @param supplyList sup list
#' #' @param demandList dem list
#' #' @param eps num val
#' #' @param iterMax num val
#' #' @param tol num val
#' #' @param method cost method
#' #' @param exp cost exponent
#' #' @param p cost parameter
#' #' @param wfr wfr param
#' #' @param Cxx cost matrix
#' #' @param Cyy cost matrix
#' #' @param Cxy cost matrix 
#' #' @export
#' #'
#' hausdorff_divergence <- function(supplyList, demandList, eps, iterMax = 1000, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
#'                                 Cxy = NULL, Cxx = NULL, Cyy = NULL){
#' 
#'     lenSup <- length(supplyList)
#'     lenDem <- length(demandList)
#' 
#'     if(is.null(Cxy)){
#'         Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
#'     }
#' 
#'     if(is.null(Cxx)){
#'         Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], method, exp, wfr, p)
#'     }
#' 
#'     if(is.null(Cyy)){
#'         Cyy <- costMatrix(demandList[[lenDem]], demandList[[lenDem]], method, exp, wfr, p)
#'     }
#' 
#' 
#' 
#'     if(supplyList[[1]] != demandList[[1]]){
#' 
#'         print("Please chose the same entropy for supply and demand.")
#' 
#'     }
#' 
#'     
#'     reg <- 0
#'     param1 <- 0
#'     param2 <- 0
#'     
#'     
#' 
#'     
#'     supply <- supplyList[[2]]
#'     demand <- demandList[[2]]
#'     
#' 
#'     res_x = sinkhornAlgorithm(Cxx, supplyList, supplyList, iterMax, eps, tol)
#'     res_y = sinkhornAlgorithm(Cyy, demandList, demandList, iterMax, eps, tol)
#' 
#'     f_x = res_x$dual_f
#'     g_y = res_y$dual_g
#' 
#'     if(supplyList[[1]] == "KL"){
#'         div <- 1
#'         reg <- supplyList[[3]]
#'         
#'     }else if(supplyList[[1]] == "TV"){
#'         div <- 2
#'         reg <- supplyList[[3]]
#'     }else if(supplyList[[1]] == "RG"){
#'         div <- 3
#'         param1 <- supplyList[[3]]
#'         param2 <- supplyList[[4]]
#'         
#'         if(param1 < 0 || param2 < param1){
#'             stop("0 <= Alpha <= Beta")
#'         }
#'         
#'     }else if(supplyList[[1]] == "Berg"){
#'         supplyList[[1]] <- "Power"
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- 0
#'     }else if(supplyList[[1]] == "Hellinger"){
#'         supplyList[[1]] <- "Power"
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- -1
#'     }else if(supplyList[[1]] == "Power"){
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- supplyList[[4]]
#'     }else{
#'         stop("Please supply a divergence")
#'     }
#'     
#'     
#'     g_xy = Hausdorff_Vec_Rcpp(Cxy, supply, f_x, reg, param1, param2, div, eps)
#'     f_xy = Hausdorff_Vec_Rcpp(t(Cxy), demand, g_y, reg, param1, param2, div, eps)
#'     
#'     
#'     
#'     print("f_x, g_y, ...")
#'     print(f_x)
#'     print(g_y)
#'     print(g_xy)
#'     print(f_xy)
#'     
#'     
#'     print("part outs")
#'     print(legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     print(grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     
#'     print("outs")
#'     
#'     print(legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     print(legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_xy,supplyList[[1]], param1, param2))
#'     
#'     print(legendre_entropy(supplyList[[3]], -g_y, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_y,supplyList[[1]], param1, param2))
#'     print(legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2))
#' 
#'     res <- sum(supplyList[[2]] * (legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) -
#'                                       legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, param2) - eps * grad_legrende(supplyList[[3]], -f_xy,supplyList[[1]], param1, param2))) +
#'         sum(demandList[[2]] *(legendre_entropy(supplyList[[3]], -g_y, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_y,supplyList[[1]], param1, param2)-
#'                                   legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2) - eps * grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2)))
#' 
#'     
#'     
#'     
#'     return(res)
#' }
