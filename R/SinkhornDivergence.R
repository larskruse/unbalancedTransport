

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
#' @description Calculating the Sinkhorn divergence for regularized unbalanced optimal transport.
#' 
#' 
#' @details  A function to calculating the Sinkhorn divergence for regularized unbalanced optimal transport.
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
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#' 
#' @param supplyList A list containing the information about the supply distribution, 
#' divergence function and parameter. The first element is the supply distribution itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation, "RG" for range constraint, "Power", "Berg" and "Hellinger") followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL divergence (\eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()}),
#' TV divergence (\eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()}) and all Power
#' divergences (\eqn{F_1 = \lambda \cdot Power()}{F1 = lambda * Power()}) the
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The "Power" divergence also needs the conjugate exponent as second parameter.
#' 
#' The divergence function associated with the range constraint needs two parameters, that define the upper and 
#' lower bound. 
#' 
#' 
#' If the cost matrix is not provided, the support of the supply distribution has to be provided as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' 
#' 
#' @param demandList A list containing the information about the demand measure in the same form as the supplyList.
#' @param eps A numeric value for the regularization parameter.
#' @param iterMax (optional) The maximum number of algorithm iterations. The default value is 10000.
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
#' @param Cxy (optional) A cost matrix for transport between the supply and demand distributions.
#' @param Cxx (optional) A cost matrix for transport between the supply and supply distributions.
#' @param Cyy (optional) A cost matrix for transport between the demand and demand distributions.
#' @return The Sinkhorn divergence.
#' @examples 
#' I <- 1000
#' J <- 1000
#' X <- seq(0,1,length.out = I)
#' Y <- seq(0,1,length.out = J)
#' p <- supplyExample
#' q <- demandExample
#'
#' supply <- list(p,X)
#' demand <- list(q,Y)
#'
#' eps <- 1e-3 
#' supplyList <- list(p, "KL", 0.04, X)
#' demandList <- list(q, "KL", 0.04, Y)
#' 
#' sinkhorn_divergence(supplyList, demandList, eps, exp = 2)
#' 
#'
#' @export
sinkhorn_divergence <- function(supplyList, demandList, eps, iterMax = 10000, tol = 1e-5, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){
    
    
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy) | is.null(Cyy) | is.null(Cxx)){

        X <- supplyList[[lenSup]]
        Y <- demandList[[lenDem]]


        if(is.null(Cxy)){
            Cxy <- costMatrix(X, Y , method, exp, p, wfr)
        }

        if(is.null(Cxx)){
            Cxx <- costMatrix(X, X, method, exp, p, wfr)
        }

        if(is.null(Cyy)){
            Cyy <- costMatrix(Y, Y , method, exp, p, wfr)
        }

    }
    
    
    if(supplyList[[2]] != demandList[[2]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy <- sinkhornAlgorithm(supplyList, demandList, eps, iterMax, tol, costMatrix = Cxy, duals = TRUE)
    res_x <- sinkhornAlgorithm(supplyList, supplyList, eps, iterMax, tol, costMatrix = Cxx, duals = TRUE)
    res_y <- sinkhornAlgorithm(demandList, demandList, eps, iterMax, tol, costMatrix = Cyy, duals = TRUE)
    
    # print(res_xy)
    
    f_xy <- res_xy$dual_f
    g_xy <- res_xy$dual_g
    
    f_x1 <- res_x$dual_f
    f_x2 <- res_x$dual_g
    
    g_y1 <- res_y$dual_f
    g_y2 <- res_y$dual_g
    
    
    # print(f_xy)
    # print(g_xy)
    # print(f_x1)
    # print(f_x2)
    # print(g_y1)
    # print(g_y2)


    if(supplyList[[2]] == "TV"){
        
        
        func <- sum(supplyList[[1]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) + sum(demandList[[1]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
        
        # print("funcs")
        # print(func)
        # 
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        supxsup <- supplyList[[1]] %*% t(supplyList[[1]])
        demxdem <- demandList[[1]] %*% t(demandList[[1]])
        
        # print(sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        # print(0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))))
        # print(0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps))))

        func <- func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)
        
        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
      
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        supxsup <- supplyList[[1]] %*% t(supplyList[[1]])
        demxdem <- demandList[[1]] %*% t(demandList[[1]])
        
        # print("funcs")
        # print(sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)-
        #                                  0.5 * (-legendre_entropy(0, -f_x1, supplyList[[2]], param1, param2))-
        #                                  0.5 * (-legendre_entropy(0, -f_x2, supplyList[[2]], param1, param2)))))
        # print(sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2) -
        #                                  0.5 * (-legendre_entropy(0, -g_y1, supplyList[[2]], param1, param2)) -
        #                                  0.5 * (-legendre_entropy(0, -g_y2, supplyList[[2]], param1, param2)))))
        # 
        # print(sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        # print(0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))))
        # print(0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps))))
        # 
        func <- sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)-
                                          0.5 * (-legendre_entropy(0, -f_x1, supplyList[[2]], param1, param2))-
                                          0.5 * (-legendre_entropy(0, -f_x2, supplyList[[2]], param1, param2)))) +
                                        sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2) -
                                            0.5 * (-legendre_entropy(0, -g_y1, supplyList[[2]], param1, param2)) -
                                            0.5 * (-legendre_entropy(0, -g_y2, supplyList[[2]], param1, param2))))
                  
        func <- func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
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
        
        # print("funcs")
        
        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[2]], param1)
        outf_xx <- -legendre_entropy(supplyList[[3]], -f_x1, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_x1, supplyList[[2]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1)
        outg_yy <- -legendre_entropy(supplyList[[3]], -g_y1, supplyList[[2]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_y1, supplyList[[2]], param1)
      
        # print(sum(supplyList[[1]]*(outf_xy-outf_xx)))
        # print(sum(demandList[[1]] * (outg_xy - outg_yy)))
        
        func <- sum(supplyList[[1]]*(outf_xy-outf_xx)) + sum(demandList[[1]] * (outg_xy - outg_yy))
        
        return(func)
    }
    
    
    
}





#' @title Regularized transport cost
#' 
#' @description Calculating the regularized transport cost.
#' 
#' @details This function calculates the regularized optimal transport costs using the \code{\link[unbalancedTransport]{sinkhornAlgorithm}}.
#' The regularized transport cost is defined as \eqn{OT_\varepsilon = sup_{f,g} -<\alpha, \phi^*(-f)> - 
#' <\beta, \phi^*(-g)> - \varepilon <\alpha \times \beta , \exp{(f+g-C)/\varepsilon}- 1>}
#' 
#' with supply and demand distributions \eqn{\alpha}{a} and \eqn{\beta}{b}, optimal dual potentials \eqn{f} and \eqn{g}, 
#' Legendre transforms of the divergence function \eqn{\phi^*}, cost matrix \eqn{C} and regularization parameter \eqn{\varepsilon}{eps}.
#' 
#' 
#' 
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#' @param supplyList A list containing the information about the supply distribution, 
#' divergence function and parameter. The first element is the supply distribution itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation, "RG" for range constraint, "Power", "Berg" and "Hellinger") followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL divergence (\eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()}),
#' TV divergence (\eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()}) and all Power
#' divergences (\eqn{F_1 = \lambda \cdot Power()}{F1 = lambda * Power()}) the
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The "Power" divergence also needs the conjugate exponent as second parameter.
#' 
#' The divergence function associated with the range constraint needs two parameters which define the upper and 
#' lower bounds. 
#' 
#' 
#' If the cost matrix is not provided, the support of the supply and demand distributions has to be provided as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' 
#' 
#' @param demandList A list containing the information about the demand measure in the same form as the supplyList.
#' @param eps A numeric value for the regularization parameter.
#' @param iterMax (optional) The maximum number of algorithm iterations. The default value is 5000
#' @param tol  (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-5.
#' @param method (optional) Determines the method that is used to compute the cost matrix:
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
#' @param Cxy (optional) A cost matrix for transport between the supply and demand distributions.
#' @return The regularized transport cost.
#' @examples 
#' I <- 1000
#' J <- 1000
#' X <- seq(0,1,length.out = I)
#' Y <- seq(0,1,length.out = J)
#' p <- supplyExample
#' q <- demandExample
#'
#' supply <- list(p,X)
#' demand <- list(q,Y)
#'
#' eps <- 1e-3 
#' supplyList <- list(p, "KL", 0.04, X)
#' demandList <- list(q, "KL", 0.04, Y)
#' 
#' regularized_ot(supplyList, demandList, eps, exp = 2)
#' 
#' 
#' @export
#' 
regularized_ot <- function(supplyList, demandList, eps, iterMax = 5000, tol = 1e-5,
                           method = "euclidean", exp = 1, p = 2,  wfr = FALSE, Cxy = NULL){
    
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, p, wfr)
    }
    
    if(supplyList[[2]] != demandList[[2]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy <- sinkhornAlgorithm(supplyList, demandList, eps, iterMax, tol = tol, costMatrix = Cxy, duals = TRUE)
    
    f_xy <- res_xy$dual_f
    g_xy <- res_xy$dual_g
    
  
    if(supplyList[[2]] == "TV"){
  
      
        
        func = sum(supplyList[[1]] * f_xy) + sum(demandList[[1]] * g_xy)
      
        
        expFun <- supplyList[[1]] %*% t(demandList[[1]]) * (1-exp(expC(f_xy,g_xy,Cxy)/eps))
        
        func <- func + sum(eps * expFun)
        return(func)
        
       
        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        
        
        func <- sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2)))+ 
            sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2)))

        func <- func + sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        
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
        
        res <- sum(supplyList[[1]] * outf_xy) + sum(demandList[[1]] * outg_xy) + 
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
        # 
        # print(f_xy)
        # print(g_xy)
        # 
        f_xy[!is.finite(f_xy) & supplyList[[1]] == 0] <- 0
        g_xy[!is.finite(g_xy) & demandList[[1]] == 0] <- 0
        
        # print(f_xy)
        # print(g_xy)
        
        func <- sum(supplyList[[1]] * f_xy) + sum(demandList[[1]] * g_xy)
        
        
        expMat <- (1-exp(expC(f_xy,g_xy,costMatrix)/eps))
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        
        expMat[!is.finite(expMat) & supxdem == 0 ] <- 0
        
        # expFun <- supplyList[[1]] %*% t(demandList[[1]]) * (1-exp(expC(f_xy,g_xy,costMatrix)/eps))
        expFun <- supxdem * expMat
        
        # print("funcs")
        # print(sum(supplyList[[1]] * f_xy))
        # print(sum(demandList[[1]] * g_xy))
        # print(func)
        # print(sum(expFun))
        # print(sum(eps * expFun))
        # 
        
        
        
        
        func <- func + sum(eps * expFun)
        return(func)
        

        
    }else if(supplyList[[2]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        expMat <- (1-exp(expC(f_xy,g_xy,costMatrix)/eps))
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        expMat[!is.finite(expMat) & supxdem == 0 ] <- 0
        
        
        
        outf_xy <- (-legendre_entropy(0, -f_xy, supplyList[[2]], param1, param2))
        outg_xy <- (-legendre_entropy(0, -g_xy, supplyList[[2]], param1, param2))
        
        outf_xy[!is.finite(outf_xy) & supplyList[[1]] == 0] <- 0
        outg_xy[!is.finite(outg_xy) & demandList[[1]] == 0] <- 0
        
        # print("func")
        # print(sum(supplyList[[1]] * outf_xy))
        # print(sum(demandList[[1]] * outg_xy))
        
        func <- sum(supplyList[[1]] * outf_xy)+ 
            sum(demandList[[1]] * outg_xy)
        
        # 
        # print(sum(supxdem * expMat))
        # 
        # print(sum(eps * (supxdem * expMat)))
        
        func <- func + sum(eps * (supxdem * expMat))
        
        # print(func)
        
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
        
        # print("funcs")
        # 
        # print(-legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1, 0))
        # print(- 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[2]], param1))
        # print(-legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1, 0))
        # print(- 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1))
        
        
        outf_xy[!is.finite(outf_xy) & supplyList[[1]] == 0] <- 0
        outg_xy[!is.finite(outg_xy) & demandList[[1]] == 0] <- 0
        
        # print(supplyList[[1]] * outf_xy)
        # print(demandList[[1]] * outg_xy)
        # 
        # print(sum(supplyList[[1]] * outf_xy))
        # print(sum(demandList[[1]] * outg_xy))
        # print(eps*(sum(supplyList[[1]]) * sum(demandList[[1]])))
                
        res = sum(supplyList[[1]] * outf_xy) + sum(demandList[[1]] * outg_xy) + 
            eps*(sum(supplyList[[1]]) * sum(demandList[[1]]))
        
        # print(res)
        
        return(res)
        
    }
    
}







#' @title Hausdorff Divergence
#' 
#' @description Calculating the Hausdorff divergence.
#' 
#' @details This function calculates the Hausdorff divergence using the dual optimal optimal potentials
#' calculated by the \code{\link[unbalancedTransport]{sinkhornAlgorithm}}.
#' The Hausdorff divergence is defined as 
#' \eqn{H_\varepsilon(\alpha, \beta) = <\alpha - \beta, \nabla F_\varepsilon(\alpha) - 
#' \nabla F_\varepsilon(\beta)>}{H_eps(a,b) = <a-b, nabla(F_eps(a)) - nabla(F_eps(b)) >}  
#' with supply and demand distributions \eqn{\alpha}{a} and \eqn{\beta}{b} and
#' \eqn{F_\varepsilon(\alpha) = -\frac{1}{2} OT_\varepsilon(\alpha,\alpha) + 
#' \frac{\varepsilon}{2}(m(\alpha)^2)}{F_eps(a) = -0.5 OT_eps(a,a) + 0.5eps * m(a)^2}.
#' 
#' 
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#' @param supplyList A list containing the information about the supply distribution, 
#' divergence function and parameter. The first element is the supply distribution itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation, "RG" for range constraint, "Power", "Berg" and "Hellinger") followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL divergence (\eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()}),
#' TV divergence (\eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()}) and all Power
#' divergences (\eqn{F_1 = \lambda \cdot Power()}{F1 = lambda * Power()}) the
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' The "Power" divergence also needs the conjugate exponent as second parameter.
#' 
#' The divergence function associated with the range constraint needs two parameters which define the upper and 
#' lower bounds. 
#' 
#' 
#' If the cost matrix is not provided, the support of the supply and demand distributions has to be provided as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' 
#' 
#' @param demandList A list containing the information about the demand measure in the same form as the supplyList.
#' @param eps A numeric value for the regularization parameter.
#' @param iterMax (optional) The maximum number of algorithm iterations. The default value is 5000
#' @param tol  (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-5.
#' @param method (optional) Determines the method that is used to compute the cost matrix:
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
#' @param Cxy (optional) A cost matrix for transport between the supply and demand distributions.
#' @param Cxx (optional) A cost matrix for transport between the supply and supply distributions.
#' @param Cyy (optional) A cost matrix for transport between the demand and demand distributions.
#' @return The Hausdorff divergence.
#' 
#' @examples 
#' I <- 1000
#' J <- 1000
#' X <- seq(0,1,length.out = I)
#' Y <- seq(0,1,length.out = J)
#' p <- supplyExample
#' q <- demandExample
#'
#' supply <- list(p,X)
#' demand <- list(q,Y)
#'
#' eps <- 1e-3 
#' supplyList <- list(p, "KL", 0.04, X)
#' demandList <- list(q, "KL", 0.04, Y)
#' 
#' hausdorff_divergence(supplyList, demandList, eps, exp = 2)
#'
#' @export
hausdorff_divergence <- function(supplyList, demandList, eps, iterMax = 1000, tol = 1e-3,
                                 method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                                 Cxy = NULL, Cxx = NULL, Cyy = NULL){

    lenSup <- length(supplyList)
    lenDem <- length(demandList)

    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, p, wfr)
    }

    if(is.null(Cxx)){
        Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], method, exp, p, wfr)
    }

    if(is.null(Cyy)){
        Cyy <- costMatrix(demandList[[lenDem]], demandList[[lenDem]], method, exp, p, wfr)
    }



    if(supplyList[[2]] != demandList[[2]]){

        print("Please chose the same entropy for supply and demand.")

    }


    reg <- 0
    param1 <- 0
    param2 <- 0



    supply <- supplyList[[1]]
    demand <- demandList[[1]]


    res_x = sinkhornAlgorithm(supplyList, supplyList, eps, iterMax, tol, costMatrix = Cxx, duals = TRUE)
    res_y = sinkhornAlgorithm(demandList, demandList, eps, iterMax, tol, costMatrix = Cyy, duals = TRUE)

    f_x = res_x$dual_f
    g_x = res_x$dual_g
    g_y = res_y$dual_g
    f_y = res_y$dual_f
    
    f_x = 0.5*(f_x+g_x)
    g_y = 0.5*(g_y+f_y)
    

    if(supplyList[[2]] == "KL"){
        div <- 1
        reg <- supplyList[[3]]

    }else if(supplyList[[2]] == "TV"){
        div <- 2
        reg <- supplyList[[3]]
    }else if(supplyList[[2]] == "RG"){
        div <- 3
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]

        if(param1 < 0 || param2 < param1){
            stop("0 <= Alpha <= Beta")
        }

    }else if(supplyList[[2]] == "Berg"){
        supplyList[[2]] <- "Power"
        div <- 4
        reg <- supplyList[[3]]
        param1 <- 0
    }else if(supplyList[[2]] == "Hellinger"){
        supplyList[[2]] <- "Power"
        div <- 4
        reg <- supplyList[[3]]
        param1 <- -1
    }else if(supplyList[[2]] == "Power"){
        div <- 4
        reg <- supplyList[[3]]
        param1 <- supplyList[[4]]
    }else{
        stop("Please supply a divergence")
    }


    g_xy = Hausdorff_Vec_Rcpp(Cxy, supply, f_x, reg, param1, param2, div, eps)
    f_xy = Hausdorff_Vec_Rcpp(t(Cxy), demand, g_y, reg, param1, param2, div, eps)


    res <- sum(supplyList[[1]] * (legendre_entropy(supplyList[[3]], -f_x, supplyList[[2]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_x, supplyList[[2]], param1, param2) -
                                      legendre_entropy(supplyList[[3]], -f_xy, supplyList[[2]], param1, param2) - eps * grad_legrende(supplyList[[3]], -f_xy,supplyList[[2]], param1, param2))) +
        sum(demandList[[1]] *(legendre_entropy(supplyList[[3]], -g_y, supplyList[[2]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_y,supplyList[[2]], param1, param2)-
                                  legendre_entropy(supplyList[[3]], -g_xy, supplyList[[2]], param1, param2) - eps * grad_legrende(supplyList[[3]], -g_xy, supplyList[[2]], param1, param2)))

   
    return(res)
}

