#' @title legendre Entropy
#' @param lambda num value
#' @param x num value
#' @param DivFun div fun value
#' @param param1 num val
#' @param param2 num val
#' @noRd
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
        
        return(exp(x/lambda))
        
    }else if(DivFun == "Power"){
        
        return((1-(x/(lambda * (1-param1))))**(param1-1))
        
    }else if(DivFun == "RG"){
        return(pmax(-param1*sign(x), param2 * sign(x)))
        
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
#' The transport costs are calculated using the Sinkhorn algorithm \code{\link[unbalancedTransport]{regularizedTransport}}.
#' 
#' @references
#' \insertRef{Sejourne2019}{unbalancedTransport}
#'
#'  
#' @param supplyList A list containing the supply measure and if the cost matrix is not provided the support of the supply distribution.
#' @param demandList A list containing the demand measure and if the cost matrix is not provided the support of the demand distribution.
#' 
#' @param DivList A list containing the information about the divergence function for the supply and demand measures. The first element is the
#' abbreviation of the divergence function. The following values are available: 
#' 
#' \itemize{
#'   \item "KL" for the Kullback-Leibner divergence \eqn{F = \lambda \cdot KL()}{F1 = lambda * KL()}
#'   \item "TV" for total variation divergence \eqn{F = \lambda \cdot TV()}{F1 = lambda * TV()}
#'   \item "RG" for the range constraint 
#'   \item "Power" for the Power divergence with conjugate exponent r: \eqn{F_1 = \lambda \cdot Power_{r}()}{F1 = lambda * Power_r()}
#'   \item "Hellinger" for the Power divergence with conjugate exponent r = -1: \eqn{F_1 = \lambda \cdot Power_{-1}()}{F1 = lambda * Power_-1()}
#'   \item "Berg" for the Power divergence with conjugate exponent r = 0: \eqn{F_1 = \lambda \cdot Power_0()}{F1 = lambda * Power_0()}
#' }
#' 
#' 
#' The next elements in DivList are the additional parameters for the divergence functions. For "KL", "TV", "Hellinger"
#' and "Berg" the only required parameter is the regularization parameter \eqn{\lambda}{lambda}. The "RG" divergence function requires two parameters
#' alpha and beta with \eqn{0 \leq \alpha \leq \beta}{0 <= alpha <= beta} that define the lower and upper bounds. Lastly, "Power" requires
#' two parameters as well. The first is the regularization parameter \eqn{\lambda}{lambda}, followed by the conjugate exponent \eqn{r}  
#' 
#' @param epsVector  A numeric value or vector holding the regularization parameter. If a vector is given, epsilon scaling will be performed.
#' @param iterMax (optional) The maximum number of algorithm iterations. The default value is 10000.
#' @param tol  (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-5.
#' @param exp (optional) The exponent that is applied to the cost matrix. Can be used to compute quadratic cost. The default value is 1.
#' @param p (optional) Parameter for the L_p cost function. Can be omitted if either "euclidean" or "maximum" is used. The default value is 2.
#' @param wfr (optional) Set to "TRUE" to calculate the cost matrix for the Wasserstein-Fisher-Rao distance
#' \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}. Default value is "FALSE".
#' @param Cxy (optional) A cost matrix for transport between the supply and demand distributions.
#' @param Cxx (optional) A cost matrix for transport between the supply and supply distributions.
#' @param Cyy (optional) A cost matrix for transport between the demand and demand distributions.
#' @return The Sinkhorn divergence.
#' 
#' 
#' 
#' @examples 
#' x <- c(0,1)
#'y <- c(0,1)
#'a <- c(1.5,0.0)
#'b <- c(0.0,1.0)
#'
#'
#'supplyList <- list(a,x)
#'demandList <- list(b,y)
#'
#'
#'eps <- 0.1
#'reg <- list("KL", 1.5 )
#'
#'sinkhorn_divergence(supplyList,demandList,reg, eps)
#'
#'
#'
#' @export
sinkhorn_divergence <- function(supplyList, demandList, DivList, epsVector, iterMax = 1000, tol = 1e-5,
                                exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){
    
    eps <- epsVector[length(epsVector)]
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy) | is.null(Cyy) | is.null(Cxx)){

        X <- supplyList[[lenSup]]
        Y <- demandList[[lenDem]]


        if(is.null(Cxy)){
            Cxy <- costMatrix(X, Y , exp, p, wfr)
        }

        if(is.null(Cxx)){
            Cxx <- costMatrix(X, X, exp, p, wfr)
        }

        if(is.null(Cyy)){
            Cyy <- costMatrix(Y, Y , exp, p, wfr)
        }

    }
    
    
    
    res_xy <- regularizedTransport(supplyList,
                                   demandList,
                                   DivList,
                                   DivList,
                                   epsVector,
                                   iterMax,
                                   tol,
                                   costMatrix = Cxy,
                                   duals = TRUE,
                                   algorithm = "sinkhorn")
    
    res_x <- regularizedTransport(supplyList,
                                  supplyList,
                                  DivList,
                                  DivList,
                                  epsVector,
                                  iterMax,
                                  tol,
                                  costMatrix = Cxx,
                                  duals = TRUE,
                                  algorithm = "sinkhorn")
    
    res_y <- regularizedTransport(demandList,
                                  demandList,
                                  DivList,
                                  DivList,
                                  epsVector,
                                  iterMax,
                                  tol,
                                  costMatrix = Cyy,
                                  duals = TRUE,
                                  algorithm = "sinkhorn")
    
    
    f_xy <- res_xy$dual_f
    g_xy <- res_xy$dual_g
    
    f_x1 <- res_x$dual_f
    f_x2 <- res_x$dual_g
    
    g_y1 <- res_y$dual_f
    g_y2 <- res_y$dual_g
    
    
    
    if(DivList[[1]] == "TV"){
        
        
        func <- sum(supplyList[[1]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) +
            sum(demandList[[1]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
        

        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        supxsup <- supplyList[[1]] %*% t(supplyList[[1]])
        demxdem <- demandList[[1]] %*% t(demandList[[1]])

        func <- func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)
        
        
    }else if(DivList[[1]] == "RG" ){
        
        param1 <- DivList[[2]]
        param2 <- DivList[[3]]
        
      
        supxdem <- supplyList[[1]] %*% t(demandList[[1]])
        supxsup <- supplyList[[1]] %*% t(supplyList[[1]])
        demxdem <- demandList[[1]] %*% t(demandList[[1]])
 
        func <- sum(supplyList[[1]] * (-legendre_entropy(0, -f_xy, DivList[[1]], param1, param2)-
                                          0.5 * (-legendre_entropy(0, -f_x1, DivList[[1]], param1, param2))-
                                          0.5 * (-legendre_entropy(0, -f_x2, DivList[[1]], param1, param2)))) +
                                        sum(demandList[[1]] * (-legendre_entropy(0, -g_xy, DivList[[1]], param1, param2) -
                                            0.5 * (-legendre_entropy(0, -g_y1, DivList[[1]], param1, param2)) -
                                            0.5 * (-legendre_entropy(0, -g_y2, DivList[[1]], param1, param2))))
                  
        func <- func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)
        
    }else{
        
        param1 <- 0
        
        if(DivList[[1]] == "Berg"){
            DivList[[1]] <- "Power"
            
            param1 <- 0
            
            
        }else if(DivList[[1]] == "Hellinger"){
            DivList[[1]] <- "Power"
            
            param1 <- -1
            
        }else if(DivList[[1]] == "Power"){
            
            param1 <- DivList[[3]]
        }
        

        outf_xy <- -legendre_entropy(DivList[[2]], -f_xy, DivList[[1]], param1) -
            0.5*eps*grad_legrende(DivList[[2]], -f_xy, DivList[[1]], param1)
        outf_xx <- -legendre_entropy(DivList[[2]], -f_x1, DivList[[1]], param1) -
            0.5*eps*grad_legrende(DivList[[2]], -f_x1, DivList[[1]], param1)
        outg_xy <- -legendre_entropy(DivList[[2]], -g_xy, DivList[[1]], param1) -
            0.5*eps*grad_legrende(DivList[[2]], -g_xy, DivList[[1]], param1)
        outg_yy <- -legendre_entropy(DivList[[2]], -g_y1, DivList[[1]], param1) -
            0.5*eps*grad_legrende(DivList[[2]], -g_y1, DivList[[1]], param1)
        
        func <- sum(supplyList[[1]]*(outf_xy-outf_xx)) + sum(demandList[[1]] * (outg_xy - outg_yy))
        
        return(func)
    }
    

}




#' @title Hausdorff Divergence
#' 
#' @description Calculating the Hausdorff divergence for regularized unbalanced optimal transport.
#' 
#' 
#' @details This function calculates the Hausdorff divergence using the dual optimal optimal potentials
#' calculated by the \code{\link[unbalancedTransport]{regularizedTransport}} function.
#' The Hausdorff divergence is defined as 
#' \eqn{H_\varepsilon(\alpha, \beta) = <\alpha - \beta, \nabla F_\varepsilon(\alpha) - 
#' \nabla F_\varepsilon(\beta)>}{H_eps(a,b) = <a-b, nabla(F_eps(a)) - nabla(F_eps(b)) >}  
#' with supply and demand distributions \eqn{\alpha}{a} and \eqn{\beta}{b} and
#' \eqn{F_\varepsilon(\alpha) = -\frac{1}{2} OT_\varepsilon(\alpha,\alpha) + 
#' \frac{\varepsilon}{2}(m(\alpha)^2)}{F_eps(a) = -0.5 OT_eps(a,a) + 0.5eps * m(a)^2}.
#' 
#' @references
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#' 
#' @param supplyList A list containing the supply measure and if the cost matrix is not provided the support of the supply distribution.
#' @param demandList A list containing the demand measure and if the cost matrix is not provided the support of the demand distribution.
#' 
#' @param DivList A list containing the information about the divergence function for the supply and demand measures. The first element is the
#' abbreviation of the divergence function. The following values are available: 
#' 
#' \itemize{
#'   \item "KL" for the Kullback-Leibner divergence \eqn{F = \lambda \cdot KL()}{F1 = lambda * KL()}
#'   \item "TV" for total variation divergence \eqn{F = \lambda \cdot TV()}{F1 = lambda * TV()}
#'   \item "RG" for the range constraint 
#'   \item "Power" for the Power divergence with conjugate exponent r: \eqn{F_1 = \lambda \cdot Power_{r}()}{F1 = lambda * Power_r()}
#'   \item "Hellinger" for the Power divergence with conjugate exponent r = -1: \eqn{F_1 = \lambda \cdot Power_{-1}()}{F1 = lambda * Power_-1()}
#'   \item "Berg" for the Power divergence with conjugate exponent r = 0: \eqn{F_1 = \lambda \cdot Power_0()}{F1 = lambda * Power_0()}
#' }
#' 
#' The next elements in DivList are the additional parameters for the divergence functions. For "KL", "TV", "Hellinger"
#' and "Berg" the only required parameter is the regularization parameter \eqn{\lambda}{lambda}. The "RG" divergence function requires two parameters
#' alpha and beta with \eqn{0 \leq \alpha \leq \beta}{0 <= alpha <= beta} that define the lower and upper bounds. Lastly, "Power" requires
#' two parameters as well. The first is the regularization parameter \eqn{\lambda}{lambda}, followed by the conjugate exponent \eqn{r}.  
#' 
#' 
#' @param epsVector A numeric value or vector holding the regularization parameter. If a vector is given, epsilon scaling will be performed.
#' @param iterMax (optional) The maximum number of algorithm iterations. The default value is 5000
#' @param tol  (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-5.
#' @param exp (optional) The exponent that is applied to the cost matrix. Can be used to compute quadratic cost. The default value is 1.
#' @param p (optional) Parameter for the L_p cost function. Can be omitted if either "euclidean" or "maximum" is used. The default value is 2.
#' @param wfr (optional) Set to "TRUE" to calculate the cost matrix for the Wasserstein-Fisher-Rao distance
#' \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}. Default value is "FALSE".
#' @param Cxy (optional) A cost matrix for transport between the supply and demand distributions.
#' @param Cxx (optional) A cost matrix for transport between the supply and supply distributions.
#' @param Cyy (optional) A cost matrix for transport between the demand and demand distributions.
#' @return The Hausdorff divergence.
#' 
#' 
#' @examples 
#' x <- c(0,1)
#' y <- c(0,1)
#' a <- c(1.5,0.0)
#' b <- c(0.0,1.0)
#' 
#' supplyList <- list(a,x)
#' demandList <- list(b,y)
#' 
#' eps <- 0.1
#' reg <- list("KL", 1.5 )
#' 
#' hausdorff_divergence(supplyList, demandList, reg, eps)
#' 
#' 
#' 
#'
#' @export
hausdorff_divergence <- function(supplyList, demandList, DivList, epsVector, iterMax = 1000, tol = 1e-3,
                                 exp = 1, p = 2,  wfr = FALSE,
                                 Cxy = NULL, Cxx = NULL, Cyy = NULL){
    
    
    eps <- epsVector[length(epsVector)]
    lenSup <- length(supplyList)
    lenDem <- length(demandList)

    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], exp, p, wfr)
    }

    if(is.null(Cxx)){
        Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], exp, p, wfr)
    }

    if(is.null(Cyy)){
        Cyy <- costMatrix(demandList[[lenDem]], demandList[[lenDem]], exp, p, wfr)
    }

    reg <- 0
    param1 <- 0
    param2 <- 0

    supply <- supplyList[[1]]
    demand <- demandList[[1]]


    res_x <- regularizedTransport(supplyList, supplyList, DivList,
                                  DivList, epsVector, iterMax, tol, costMatrix = Cxx,
                                  duals = TRUE, algorithm = "sinkhorn")
    
    res_y <- regularizedTransport(demandList, demandList, DivList,
                                  DivList, epsVector, iterMax, tol, costMatrix = Cyy,
                                  duals = TRUE, algorithm = "sinkhorn")

    f_x <- res_x$dual_f
    g_x <- res_x$dual_g
    g_y <- res_y$dual_g
    f_y <- res_y$dual_f

    f_x <- 0.5*(f_x+g_x)
    g_y <- 0.5*(f_y+g_y)


    if(DivList[[1]] == "KL"){
        div <- 1
        reg <- DivList[[2]]

    }else if(DivList[[1]] == "TV"){
        div <- 2
        reg <- DivList[[2]]
    }else if(DivList[[1]] == "RG"){
        div <- 3
        param1 <- DivList[[2]]
        param2 <- DivList[[3]]

        if(param1 < 0 || param2 < param1){
            stop("0 <= Alpha <= Beta")
        }

    }else if(DivList[[1]] == "Berg"){
        DivList[[1]] <- "Power"
        div <- 4
        reg <- DivList[[2]]
        param1 <- 0
    }else if(DivList[[1]] == "Hellinger"){
        DivList[[1]] <- "Power"
        div <- 4
        reg <- DivList[[2]]
        param1 <- -1
    }else if(DivList[[1]] == "Power"){
        div <- 4
        reg <- DivList[[2]]
        param1 <- DivList[[3]]
    }else{
        stop("Please supply a divergence")
    }


    g_xy = Hausdorff_Vec_Rcpp(Cxy, supply, f_x, reg, param1, param2, div, eps)
    f_xy = Hausdorff_Vec_Rcpp(t(Cxy), demand, g_y, reg, param1, param2, div, eps)
    
   
    res <- sum(supplyList[[1]] * (legendre_entropy(reg, -f_x, DivList[[1]], param1,param2) +
                                      eps * grad_legrende(reg, -f_x, DivList[[1]], param1,param2)-
                                      legendre_entropy(reg, -f_xy, DivList[[1]], param1, param2) -
                                      eps * grad_legrende(reg, -f_xy,DivList[[1]], param1, param2))) +
        sum(demandList[[1]] *(legendre_entropy(reg, -g_y, DivList[[1]], param1, param2) +
                                  eps * grad_legrende(reg, -g_y,DivList[[1]], param1, param2)-
                                  legendre_entropy(reg, -g_xy, DivList[[1]], param1, param2) -
                                  eps * grad_legrende(reg, -g_xy, DivList[[1]], param1, param2)))

    return(res)
}

