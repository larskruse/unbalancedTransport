
#' @title Regularized Unbalanced Optimal Transport
#' 
#' @description  Solving regularized unbalanced optimal transport problems.
#' 
#' 
#' @details  This functions solves the regularized unbalanced optimal transport problem using either the Scaling 
#' or Sinkhorn algorithm described in the papers given below.
#' The algorithms minimize the regularized unbalanced optimal transport problem given by
#' 
#' \eqn{\min_{b \in \!R^{X \times Y}_+} \cdot F_1(P_X b) + \cdot F_2(P_Y b) + \varepsilon KL(b|K)}{min_b F1(P_X b) +  F_2(P_Y b) + epsilon * KL(b|K)}
#' 
#' where F_X and F_Y are divergence functions, \eqn{P_X b} and \eqn{P_Y b} the marginals of \eqn{b}, 
#' \eqn{KL(\cdot | \cdot)}{KL( | )} the Kullback Leibner divergence and \eqn{K = \exp(-c(x,y)/eps)}{K = exp(-c(x,y)/eps)} the kernel associated with the cost matrix.
#' 
#' 
#' 
#' The following divergence functions are currently implemented: 
#' \itemize{
#'   \item Kullback-Leibner divergence ("KL")
#'   \item Total variation ("TV")
#'   \item The divergence associated with the range constraint ("RG")
#'   \item The Power divergence ("Power") and its special cases the Hellinger ("Hellinger") and Berg ("Berg") divergence. (Only available for the Sinkhorn algorithm.)
#' }
#' 
#' Instead of providing a cost matrix to the function, it is possible to provide the support of the supply and demand measures and 
#' have the function compute the cost matrix itself. 
#' 
#' 
#' @references
#' \insertRef{Chizat2016}{unbalancedTransport}
#' 
#' \insertRef{Schmitzer2016}{unbalancedTransport}
#' 
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#' 
#'
#' @param supplyList A list containing the supply measure and, if the cost matrix is not provided, the support of the supply distribution.
#' @param demandList A list containing the demand measure and, if the cost matrix is not provided, the support of the demand distribution.
#' 
#' @param supplyDivList A list containing the information about the divergence function for the supply measure. The first element is the
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
#' The Power divergence and its special cases are only implemented for the Sinkhorn algorithm. 
#' 
#' The next elements in supplyDivList are the additional parameters for the divergence functions. For "KL", "TV", "Hellinger"
#' and "Berg" the only required parameter is the regularization parameter \eqn{\lambda}{lambda}. The "RG" divergence function requires two parameters
#' alpha and beta with \eqn{0 \leq \alpha \leq \beta}{0 <= alpha <= beta} that define the lower and upper bounds. Lastly, "Power" requires
#' two parameters as well. The first is the regularization parameter \eqn{\lambda}{lambda}, followed by the conjugate exponent \eqn{r}.  
#' The default value is list("KL",0.5).
#' 
#' @param demandDivList A list containing the information about the divergence
#' function for the demand measure in the same form as the supplyDivList. The default value is list("KL",0.5).
#' @param epsVector A numeric value or vector holding the regularization parameter. If a vector is given, epsilon scaling will be performed.
#' @param maxIteration (optional) The maximum number of algorithm iterations. The default value is 5000
#' @param tol (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value the algorithm is 
#' considered to be converged. The default value is 1e-5.  
#' @param exp (optional) The exponent applied to the cost function. Default value is 1.
#' @param p (optional) The parameter for calculating the L_q cost. Can be a positive real number or Inf. Default value is 2 calculating the euclidean distance.
#' @param wfr (optional) Set to "TRUE" to calculate the cost matrix for the Wasserstein-Fisher-Rao distance
#' \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}. Default value is "FALSE".
#' @param costMatrix (optinal) A cost matrix for transport between the supply and demand points.
#' @param duals (optinal) Set to "TRUE" to return the value of the dual solutions as well. Default value is "FALSE"
#' @param algorithm (optinal) Set to "sinkhorn" to use the Sinkhorn algorithm. Default value is "scaling". Only 
#' these two algorithm are available.
#' 
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
#' maxIter <- 1000
#' eps <- 1e-3
#' 
#' 
#' suppyDiv <- list("KL", 0.04)
#' demandDiv <- list("KL", 0.04)
#' res <- regularizedTransport(supply, demand, suppyDiv, demandDiv,
#'                             maxIteration = maxIter, epsVector = eps, exp = 2)
#' plot1DTransport(res$TransportPlan, supply, demand)
#'
#'
#' @export
regularizedTransport <- function(supplyList, demandList, supplyDivList, demandDivList, epsVector,
                                 maxIteration = 5000, tol = 1e-5, exp = 1, p = 2, wfr = FALSE,
                                 costMatrix = NULL, duals = FALSE, algorithm = "scaling"){
    
    

    
    supplyReg <- 0
    demandReg <- 0
    
    supplyAlpha <- 0
    supplyBeta <- 0
    demandAlpha <- 0
    demandBeta <- 0
    
    
    if(supplyDivList[[1]] == "KL"){
        Div1 <- 1
        supplyReg <- supplyDivList[[2]]
        
    }else if(supplyDivList[[1]] == "TV"){
        Div1 <- 2
        supplyReg <- supplyDivList[[2]]
    }else if(supplyDivList[[1]] == "RG"){
        Div1 <- 3
        supplyAlpha <- supplyDivList[[2]]
        supplyBeta <- supplyDivList[[3]]
        
        if(supplyAlpha < 0 || supplyBeta < supplyAlpha){
            stop("0 <= Alpha <= Beta")
        }
        
    }else if(supplyDivList[[1]] == "Berg"){
        Div1 <- 4
        supplyReg <- supplyDivList[[2]]
        supplyAlpha <- 0
    }else if(supplyDivList[[1]] == "Hellinger"){
        Div1 <- 4
        supplyReg <- supplyDivList[[2]]
        supplyAlpha <- -1
    }else if(supplyDivList[[1]] == "Power"){
        Div1 <- 4
        supplyReg <- supplyDivList[[2]]
        supplyAlpha <- supplyDivList[[3]]
    }else{
        stop("Please supply a divergence")
    }
    
    
    
    
    if(demandDivList[[1]] == "KL"){
        Div2 <- 1
        demandReg <- demandDivList[[2]]
    }else if(demandDivList[[1]] == "TV"){
        Div2 <- 2
        demandReg <- demandDivList[[2]]
    }else if(demandDivList[[1]] == "RG"){
        Div2 <- 3
        demandAlpha <- demandDivList[[2]]
        demandBeta <- demandDivList[[3]]
        
        if(demandAlpha < 0 || demandBeta < demandAlpha){
            stop("0 <= Alpha <= Beta")
        }
    }else if (demandDivList[[1]] == "Berg"){
        Div2 <- 4
        demandReg <- demandDivList[[2]]
        demandAlpha <- 0
    }else if(demandDivList[[2]] == "Hellinger"){
        Div2 <- 4
        demandReg <- demandDivList[[2]]
        demandAlpha <- -1
    }else if(demandDivList[[1]] == "Power"){
        Div2 <- 4
        demandReg <- demandDivList[[2]]
        demandAlpha <- demandDivList[[3]]
    }else{
        stop("Please supply a divergence")
    }
    
    
    if((Div1 > 3 | Div2 > 3) & algorithm == "scaling"){
        
        print("These divergence functions are not implemeted for the Scaling algorithm.")
        print("The Sinkhorn algorithm is used instead of the Scaling algorithm.")
        algorithm <- "sinkhorn"
        
    }
    
    supply <- supplyList[[1]]
    demand <- demandList[[1]]
    
    
    if(algorithm == "sinkhorn"){
        
        if(!is.null(costMatrix)){
            
            res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                 supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                 Div2, maxIteration, epsVector, tol)
            
            
        }else if((is.null(costMatrix) & length(supply) <= 1000)){
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], exp, p, wfr)
            res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                 supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                 Div2, maxIteration, epsVector, tol)
            
        }else{
            
            # res <- Sinkhorn_SpaceOptim_Rcpp(supplyList[[2]], demandList[[2]], supply, demand, supplyReg, supplyAlpha,
            #                      supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
            #                      Div2, maxIteration, epsVector, tol)
            # 
            # 
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], exp, p, wfr)
            res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                 supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                 Div2, maxIteration, epsVector, tol)
            
        }
        
        TransportPlan <- res$TransportPlan*(supply %*% t(demand))
        
        
        if(Div1 == Div2 & supplyReg == demandReg & supplyAlpha == demandAlpha & supplyBeta == demandBeta){
            
            cost <- regularized_ot_intern(supplyList, demandList, supplyDivList, demandDivList, res$dual_f, res$dual_g, epsVector[length(epsVector)], costMatrix)
            
        }else{
            
            cost = -epsVector[length(epsVector)]*(sum((res$TransportPlan-1)*(supply %*% t(demand))))-
                sum(supply * legendre_entropy(supplyReg, -res$dual_f, supplyList[[2]], supplyAlpha, supplyBeta))-
                sum(demand * legendre_entropy(demandReg, -res$dual_g, demandList[[2]], demandAlpha, demandBeta))
            
            
        }
        
        
        
        if(duals){
            returnList <- list("TransportPlan" = TransportPlan, "TransportCost" = cost,
                               "dual_f" = res$dual_f, "dual_g" = res$dual_g, "converge" = res$converge)
            return(returnList)
        }else{
            returnList <- list("TransportPlan" = TransportPlan, "TransportCost" = cost, "converge" = res$converge)
            return(returnList)
        }
    
        
    }else if(algorithm == "scaling"){
        
        if(is.null(costMatrix)){
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], exp, p, wfr)
            
        }
        
        res <- StabilizedScaling_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                             supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                             Div2, maxIteration, epsVector, tol)
        
        TransportPlan <- res$TransportPlan
        
        
        if(Div1 == Div2 & supplyReg == demandReg & supplyAlpha == demandAlpha & supplyBeta == demandBeta){
            
            cost <- regularized_ot_intern(supplyList, demandList, supplyDivList, demandDivList, res$dual_f, res$dual_g, epsVector[length(epsVector)], costMatrix)
            
        }else{
            
            outf_xy <- legendre_entropy(supplyReg, -res$dual_f, supplyList[[2]], supplyAlpha, supplyBeta)
            outg_xy <- legendre_entropy(demandReg, -res$dual_g, demandList[[2]], demandAlpha, demandBeta)
            outf_xy[!is.finite(outf_xy) & supplyList[[1]] == 0] <- 0
            outg_xy[!is.finite(outg_xy) & demandList[[1]] == 0] <- 0
            
            
            cost = -epsVector[length(epsVector)]*(sum((TransportPlan-1)*(supply %*% t(demand))))-
                sum(supply * outf_xy)-
                sum(demand * outg_xy)
            
        }
        

        if(duals){
            returnList <- list("TransportPlan" = TransportPlan, "TransportCost" = cost,
                               "dual_f" = res$dual_f, "dual_g" = res$dual_g, "converge" = res$converge)
            return(returnList)
        }else{
            returnList <- list("TransportPlan" = TransportPlan, "TransportCost" = cost, "converge" = res$converge)
            return(returnList)
        }
        
    }
    
    
}
