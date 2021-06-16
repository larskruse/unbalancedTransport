

#' The Sinkhorn Algorithm
#' 
#' Solving regularized unbalanced optimal transport using the Sinkhorn algorithm.
#' 
#' This functions uses to Sinkhorn algorithm to solve the regularized unbalanced optimal transport problem. 
#' The algorithm minimizes the optimization problem
#' 
#' \eqn{\min_{r \in \!R^{X \times Y}_+} \cdot F_1(P_X r) + \cdot F_2(P_Y r) + \varepsilon KL(r|K)}{min_r F1(P_X r) +  F_2(P_Y r) + epsilon * KL(r|K)}
#' 
#' where F_X and F_Y are divergence functions, \eqn{P_X r} and \eqn{P_Y r} the marginals of \eqn{r}, 
#' \eqn{KL(\cdot | \cdot)} the Kullback Leibner divergence and \eqn{K = \exp(-c(x,y)/eps)}{K = exp(-c(x,y)/eps)} the kernel associated with the cost matrix.
#' 
#' 
#' Four divergence functions are available to use with this algorithm: 
#' \itemize{
#'   \item Kullback-Leibner divergence ("KL")
#'   \item Total variation ("TV")
#'   \item The divergence associated with the range constraint ("RG")
#'   \item The Power divergence ("Power")
#' }
#' 
#' Two special cases of the Power divergence are directly available: 
#' The Berg divergence ("Berg") uses the conjugate exponent \eqn{r = -1} and the Hellinger divergence ("Hellinger") with \eqn{r = 0}.
#'  
#' 
#' Instead of providing a cost matrix to the algorithm, it is possible to provide the support of the supply and demand measures and 
#' have the function compute the cost matrix automatically.
#' 
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
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
#' @param demandList A list containing the information about the demand measure. It has to have the same structure as the supplyList.
#' 
#' @param maxIteration (optional) The maximum number of iterations.
#' @param epsVector A numeric vector for the regularization parameter. In order to use epsilon scaling, a declining vector of values can be used.
#' @param tol (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value. The algorithm
#' is terminated as it is already converged. The default value is 1e-8.
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
#' @param costMatrix (optional) Instead of having the algorithm compute the cost matrix, a custom cost matrix can be passed to the algorithm. 
#' @param duals (optional) Set to TRUE to return the dual potentials as well.
#'
#' @return A list containing the transport plan ("TransportPlan") , the transport cost ("TransportCost"),
#'  the maximum change of the dual potential in the last iteration ("converge") and if "duals" is set to TRUE, the dual potentials.
#'  
#' @examples 
#'
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
#' maxIter <- 2000
#' eps <- 1e-3
#' 
#' suppList <- list(p, "KL", 0.04, X)
#' demList <- list(q, "KL", 0.04, Y)
#' res <- sinkhornAlgorithm(suppList, demList, eps, maxIter,  exp = 2)
#' plot1DTransport(res$TransportPlan, supply, demand)
#' plotGridTransport(res$TransportPlan)
#'
#'
#' @export
#'
sinkhornAlgorithm <- function(supplyList, demandList, epsVector, maxIteration = 2000, tol = 1e-8, method = "euclidean",
                                     exp = 1, p = 2,  wfr = FALSE, costMatrix = NULL, duals = FALSE){
    
    if(is.null(costMatrix)){
        
        costMatrix <- costMatrix(supplyList[[length(supplyList)]], demandList[[length(demandList)]], method, exp, p, wfr)
        
    }
    

    
    
    supplyReg <- 0
    demandReg <- 0
    
    supplyAlpha <- 0
    supplyBeta <- 0
    demandAlpha <- 0
    demandBeta <- 0
    

    if(supplyList[[2]] == "KL"){
        Div1 <- 1
        supplyReg <- supplyList[[3]]
        
    }else if(supplyList[[2]] == "TV"){
        Div1 <- 2
        supplyReg <- supplyList[[3]]
    }else if(supplyList[[2]] == "RG"){
        Div1 <- 3
        supplyAlpha <- supplyList[[3]]
        supplyBeta <- supplyList[[4]]
        
        if(supplyAlpha < 0 || supplyBeta < supplyAlpha){
            stop("0 <= Alpha <= Beta")
        }
            
    }else if(supplyList[[2]] == "Berg"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- 0
    }else if(supplyList[[2]] == "Hellinger"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- -1
    }else if(supplyList[[2]] == "Power"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- supplyList[[4]]
    }else{
        stop("Please supply a divergence")
    }
    
    
    
    
    if(demandList[[2]] == "KL"){
        Div2 <- 1
        demandReg <- demandList[[3]]
    }else if(demandList[[2]] == "TV"){
        Div2 <- 2
        demandReg <- demandList[[3]]
    }else if(demandList[[2]] == "RG"){
        Div2 <- 3
        demandAlpha <- demandList[[3]]
        demandBeta <- demandList[[4]]
        
        if(demandAlpha < 0 || demandBeta < demandAlpha){
            stop("0 <= Alpha <= Beta")
        }
    }else if (demandList[[2]] == "Berg"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- 0
    }else if(demandList[[2]] == "Hellinger"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- -1
    }else if(demandList[[2]] == "Power"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- demandList[[4]]
    }else{
        stop("Please supply a divergence")
    }
    
    supply <- supplyList[[1]]
    demand <- demandList[[1]]


    res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                   supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                   Div2, maxIteration, epsVector, tol)
    
    
    
    
    TransportPlan <- res$TransportPlan*(supply %*% t(demand))
    

    
    # calcuated transport cost
    
    if(Div1 == Div2 & supplyReg == demandReg & supplyAlpha == demandAlpha & supplyBeta == demandBeta){
        
        cost <- regularized_ot_intern(supplyList, demandList, res$dual_f, res$dual_g, epsVector[length(epsVector)], costMatrix)
    
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
    
    
}
