
#' The Scaling Algorithm
#' 
#' This functions uses to scaling algorithm to solve the regularized unbalanced optimal transport problem. 
#' 
#' The algorithm minimizes the regularized unbalanced optimal transport problem given by
#' 
#' \eqn{\min_{r \in \!R^{X \times Y}_+} \cdot F_1(P_X r) + \cdot F_2(P_Y r) + \varepsilon KL(r|K)}{min_r F1(P_X r) +  F_2(P_Y r) + epsilon * KL(r|K)}
#' 
#' where F_X and F_Y are divergence functions, \eqn{P_X r} and \eqn{P_Y r} the marginals of \eqn{r}, 
#' \eqn{KL(\cdot | \cdot)} the Kullback Leibner divergence and \eqn{K = \exp(-c(x,y)/eps)}{K = exp(-c(x,y)/eps)} the kernel associated with the cost matrix.
#' 
#' 
#' The algorithm uses a stabilization method called log-domain stabilization to handle 
#' small values for \eqn{\varepsilon}{epsilon} and epsilon scaling to accelerate the convergence 
#' for those small values.
#' 
#' Three divergence functions are available to use with this algorithm: 
#' \itemize{
#'   \item Kullback-Leibner divergence ("KL")
#'   \item Total variation ("TV")
#'   \item The divergence associated with the range constraint ("RG")
#' }
#' 
#' Instead of providing a cost matrix to the algorithm, it is possible to provide the support of the supply and demand measures and 
#' have the function compute the cost matrix. Four additional parameters can be provided to define the cost matrix to compute. 
#' 
#' 
#' \insertRef{Chizat2016}{unbalancedTransport}
#' 
#' \insertRef{Schmitzer2016}{unbalancedTransport}
#' 
#' 
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
#' 
#'
#' @export
#'
scalingAlgorithm <- function(supplyList, demandList,
                             maxIteration, epsVector, tol = 1e-10, method = "euclidean",
                             exp = 1, p = 2,  wfr = FALSE, costMatrix = NULL){
    
    
    if(is.null(costMatrix)){
        
        costMatrix <- costMatrix(supplyList[[length(supplyList)]], demandList[[length(demandList)]], method, exp, p, wfr)
        
    }
    

    # Using either Kullback-Leibler divergence or total variation
    if(supplyList[[2]] == "KL"){
        Div1 <- 1
    }else if(supplyList[[2]] == "TV"){
        Div1 <- 2
    }else if(supplyList[[2]] == "RG"){
        Div1 <- 3
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for supplyList")
    }

    if(demandList[[2]] == "KL"){
        Div2 <- 1
    }else if(demandList[[2]] == "TV"){
        Div2 <- 2
    }else if(demandList[[2]] == "RG"){
        Div2 <- 3
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for demandList")
    }

    supply <- supplyList[[1]]
    demand <- demandList[[1]]


    if(Div1 != 3){
        supplyReg <- supplyList[[3]]
        supplyAlpha <- 0
        supplyBeta <- 0

    }else{
        supplyReg <- 0
        supplyAlpha <- supplyList[[3]]
        supplyBeta <- supplyList[[4]]

        if(supplyAlpha < 0 || supplyBeta < supplyAlpha){

            stop("0 <= Alpha <= Beta")

        }
    }


    if(Div2 != 3){
        demandReg <- demandList[[3]]
        demandAlpha <- 0
        demandBeta <- 0

    }else{
        demandReg <- 0
        demandAlpha <- demandList[[3]]
        demandBeta <- demandList[[4]]
        if(demandAlpha < 0 || demandBeta < demandAlpha){

            stop("0 <= Alpha <= Beta")

        }
    }
    
    res <- StabilizedScaling_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                  supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                  Div2, maxIteration, epsVector, tol)

    cost = -epsVector[length(epsVector)]*(sum( (res$TransportPlan-1)*(supply %*% t(demand))))-
        sum(supply * legendre_entropy(supplyReg, res$dual_f, Div1, supplyAlpha, supplyBeta))-
        sum(demand*legendre_entropy(demandReg, res$dual_g, Div2, demandAlpha, demandBeta))
    print("cost:")
    print(cost)

    return(res)

}
