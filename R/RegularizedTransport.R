
#' @title Regularized Unbalanced Optimal Transport
#' 
#' Solving regularized unbalanced optimal transport problems.
#' 
#' This functions solves the regularized unbalanced optimal transport problem using either the Scaling 
#' or Sinkhorn algorithm described in the papers given below.
#' The algorithms minimize the regularized unbalanced optimal transport problem given by
#' 
#' \eqn{\min_{r \in \!R^{X \times Y}_+} \cdot F_1(P_X r) + \cdot F_2(P_Y r) + \varepsilon KL(r|K)}{min_r F1(P_X r) +  F_2(P_Y r) + epsilon * KL(r|K)}
#' 
#' where F_X and F_Y are divergence functions, \eqn{P_X r} and \eqn{P_Y r} the marginals of \eqn{r}, 
#' \eqn{KL(\cdot | \cdot)} the Kullback Leibner divergence and \eqn{K = \exp(-c(x,y)/eps)}{K = exp(-c(x,y)/eps)} the kernel associated with the cost matrix.
#' 
#' 
#' 
#' The following divergence functions are currently implemented: 
#' \itemize{
#'   \item Kullback-Leibner divergence ("KL")
#'   \item Total variation ("TV")
#'   \item The divergence associated with the range constraint ("RG")
#'   \item The Power divergence ("Power") and its special cases the Hellinger ("Hellinger") and Berg ("Berg") divergence. (Only for the Sinkhorn algorithm.)
#' }
#' 
#' Instead of providing a cost matrix to the function, it is possible to provide the support of the supply and demand measures and 
#' have the function compute the cost matrix itself. 
#' 
#' 
#' \insertRef{Chizat2016}{unbalancedTransport}
#' \insertRef{Schmitzer2016}{unbalancedTransport}
#' \insertRef{Sejourne2019}{unbalancedTransport}
#' 
#'
#' @param supplyList A list containing the information about the supply distribution. 
#' 
#' 
#' divergence function and parameter. The first element is the supply distribution itself,
#' the second the abbreviation of the divergence function ("KL" for Kullback Leibner, 
#' "TV" for total variation and "RG" for range constraint) followed by the parameters 
#' needed for the divergence function: 
#' 
#' In case of the KL (\eqn{F_1 = \lambda \cdot KL()}{F1 = lambda * KL()})
#' and TV divergence (\eqn{F_1 = \lambda \cdot TV()}{F1 = lambda * TV()})
#' regularization parameter \eqn{\lambda}{lambda} has to be provided.
#' 
#' 
#' The divergence function associated with the range constraint needs two parameters, that define the upper and 
#' lower bound. 
#' 
#' If the cost matrix is not provided, the support of the supply distribution has to be provided as last element in the list. This can
#' be omitted if a cost matrix is passed as argument.
#' 
#' 
#' @param supplyList a
#' @param demandList a
#' @param supplyDivList a
#' @param demandDivList a
#' @param epsVector a
#' @param maxIteration a
#' @param tol a
#' @param p a
#' @param q a
#' @param wfr a
#' @param costMatrix a
#' @param duals a
#' @param algorithm a
#'
#' @export
regularizedTransport <- function(supplyList, demandList, supplyDivList = list("KL", 0.5),
                                    demandDivList = list("KL", 0.5), epsVector = c(0.5), maxIteration = 5000,
                                    tol = 1e-8, p = 1, q = 2, wfr = FALSE, costMatrix = NULL, duals = FALSE,
                                    algorithm = "scaling"){
    
    

    
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
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], p, q, wfr)
            res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                 supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                 Div2, maxIteration, epsVector, tol)
            
        }else{
            
            # res <- Sinkhorn_SpaceOptim_Rcpp(supplyList[[2]], demandList[[2]], supply, demand, supplyReg, supplyAlpha,
            #                      supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
            #                      Div2, maxIteration, epsVector, tol)
            # 
            # 
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], p, q, wfr)
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
            
            costMatrix <- costMatrix(supplyList[[2]], demandList[[2]], p, q, wfr)
            
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
