

#' @title The Sinkhorn Algorihm for unbalanced optimal transport
#'
#' @param costMatrix A numeric matrix.
#' @param supplyList A supply list containing the divergence to use (either "KL" or "TV"),
#'  a numeric supply vector, the reference measure as numeric vector and the
#'  value for the lambda parameter.
#' @param demandList A demand list containing the divergence to use (either "KL" or "TV"),
#'  a numeric demand vector, the reference measure as numeric vector and the
#'  value for the lambda parameter.
#' @param maxIteration A numeric value for the maximum number of iterations.
#' The default value is 20000.
#' @param eps vector of epsilon values to use
#' @param tol tolerance for break
#'
#' @export
#'
sinkhornAlgorithmFromCost <- function(costMatrix, supplyList, demandList,
                                     maxIteration, eps, tol = 1e-8){
    
    supplyReg <- 0
    demandReg <- 0
    
    supplyAlpha <- 0
    supplyBeta <- 0
    demandAlpha <- 0
    demandBeta <- 0
    

    if(supplyList[[1]] == "KL"){
        Div1 <- 1
        supplyReg <- supplyList[[3]]
        
    }else if(supplyList[[1]] == "TV"){
        Div1 <- 2
        supplyReg <- supplyList[[3]]
    }else if(supplyList[[1]] == "RG"){
        Div1 <- 3
        supplyAlpha <- supplyList[[3]]
        supplyBeta <- supplyList[[4]]
        
        if(supplyAlpha < 0 || supplyBeta < supplyAlpha){
            stop("0 <= Alpha <= Beta")
        }
            
    }else if(supplyList[[1]] == "Berg"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- 0
    }else if(supplyList[[1]] == "Hellinger"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- -1
    }else if(supplyList[[1]] == "Power"){
        Div1 <- 4
        supplyReg <- supplyList[[3]]
        supplyAlpha <- supplyList[[4]]
    }else{
        stop("Please supply a divergence")
    }
    
    
    
    
    if(demandList[[1]] == "KL"){
        Div2 <- 1
        demandReg <- demandList[[3]]
    }else if(demandList[[1]] == "TV"){
        Div2 <- 2
        demandReg <- demandList[[3]]
    }else if(demandList[[1]] == "RG"){
        Div2 <- 3
        demandAlpha <- demandList[[3]]
        demandBeta <- demandList[[4]]
        
        if(demandAlpha < 0 || demandBeta < demandAlpha){
            stop("0 <= Alpha <= Beta")
        }
    }else if (demandList[[1]] == "Berg"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- 0
    }else if(demandList[[1]] == "Hellinger"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- -1
    }else if(demandList[[1]] == "Power"){
        Div2 <- 4
        demandReg <- demandList[[3]]
        demandAlpha <- demandList[[4]]
    }else{
        stop("Please supply a divergence")
    }
    
    supply <- supplyList[[2]]
    demand <- demandList[[2]]
    
    
    print(supplyAlpha)
    print(demandAlpha)




    res <- Sinkhorn_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                   supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                   Div2, maxIteration, eps, tol)
    
    
    
    
    TransportPlan <- res$TransportPlan*(supply %*% t(demand))
    
    returnList <- list("TransportPlan" = TransportPlan, "dual_f" = res$dual_f, "dual_g" =  res$dual_g)
    
    
    # res2 <- Sinkhorn_Eigen_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
    #                             supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
    #                             Div2, maxIteration, eps)
    # 
    # TransportPlan2 <- res2$TransportPlan*supply %*% t(demand)
    # 
    # returnList2 <- list("TransportPlan" = TransportPlan2, "dual_f" = res2$dual_f, "dual_g" =  res2$dual_g)
    
    
    return(returnList)
    
}



#' The log-domain stabilized Scaling Algorithm
#'
#' @param supplyList A supply list containing the divergence to use (either "KL" or "TV"),
#'  a numeric supply vector, the reference measure as numeric vector and the
#'  value for the lambda parameter.
#' @param demandList A demand list containing the divergence to use (either "KL" or "TV"),
#'  a numeric demand vector, the reference measure as numeric vector and the
#'  value for the lambda parameter.
#' @param maxIteration A numeric value for the maximum number of iterations.
#' The default value is 20000.
#' @param eps vector of epsilon values to use
#' @param method distance method
#' @param exp exponent for distance
#' @param p p norm value
#' @param wfr wfr distance?
#'
#' @export
#'
sinkhornAlgorithm <- function(supplyList, demandList, maxIteration, eps, method = "euclidean", exp = 1, p = 2,  wfr = FALSE){
    
    costMatrix <- costMatrix(supplyList[[length(supplyList)]], demandList[[length(demandList)]], method, exp, wfr, p)
    
    res <- scalingAlgorithmFromCost(costMatrix, supplyList, demandList, maxIteration, eps)
    
    return(res)
    
    
}





