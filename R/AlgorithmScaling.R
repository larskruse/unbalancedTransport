
#' The log-domain stabilized Scaling Algorithm
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
#' @param epsVector vector of epsilon values to use
#'
#' @export
#'
scalingAlgorithm <- function(costMatrix, supplyList, demandList,
                             maxIteration, epsVector){

    # Using either Kullback-Leibler divergence or total variation
    if(supplyList[[1]] == "KL"){
        Div1 <- 1
    }else if(supplyList[[1]] == "TV"){
        Div1 <- 2
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for supplyList")
    }

    if(demandList[[1]] == "KL"){
        Div2 <- 1
    }else if(demandList[[1]] == "TV"){
        Div2 <- 2
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for demandList")
    }

    supply <- supplyList[[2]]
    demand <- demandList[[2]]

    supplyRefMeasure <- supplyList[[3]]
    demandRefMeasure <- demandList[[3]]

    supplyReg <- supplyList[[4]]
    demandReg <- demandList[[4]]

    res <- StabilizedScaling_Rcpp(costMatrix, supply, demand, supplyRefMeasure,
                                demandRefMeasure, supplyReg, demandReg, Div1,
                                Div2, maxIteration, epsVector)

    transportPlan <- res$TransportMap

    transport <- list(transportPlan = transportPlan)
    return(transport)

}

