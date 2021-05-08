
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
scalingAlgorithmFromCost <- function(costMatrix, supplyList, demandList,
                             maxIteration, epsVector){

    # Using either Kullback-Leibler divergence or total variation
    if(supplyList[[1]] == "KL"){
        Div1 <- 1
    }else if(supplyList[[1]] == "TV"){
        Div1 <- 2
    }else if(supplyList[[1]] == "RG"){
        Div1 <- 3
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for supplyList")
    }

    if(demandList[[1]] == "KL"){
        Div2 <- 1
    }else if(demandList[[1]] == "TV"){
        Div2 <- 2
    }else if(demandList[[1]] == "RG"){
        Div2 <- 3
    }else{
        stop("Please use 'KL' or 'TV' as divergence Function parameter for demandList")
    }

    supply <- supplyList[[2]]
    demand <- demandList[[2]]


    if(Div1 != 3){
        supplyReg <- supplyList[[3]]
        supplyAlpha <- 0
        supplyBeta <- 0

    }else{
        supplyReg <- 0
        supplyAlpha <- supplyList[[3]]
        supplyBeta <- supplyList[[4]]
    }


    if(Div2 != 3){
        demandReg <- demandList[[3]]
        demandAlpha <- 0
        demandBeta <- 0

    }else{
        demandReg <- 0
        demandAlpha <- demandList[[3]]
        demandBeta <- demandList[[4]]
    }

    print(supplyReg)
    print(demandReg)
    print(Div1)
    print(Div2)
    res <- StabilizedScaling_Rcpp(costMatrix, supply, demand, supplyReg, supplyAlpha,
                                  supplyBeta, demandReg, demandAlpha, demandBeta, Div1,
                                  Div2, maxIteration, epsVector)

    print(res$TransportMap[1:10,1:10])


    transportPlan <- res$TransportMap

    transport <- list(transportPlan = transportPlan)
    return(transport)

}



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
#' @param method distance method
#' @param exp exponent for distance
#' @param p p norm value
#' @param wfr wfr distance?
#'
#' @export
#'
scalingAlgorithm <- function(supplyList, demandList, maxIteration, epsVector, method = "euclidean", exp = 1, p = 2,  wfr = FALSE){

    costMatrix <- costMatrix(supplyList[[1]], demandList[[1]], method, exp, wfr, p)

    res <- scalingAlgorithmFromCost(costMatrix, supplyList, demandList, maxIteration, epsVector)


}





