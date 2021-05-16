
#' Extendes the UOT Problemt to the balanced case and solves it using the known methods
#'
#' @param supplyList list of supply values
#' @param demandList list of demand values
#' @param costMethod method to use for the cost matrix
#' @param costExp exponent for the cost matrix
#' @param costP parameter to use for the minkowsi cost matrix
#' @param wfr using wfr cost matrix??
#'
#' @export
#'
BalancedExtensionSolver <- function(supplyList, demandList, costMethod, costExp, costP, wfr = FALSE){

    costMatrix <- costMatrix(supplyList[[1]], demandList[[1]], costMethod, costExp, wfr, costP)

    supply <- c(supplyList[[2]], sum(demandList[[2]]))
    demand <- c(demandList[[2]], sum(supplyList[[2]]))

    costMatrix <- cbind(costMatrix, supplyList[[3]])
    costMatrix <- rbind(costMatrix, c(demandList[[3]],0))



    res <- transport::transport(supply, demand, costMatrix, method = "revsimplex")

    print(res)

    cost <- 0
    for (i in 1:nrow(res)){
        cost = cost + res[i,3]*costMatrix[res[i,1], res[i,2]]
    }

    import <-rep(0, length(demandList[[2]]))
    export <- rep(0, length(supplyList[[2]]))

    expTransport <- res[res$to > length(demandList[[2]]) & res$from <= length(supplyList[[2]]) ,]
    impTransport <- res[res$from > length(supplyList[[2]]) & res$to <= length(demandList[[2]]),]



    import[impTransport$to] <- impTransport$mass
    export[expTransport$from] <- expTransport$mass


    res <- res[res$to <= length(demandList[[2]]) & res$from <= length(supplyList[[2]]),]



    if(length(res$from) > 0){

        transportPlan <- matrix(0,length(supplyList[[2]]),length(demandList[[2]]))
        transportPlan[cbind(res$from,res$to)] <- res$mass

    }else{
        transportPlan <- matrix(0,length(supply),length(demand))

    }


    return(list(cost, transportPlan, import, export))

}




