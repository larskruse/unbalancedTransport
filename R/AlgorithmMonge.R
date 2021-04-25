#' @importFrom stats runif





#' @title Monge Matrix Generator
#'
#' Generates a Monge matrix of given dimensions.
#'
#'
#' @param dimx number of columns
#' @param dimy number of rows
#'
#' @return A Monge matrix of size dimx times dimy.
#' @export
#'
genMongeMat <- function(dimx, dimy){
    a <- round(runif(dimx, min = 0, max = 1),3)
    b <- round(runif(dimy, min = 0, max = 1),3)
    a <- sort(a)
    b <- sort(b)


    mat <- matrix(rep(0, dimy*dimx), nrow = dimx)

    for(i in 1:dimx){
        for(j in 1:dimy){
            mat[i,j] <- abs(a[i]-b[j])^2


        }

    }

    return(mat)

}



#' checkMongeProperty
#'
#' Checking if a given matrix fulfills the Monge property
#'
#' @param mat A numeric matrix.
#'
#' @return '0' if the matrix does not fulfill the Monge property, '1' if it does
#' @export
#'
checkMongeProperty <- function(mat){
    dim <- dim(mat)

    for(i in 1:(dim[1]-1)){
        for(j in 1:(dim[2]-1)){

            if(mat[i,j]+mat[i+1,j+1] > mat[i+1,j]+mat[i,j+1]){

                print("Matrix does not fulfill the Monge property")
                print("Failed at:")
                cat("\n i Coordiante: " + i)
                cat("\n j Coordiante: " + j )
                return(0)

            }

        }

    }

    print("Matrix fulfills the Monge property")
    return(1)

}


#' pathToPlan
#'
#' Computes the transport plan from a given transport path.
#'
#' @param iList A numeric vector of i indexes
#' @param jList A numeric vector of j indexes
#' @param wList A numeric vector of weights
#' @param dimX  The number of rows in the transport matrix
#' @param dimY The number of columns in the transport matrix
#'
#' @return A dimX times dimY transport matrix
#' @export
#'
pathToPlan <- function(iList, jList, wList, dimX, dimY){
    plan <- matrix(0L, nrow = dimX, ncol = dimY)

    for(i in 1:length(iList)){
        plan[iList[i], jList[i]] <- wList[i]

    }

    return(plan)

}



#' Monge Unbalanced Optimal Transport Algorithm
#'
#' This calculates the optimal transport cost and optimal transport path of an
#' unbalanced optimal transport problem. The cost matrix has to fulfill the Monge
#' property.
#'
#' @param costMatrix A numeric matrix fulfilling the Monge property
#' @param supply A numeric supply vector
#' @param demand A numeric demand vector
#' @param constructionCost A numeric construction cost value
#' @param destructionCost A numeric destruction cost value
#'
#' @return A list containing the optimal transport cost, plan, import vector
#' and export vector
#' @export
#'
#' @examples
#' # generating random data
#' mongeMatrix <- genMongeMat(30,30)
#'
#' supply <- runif(30,0,1)
#' demand <- runif(30,0,1)
#'
#' createCost <- runif(1,0,0.2)
#' destructCost <- runif(1,0,0.3)
#'
#' # calculating the optimal transport cost and plan
#' res <- mongeTransport(mongeMatrix, supply, demand, createCost, destructCost)
#'
#' # plotting the transport plan
#' plotUOTP(res$transportPlan, res$import, res$export)
#'

mongeTransport <- function(costMatrix, supply, demand, constructionCost, destructionCost){

    # computing the optimal transport cost, path and import / export vectors
    transportResult <- Monge_Rcpp(costMatrix, supply, demand, constructionCost, destructionCost)

    # computing the transport plan from the path
    transportPlan <- pathToPlan(transportResult$iList, transportResult$jList, transportResult$weightList, length(supply), length(demand))

    transport <- list(cost = transportResult$transportCost, transportPlan = transportPlan, import = transportResult$import, export = transportResult$export)
    return(transport)
}


