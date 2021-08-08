
#' Coordinates by distance matrix
#'
#' Calculating the coordinates of the supply and demand points using
#' the eigenvalue decomposition of the distance matrix.
#'
#' @param c A numeric distance matrix.
#' @return The coordinates based on the distance matrix.
#' @noRd

positionsByEigenvalue <- function(c){
        dimC <- dim(c)[1]
        m <- matrix(rep(0, length(c)), dimC)
        for (i in 1:dimC){
                for (j in 1:dimC){
                        m[i,j] = (c[1,j]^2+c[i,1]^2-c[i,j]^2) / 2
                }
        }
        eig <- eigen(m)
        values <- eig$values
        values[abs(values) < 1e-5] <- 0

        sqrtValues <- sqrt(values)

        sqrtValues[values < 1e-5] <- 0


        coordinates <- eig$vectors %*% diag(sqrtValues)
        coordinates[abs(coordinates) < 1e-5] <- 0

        usedDimensions <- c()

        for(i in 1:dimC){
                if(!is.null(dim(coordinates)) & (any(coordinates[,i] != 0))){
                        usedDimensions <- append(usedDimensions, i)
                }
        }

        return(coordinates[,usedDimensions])
}

# 
# 
# #' Wasserstein-Fisher-Rao Cost Matrix
# #'
# #' Calculating the Wasserstein-Fisher-Rao cost matrix for the absolute value
# #' distance function between two 1D sets of points.
# #'
# #' @param X A numeric vector.
# #' @param Y A numeric vector.
# #' @return The WFR cost matrix.
# #' @noRd
# #'
# wfrCost <- function(X,Y){
#         C <- matrix(rep(0, length(X)*length(Y)), nrow = length(X))
# 
#         for(i in 1:length(X)){
#                 for(j in 1:length(Y)){
#                         C[i,j] <- -2*log(cos(pi*(0.5*min(1,abs(X[i]-Y[j])/0.2))))
#                         C[i,j] <- Re(C[i,j])
#                 }
#         }
#         return(C)
# }



#' Cost Matrix
#' @description Calculates the cost matrix between points.
#'
#' @param x A numeric matrix. Each row corresponds to the coordinates of one point in the first point cloud.
#' @param y A numeric matrix. Each row corresponds to the coordinates of one point in the second point cloud
#' @param p exponent to apply to the distance
#' @param wfr computing the wasserstein fisher rao distance
#' @param q parameter for the minkowski metric. standard p = 2 give the minkowski metric.
#' @return The distance matrix between the points. The rows correspond to the points in x, the columns to the
#' points in y
#' @noRd
costMatrix <- function(x, y, p = 1,q = 2, wfr = FALSE){

        if(is.null(ncol(x))){
                x <- matrix(x, ncol = 1)

        }

        if(is.null(ncol(y))){
                y <- matrix(y, ncol = 1)

        }

        if(ncol(x) != ncol(y)){
                print("Unequal dimensions.")
        }


        cMatrix <- matrix(rep(0, nrow(x)*nrow(y)), ncol = nrow(y))

        if(q == Inf){
            
            for(i in 1:nrow(x)){
                
                cMatrix[i,] <- max(abs(t(t(y) - x[i,])))^exp
                
            }

        }else{
            for (i in 1:nrow(x)){
                
                cMatrix[i,] <- ((rowSums((abs(t(t(y) - x[i,])))^q))^(1/q))^p
            }
                
        }

        if(wfr){

                cMatrix <- -2*log(cospi(pmin(cMatrix/pi, 1/2)))
        }

        return(cMatrix)

}

