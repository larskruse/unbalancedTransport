
#' Coordinates by distance matrix
#'
#' Calculating the coordinates of the supply and demand points using
#' the eigenvalue decomposition of the distance matrix.
#'
#' @param c A numeric distance matrix.
#' @return The coordinates based on the distance matrix.
#'

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
        values[abs(values) < 0.00001] <- 0

        sqrtValues <- sqrt(values)

        sqrtValues[values < 0.00001] <- 0


        coordinates <- eig$vectors %*% diag(sqrtValues)
        coordinates[abs(coordinates) < 0.00001] <- 0

        usedDimensions <- c()

        for(i in 1:dimC){
                if(!is.null(dim(coordinates)) & (any(coordinates[,i] != 0))){
                        usedDimensions <- append(usedDimensions, i)
                }
        }

        return(coordinates[,usedDimensions])
}



#' Wasserstein-Fisher-Rao Cost Matrix
#'
#' Calculating the Wasserstein-Fisher-Rao cost matrix for the absolute value
#' distance function between two 1D sets of points.
#'
#' @param X A numeric vector.
#' @param Y A numeric vector.
#' @return The WFR cost matrix.
#' @export
#'
wfrCost <- function(X,Y){
        C <- matrix(rep(0, length(X)*length(Y)), nrow = length(X))

        for(i in 1:length(X)){
                for(j in 1:length(Y)){
                        C[i,j] <- -2*log(cos(pi*(0.5*min(1,abs(X[i]-Y[j])/0.2))))
                        C[i,j] <- Re(C[i,j])
                }
        }
        return(C)
}



#' Quadratic Cost Matrix
#'
#' Calculating the quadratic cost matrix.
#'
#' @param X A numeric vector.
#' @param Y A numeric vector.
#' @return The quadratic cost matrix.
#' @export

quadCost <- function(X,Y){
        C <- matrix(rep(0, length(X)*length(Y)), nrow = length(X))

        for(i in 1:length(X)){
                for(j in 1:length(Y)){
                        C[i,j] <- abs(X[i]-Y[j])^2
                }
        }
        return(C)
}




#' Value in interval
#'
#' This checks if the value is in a given interval
#'
#' @param x A numeric value.
#' @param s1 The left boundary of a 1D interval.
#' @param s2 The right boundary of a 1D interval.
#' @return 1 if x is in [s1,s2] and 0 if it is not.
#'
inseg <- function(x,s1,s2){
        if(s1 <= x & x <= s2){
                return(1)
        }else{
                return(0)
        }

}

#' Example supply function
#'
#'
#' @param x A numeric value.
#' @return The function value r = p(x).
#' @export
#'
fp <- function(x){
        r <- 2*inseg(x,0.0,0.2)
        r <- r + 40*(x-0.9)*inseg(x,0.9,0.95)
        r <- r + 40*(1-x)*inseg(x,0.95,1)

        return(r)

}

#' Example demand function
#'
#'
#' @param y A numeric value.
#' @return The function value r = ´q(y).
#' @export
#'
fq <- function(y){
        r <- 10*(y-0.2)*inseg(y,0.2,0.4)
        r <- r + 1.3*(sqrt(as.complex(1-(y-0.7)^2/0.04)))*inseg(y,0.5,0.9)
        r <- Re(r)
        return(r)

}

#' Distance matrix
#'
#' Calculating the distance matrix for point clouds in 2D.
#'
#' @param x The first coordinates.
#' @param y The second coordinates.
#' @param method Determines which distance function to use for the computation. The default value is "euclidean" but every
#' method given for stats::dist can be used.
#' @return The distance matrix.
#' @export
createCostMatrix <- function(x,y,method = "euclidean"){
        costM <-matrix(rep(0, length(x)*length(x)), nrow = length(x))
        for (i in 1:length(x)){
                for (j in 1:length(x)){
                        costM[i,j] <- stats::dist(rbind(c(x[i],y[i]), c(x[j] ,y[j])), method = method)
                }
        }
        return(costM)
}


#' Cost Matrix
#' @description Calculates the cost matrix between points.
#'
#' @param x A numeric matrix. Each row corresponds to the coordinates of one point in the first point cloud.
#' @param y A numeric matrix. Each row corresponds to the coordinates of one point in the second point cloud
#' @param method Determines which distance function to use for the computation.  Currently implemented are:
#' 1. Euclidean
#' 2. Minkowski
#' 3. Maximum
#' @param exp exponent to apply to the distance
#' @param wfr computing the wasserstein fisher rao distance
#' @param p parameter for the minkowski metric. standard p = 2 give the minkowski metric.
#' @return The distance matrix between the points. The rows correspond to the points in x, the columns to the
#' points in y
#' @export
costMatrix <- function(x, y, method = "euclidean", exp = 1,  wfr = FALSE, p = 2){


        if(ncol(x) != ncol(y)){
                print("Unequal dimensions.")
        }

        if (method == "euclidean"){
                method <- "minkowski"
                p <- 2
        }

        cMatrix <- matrix(rep(0, nrow(x)*nrow(y)), ncol = nrow(x))

        if(method == "minkowski"){

                for (i in 1:nrow(x)){

                        cMatrix[i,] <- ((rowSums((t(t(y) - x[i,]))^p))^(1/p))^exp
                }

        }else if(method == "maximum"){

                for(i in 1:nrow(x)){

                        cMatrix[i,] <- max(t(t(y) - x[i,]))^exp

                }

        }else{
                print("Please specify a method.")
        }


        if(wfr){

                cMatrix <- -2*log(cos(min(cMatrix, pi/2)))
        }

        return(cMatrix)

}
#' costMatrix(x, y, method = "maximum")
#
#
# x <- matrix(c(1,1,1,0,1,2), ncol = 2)
# x
# y <- matrix(c(3,2,3,0,1,2), ncol = 2)
# y
#
# plot(x, xlim = c(1,3))
# points(y)




