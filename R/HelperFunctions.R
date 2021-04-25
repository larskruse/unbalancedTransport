
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
#'
quadCost <- function(X,Y){
        C <- matrix(rep(0, length(X)*length(Y)), nrow = length(X))

        for(i in 1:length(X)){
                for(j in 1:length(Y)){
                        C[i,j] <- abs(X[i]-Y[j])^2
                }
        }
        return(C)
}




#' Interval check
#'
#' This checks if a value is in a given interval
#'
#' @param x A numeric value.
#' @param s1 The left boundary of a 1D interval.
#' @param s2 The right boundary of a 1D interval.
#' @return 1 if x is in [s1,s2] and 0 if it is not.
#' @export
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


