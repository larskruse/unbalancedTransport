#' The Scaling Algorithm Example
#'
#' The examples given in the [paper]
#'
#' @noRd
#'
ExampleSink <- function(){
    

    # number of supply and demand points
    I <- 300
    J <- 300
    
    # Discretization of [0,1]
    X <- seq(0,1,length.out = I)
    Y <- seq(0,1,length.out = J)
    
    
    # supply measure
    p <- sapply(X, fp)
    p[1] <- 0
    p[I] <- 0
    
    # demand measure
    q <- sapply(Y, fq)
    q <- q*sum(p)/sum(q)
    
    # number of iterations
    iterMax <- 6000
    
    
    # vector of epsilon values
    epsvec <- 10^(-3)
    
    cat("Number of iterations: ", iterMax, "\n")
    
    cat("Epsilon values: ", epsvec, "\n")
    
    cat("The KL divergence with parameter 0.5 is used for both supply and demand. \n")
    
    cat("The quadratic cost is used. \n")
    
    
    # supply/demand with discretization and reference measure
    supplyList <- list(p, X)
    demandList <- list(q, Y)
    
    
    #compute quadrature cost matrix
    C <-  costMatrix(X,Y, exp = 2)
    
    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
    supply <- list("KL", p, 0.5)
    demand <- list("KL", q, 0.5)
    
    
    # compute and plot the transport plan
    
    start_time <- Sys.time()
    res <- sinkhornAlgorithmFromCost(C, supply, demand, iterMax, epsvec)
    end_time <- Sys.time()
    
    print(end_time - start_time)
    
    #cat("Transport cost: ", res$TransportCost, "\n\n\n")
    
    print(res$TransportPlan[1:10,1:10])
    
    
    
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    plotUOTP(res$TransportPlan)
    
    
}







#' example 2
#' @noRd
ExampleSink2 <- function(){    
    xx <- c(0,0.1,0.2,3,3.5,4)
    xy <- c(0,1,2,1.7,1,0.2)
    xy <- c(0,0,0,0,0,0)
    # supply measure
    p <- c(2,2,2,0,0,0)
    q <- c(0,0,0,1,2,3)
    
    # number of iterations
    iterMax <- 100
    
    # vector of epsilon values
    epsvec <- 10^(-7)
    
    #compute quadrature cost matrix
    C <- createCostMatrix(xx,xy)
    
    
    
    supply <- list("TV", p, 10)
    demand <- list("TV", q, 10)
    
    
    # compute and plot the transport plan
    res <- sinkhornAlgorithmFromCost(C, supply, demand, iterMax, epsvec)
    #transportP <- res$TransportPlan
    
    print(res)
    

}
