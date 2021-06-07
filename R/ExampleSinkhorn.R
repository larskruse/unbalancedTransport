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
    iterMax <- 600
    
    
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
    supply <- list(p, "KL", 1.)
    demand <- list(q, "KL", 1.)
    
    
    # compute and plot the transport plan
    
    start_time <- Sys.time()
    res <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, costMatrix = C)
    end_time <- Sys.time()
    
    print(end_time - start_time)
    
    #cat("Transport cost: ", res$TransportCost, "\n\n\n")
    
    print(res$TransportPlan[1:10,1:10])
    
    
    
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    gridPlotTransport(res$TransportPlan)
    
    
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
    epsvec <- 10^(-3)
    
    #compute quadrature cost matrix
    C <- createCostMatrix(xx,xy)
    
    
    
    supply <- list(p, "TV", 10)
    demand <- list(q, "TV", 10)
    
    
    # compute and plot the transport plan
    res <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, costMatrix = C)
    #transportP <- res$TransportPlan
    
    print(res)
    

}


#' example 2
#' @noRd
ExampleSink3 <- function(){    
    iterMax <- 1000
    
    
    cost <- c(0)
    costS <- c(0)
    # vector of epsilon values
    epsvec <- 10^(-3)
    
   
    nn <- 300
    pa1 <- seq(0,0.2, length.out=nn)
    pa1
    a1 <- rep(1,nn)
    a1[1] <- 0
    a1[nn] <- 0
    a1 <- a1/sum(a1)
    a1 <- 0.65*a1
    a1
    
    pa2 <- seq(0.9,1, length.out = nn)
    a2 <- 0.05-abs(pa2-0.95)
    a2[0] <- 0
    a2[nn] <- 0
    a2 = a2/sum(a2)
    a2 <- a2*0.35
    a2
    
    p = c(a1,a2)
    X = c(pa1, pa2)
    
    qa1 <- seq(0.2,0.4, length.out=nn)
    
    b1 <- seq(0,1,length.out = nn)
    b1[0] <- 0
    b1[nn] <- 0
    b1 <- b1/sum(b1)
    b1 <- b1*0.45
    
    
    qa2 <- seq(0.5,0.9, length.out = nn)
    
    b2 <- sqrt(abs(1-((qa2-0.7)/0.2)^2))
    b2[1] <- 0
    b2[nn] <- 0
    b2 <- b2/sum(b2)
    b2 <- 0.55*b2
    b2
    q <- c(b1, b2)
    Y <- c(qa1, qa2)
    

    # supply/demand with discretization and reference measure
    supplyList <- list(p, X)
    demandList <- list(q, Y)
    
    
    #compute quadrature cost matrix

    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
    # supply <- list(p, "RG", 0.7, 1.3, X)
    # demand <- list(q, "RG", 0.7, 1.3, Y)
    # 
    # supply <- list(p, "TV", 0.05, X)
    # demand <- list(q, "TV", 0.05, Y)

    # res <- scalingAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # 
    # cost <- c(ress$TransportCost)
    # costS <- c(res$TransportCost)
    # plot1DTransport(t(ress$TransportPlan), supplyList, demandList)
    # plot1DTransport(t(res$TransportPlan), supplyList, demandList)
    # gridPlotTransport(ress$TransportPlan)
    # gridPlotTransport(res$TransportPlan)
    # 
    # supply <- list(p, "KL", 0.037, X)
    # demand <- list(q, "KL", 0.037, Y)
    # res <- scalingAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # 
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # cost <- c(cost, ress$TransportCost)
    # costS <- c(costS, res$TransportCost)
    # 
    # plot1DTransport(t(ress$TransportPlan), supplyList, demandList)
    # plot1DTransport(t(res$TransportPlan), supplyList, demandList)
    # gridPlotTransport(ress$TransportPlan)
    # gridPlotTransport(res$TransportPlan)
    # 
    # 
    # 
    # supply <- list(p, "KL", 0.136, X)
    # demand <- list(q, "KL", 0.136, Y)
    # res <- scalingAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # 
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-10, exp = 2)
    # cost <- c(cost, ress$TransportCost)
    # costS <- c(costS, res$TransportCost)
    # 
    # plot1DTransport(t(ress$TransportPlan), supplyList, demandList)
    # plot1DTransport(t(res$TransportPlan), supplyList, demandList)
    # gridPlotTransport(ress$TransportPlan)
    # gridPlotTransport(res$TransportPlan)
    
    # supply <- list(p, "KL", 0.1, X)
    # demand <- list(q, "KL", 0.1, Y)
    #res <- scalingAlgorithm(supply, demand, iterMax, c(epsvec), tol = 1e-15, exp = 2)
    
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-15, exp = 2)
    
    
    # cost <- c(cost, ress[[2]])
    # costS <- c(costS, res$TransportCost)
    # 
    # plot1DTransport((ress$TransportPlan), supplyList, demandList)
    #plot1DTransport((res$TransportPlan), supplyList, demandList)
    # gridPlotTransport(ress$TransportPlan)
    #gridPlotTransport(res$TransportPlan)
    # 
    # 
    # print(cost)
    # print(costS)
    # 
    
    # supply <- list(p, "TV", 0.01, X)
    # demand <- list(q, "TV", 0.01, Y)
    # 
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-15, exp = 2)
    # plot1DTransport((ress$TransportPlan), supplyList, demandList)
    # 
    # gridPlotTransport(ress$TransportPlan)
    # 
    # supply <- list(p, "TV", 0.03, X)
    # demand <- list(q, "TV", 0.03, Y)
    # 
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-15, exp = 2)
    # plot1DTransport((ress$TransportPlan), supplyList, demandList)
    # 
    # gridPlotTransport(ress$TransportPlan)
    # 
    # 
    # supply <- list(p, "TV", 0.13, X)
    # demand <- list(q, "TV", 0.13, Y)
    # 
    # ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-15, exp = 2)
    # plot1DTransport((ress$TransportPlan), supplyList, demandList)
    # 
    # gridPlotTransport(ress$TransportPlan)
    
    supply <- list(p, "TV", 0.13, X)
    demand <- list(q, "TV", 0.13, Y)
    
    ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec, tol = 1e-8, exp = 2)
    plot1DTransport((ress$TransportPlan), supplyList, demandList)
    
    gridPlotTransport(ress$TransportPlan)
    
    print(ress$TransportCost)
    print(ress$RegCost)
    
    return(cost)
    
    
}
