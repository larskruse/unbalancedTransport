#' The Scaling Algorithm Example
#'
#' The examples given in the [paper]
#'
#' @noRd
#'
ExampleScaling <- function(){


    print("This function uses the example distributions used in the scaling paper.")

    # number of supply and demand points
    I <- 2000
    J <- 2000

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
    iterMax <- 4000
    
    # vector of epsilon values
    epsvec <- seq(-1,-7, length.out = 20)
    epsvec <- 10^(epsvec)
    #epsvec <- c(10^(-3))

    cat("Number of iterations: ", iterMax, "\n")

    cat("Epsilon values: ", epsvec, "\n")

    cat("The KL divergence with parameter 0.5 is used for both supply and demand. \n")

    cat("The quadratic cost is used. \n")


    # supply/demand with discretization and reference measure
    supplyList <- list(p, X)
    demandList <- list(q, Y)


    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
   

    supply <- list(p, "RG", 0.7, 1.2, X)
    demand <- list(q, "RG", 0.7, 1.2, Y)

    # compute and plot the transport plan
    res <- scalingAlgorithm(supply, demand, iterMax, epsvec, exp = 2, tol = 10e-10)
    
    #ress <- sinkhornAlgorithm(supply, demand, iterMax, 10e-3, exp = 2, tol = 10e-10)

    
    #plot1DTransport(ress$TransportPlan, supplyList, demandList)
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    #gridPlotTransport(ress$TransportPlan)
    gridPlotTransport(res$TransportPlan)
    
    return(res)
}

#aa <- ExampleScaling()


#print(ress$TransportCost)
    
    # num <- c(1.5,1.25,1,0.75,0.5,0.25,0.05,0.001)
    # 
    # for(i in 1:8){
    #     supply <- list(p, "TV", num[i], X)
    #     demand <- list(q, "TV", num[i], Y)
    #     ress <- sinkhornAlgorithm(supply, demand, iterMax, 10e-2, exp = 2, tol = 10e-5)
    # 
    # 
    #     plot1DTransport(ress$TransportPlan, supplyList, demandList)
    # 
    # 
    # }
    # 
    # for(i in 1:8){
    #     supply <- list(p, "TV", num[i], X)
    #     demand <- list(q, "TV", num[i], Y)
    #     ress <- sinkhornAlgorithm(supply, demand, iterMax, 10e-2, exp = 2, tol = 10e-5)
    # 
    # 
    #     gridPlotTransport(ress$TransportPlan)
    # 
    # }
    # 
    # supply <- list(p, "TV", 1.5, X)
    # demand <- list(q, "TV", 1.5, Y)
    # 
    # for(i in 1:5){
    #     ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec[2*i], exp = 2, tol = 10e-5)
    #     plot1DTransport(ress$TransportPlan, supplyList, demandList)
    #     
    #     
    # }
    # 
    # for(i in 1:5){
    #     ress <- sinkhornAlgorithm(supply, demand, iterMax, epsvec[2*i], exp = 2, tol = 10e-5)
    #     gridPlotTransport(ress$TransportPlan)
    #     
    # }
    # 
    # 
    # for(i in 1:5){
    #     supply <- list(p, "TV", num[i], X)
    #     demand <- list(q, "TV", num[i], Y)
    #     ress <- scalingAlgorithm(supply, demand, iterMax, c(10e-2), exp = 2, tol = 10e-15)
    #     
    #     
    #     plot1DTransport(ress$TransportPlan, supplyList, demandList)
    #     
    #     
    # }
    # 
    # for(i in 1:5){
    #     supply <- list(p, "TV", num[i], X)
    #     demand <- list(q, "TV", num[i], Y)
    #     ress <- scalingAlgorithm(supply, demand, iterMax, c(10e-2), exp = 2, tol = 10e-15)
    #     
    #     
    #     gridPlotTransport(ress$TransportPlan)
    #     
    # }
    # 
    # 
    # supply <- list(p, "TV", 0.2, X)
    # demand <- list(q, "TV", 0.2, Y)
    # 
    # for(i in 1:4){
    #     ress <- scalingAlgorithm(supply, demand, iterMax, c(epsvec[2*i]), exp = 2, tol = 10e-15)
    #     plot1DTransport(ress$TransportPlan, supplyList, demandList)
    #     
    #     
    # }
    # 
    # for(i in 1:4){
    #     ress <- scalingAlgorithm(supply, demand, iterMax, c(epsvec[2*i]), exp = 2, tol = 10e-15)
    #     gridPlotTransport(ress$TransportPlan)
    #     
    # }
    # 
    # 
    # 
    # return(0)
    
# }


# 
#     cat("Number of iterations: ", iterMax, "\n")
# 
#     cat("Epsilon values: ", epsvec, "\n")
# 
#     cat("The TV divergence with parameter 0.05 is used for both supply and demand vectors.\n")
# 
#     cat("The quadratic cost is used.\n")
# 
#     #compute quadrature cost matrix
#     C <-  costMatrix(X, Y, exp = 2)
#     # use Total Variation with lambda = 0.05 for supply and demand
#     
#     supply <- list(p, "TV", 0.05, X)
#     demand <- list(q, "TV", 0.05, Y)
# 
#     
# 
#     # compute and plot the transport plan
#     res <- scalingAlgorithm(supply, demand, iterMax, epsvec, exp = 2)
#     cat("Transport cost: ", res$TransportCost, "\n\n\n")
#     plot1DTransport(res$TransportPlan, supplyList, demandList)
#     gridPlotTransport(res$TransportPlan)
# 
# 
#     cat("Number of iterations: ", iterMax, "\n")
# 
#     cat("Epsilon values: ", epsvec, "\n")
# 
#     cat("The KL divergence with parameter 0.5 is used for both supply and demand vectors.\n")
# 
#     cat("The WFR distance matrix is used.\n")
# 
#     #compute Wasserstein-Fisher-Rao cost matrix
#     C <- wfrCost(X,Y)
#     # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
#     supply <- list(p, "KL", 0.5)
#     demand <- list(q, "KL", 0.5)
#     # compute and plot the transport plan
#     res <- scalingAlgorithm(supply, demand, 20000, epsvec, costMatrix = C)
#     cat("Transport cost: ", res$TransportCost, "\n\n\n")
#     plot1DTransport(res$TransportPlan, supplyList, demandList)
#     gridPlotTransport(res$TransportPlan)
# 
#     cat("Number of iterations: ", iterMax, "\n")
# 
#     cat("Epsilon values: ", epsvec, "\n")
# 
#     cat("The RG divergence with parameters 0.7 and 1.2 is used for both supply and demand vectors.\n")
# 
# 
#     C <-  costMatrix(X, Y, exp = 2)
# 
#     supply <- list(p, "RG", 0.7, 1.2, X)
#     demand <- list(q, "RG", 0.7, 1.2, Y)
# 
# 
#     #compute and plot the transport plan
#     res <- scalingAlgorithm(supply, demand, iterMax, epsvec, exp = 2)
#     cat("Transport cost: ", res$TransportCost, "\n\n\n")
#     plot1DTransport(res$TransportPlan, supplyList, demandList)
#     gridPlotTransport(res$TransportPlan)
# 
# 
# }

#' The Scaling Algorithm Example
#'
#' Another scaling algorithm example
#'
#' @noRd
ExampleScaling2 <- function(){


    xx <- c(0,0,0,3,3.5,4)
    xy <- c(0,1,2,1.7,1,0.2)

    # supply measure
    p <- c(2,2,2,0,0,0)
    q <- c(0,0,0,1,2,3)

    # number of iterations
    iterMax <- 10000

    # vector of epsilon values
    epsvec <- seq(-1,-7, length.out = 20)
    epsvec <- 10^(epsvec)
    #epsvec <- c(2)

    #compute quadrature cost matrix
    C <- createCostMatrix(xx,xy)


    cat("Number of iterations: ", iterMax, "\n")

    cat("Epsilon values: ", epsvec, "\n")

    cat("The KL divergence with parameter 0.5 is used for both supply and demand vectors.")

    cat("The euclidean distance matrix is used. \n")


    supply <- list(p, "TV", 0.1, 20)
    demand <- list(q, "TV", 0.1, 20)


    # compute and plot the transport plan
    res <- scalingAlgorithm(supply, demand, iterMax, epsvec, costMatrix = C)
    transportP <- res$TransportPlan
    #

    print(res)

    gridPlotTransport(transportP)
    #
    plotTransportByCost(C, transportP, p, q)

}
