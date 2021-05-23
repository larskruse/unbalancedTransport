#' The Scaling Algorithm Example
#'
#' The examples given in the [paper]
#'
#' @export
#'
ExampleScaling <- function(){


    print("This function uses the example distributions used in the scaling paper.")

    # number of supply and demand points
    I <- 1000
    J <- 1000

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
    iterMax <- 10000


    # vector of epsilon values
    epsvec <- seq(-1,-7, length.out = 20)
    epsvec <- 10^(epsvec)

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
    res <- scalingAlgorithmFromCost(C, supply, demand, iterMax, epsvec)


    cat("Transport cost: ", res$TransportCost, "\n\n\n")

    plot1DTransport(res$TransportPlan, supplyList, demandList)
    plotUOTP(res$TransportPlan)

    cat("Number of iterations: ", iterMax, "\n")

    cat("Epsilon values: ", epsvec, "\n")

    cat("The TV divergence with parameter 0.05 is used for both supply and demand vectors.\n")

    cat("The quadratic cost is used.\n")

    #compute quadrature cost matrix
    C <-  costMatrix(X, Y, exp = 2)
    # use Total Variation with lambda = 0.05 for supply and demand
    supply <- list("TV", p, 0.05)
    demand <- list("TV", q, 0.05)


    # compute and plot the transport plan
    res <- scalingAlgorithmFromCost(C, supply, demand, iterMax, epsvec)
    cat("Transport cost: ", res$TransportCost, "\n\n\n")
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    plotUOTP(res$TransportPlan)


    cat("Number of iterations: ", iterMax, "\n")

    cat("Epsilon values: ", epsvec, "\n")

    cat("The KL divergence with parameter 0.5 is used for both supply and demand vectors.\n")

    cat("The WFR distance matrix is used.\n")

    #compute Wasserstein-Fisher-Rao cost matrix
    C <- wfrCost(X,Y)
    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
    supply <- list("KL", p, 0.5)
    demand <- list("KL", q, 0.5)
    # compute and plot the transport plan
    res <- scalingAlgorithmFromCost(C, supply, demand, 20000, epsvec)
    cat("Transport cost: ", res$TransportCost, "\n\n\n")
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    plotUOTP(res$TransportPlan)

    cat("Number of iterations: ", iterMax, "\n")

    cat("Epsilon values: ", epsvec, "\n")

    cat("The RG divergence with parameters 0.7 and 1.2 is used for both supply and demand vectors.\n")


    C <-  costMatrix(X, Y, exp = 2)

    supply <- list("RG", p, 0.7, 1.2)
    demand <- list("RG", q, 0.7, 1.2)


    #compute and plot the transport plan
    res <- scalingAlgorithmFromCost(C, supply, demand, iterMax, epsvec)
    cat("Transport cost: ", res$TransportCost, "\n\n\n")
    plot1DTransport(res$TransportPlan, supplyList, demandList)
    plotUOTP(res$TransportPlan)


}

#' The Scaling Algorithm Example
#'
#' Another scaling algorithm example
#'
#' @export
#'
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


    supply <- list("RG", p, 0.1, 20)
    demand <- list("RG", q, 0.1, 20)



    # compute and plot the transport plan
    res <- scalingAlgorithmFromCost(C, supply, demand, iterMax, epsvec)
    transportP <- res$TransportPlan
    #

    print(res)

    plotUOTP(transportP)
    #
    plotTransportByCost(C, transportP, p, q)

}
