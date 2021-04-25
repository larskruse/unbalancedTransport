#' ex2
#'
#' @export
#'
ex2 <- function(){
    # number of supply and demand points
    I <- 1000
    J <- 1000

    # discretisation of [0,1]
    X <- seq(0,1,length.out = I)
    Y <- seq(0,1,length.out = J)

    # reference Measures
    dx <- rep(1,I)/I
    dy <- rep(1,J)/J

    # supply measure
    p <- sapply(X, fp)
    p[1] <- 0
    p[I] <- 0

    # demand measure
    q <- sapply(Y, fq)
    q <- q*sum(p)/sum(q)

    # number of iterations
    iterMax <- 20000

    # vector of epsilon values
    epsvec <- seq(-1,-7, length.out = 20)
    epsvec <- 10^(epsvec)


    # supply/demand with discretization and reference measure
    supplyList <- list(p, X, dx)
    demandList <- list(q, Y, dy)


    #compute quadrature cost matrix
    C <-  quadCost(X,Y)
    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
    supply <- list("KL", p, dx, 0.5)
    demand <- list("KL", q, dy, 0.5)
    # compute and plot the transport plan
    map <- scalingAlgorithm(C, supply, demand, iterMax, epsvec)
    plot1DTransport(map$map, supplyList, demandList, X)
    plot(map)

    #compute quadrature cost matrix
    C <-  quadCost(X,Y)
    # use Total Variation with lambda = 0.05 for supply and demand
    supply <- list("TV", p, dx, 0.05)
    demand <- list("TV", q, dy, 0.05)
    # compute and plot the transport plan
    map <- scalingAlgorithm(C, supply, demand, iterMax, epsvec)
    plot1DTransport(map$map, supplyList, demandList, X)
    plot(map)


    #compute Wasserstein-Fisher-Rao cost matrix
    C <- wfrCost(X,Y)
    # use Kulback-Leibner divergence with lambda = 0.5 for supply and demand
    supply <- list("KL", p, dx, 0.5)
    demand <- list("KL", q, dy, 0.5)
    # compute and plot the transport plan
    map <- scalingAlgorithm(C, supply, demand, 20000, epsvec)
    plot1DTransport(map$map, supplyList, demandList, X)
    plotUOTP(map)

}




