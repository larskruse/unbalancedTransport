#' The Monge cost matrix example
#'
#' Using the same values as for the scaling algorithm example given in [paper]
#'
#' @export
ExampleMonge <- function(){
    # number of supply and demand points
    I <- 300
    J <- 300

    # discretisation of [0,1]
    X <- seq(0,1,length.out = I)
    Y <- seq(0,1,length.out = J)

    # reference Measures
    dx <- rep(1,I)
    dy <- rep(1,J)

    # supply measure
    p <- sapply(X, fp)
    p[1] <- 0
    p[I] <- 0

    # demand measure
    q <- sapply(Y, fq)
    q <- q*sum(p)/sum(q)


    # supply/demand with discretization and reference measure
    supplyList <- list(p, X, dx)
    demandList <- list(q, Y, dy)


    #compute quadrature cost matrix and check Monge property
    C <-  quadCost(X,Y)
    checkMongeProperty(C)

    # compute and plot the transport plan
    res <- mongeTransport(C, p, q, 1.1, 1.1)
    plot1DTransport(res$transportPlan, supplyList, demandList, X)
    plotUOTP(res$transportPlan, res$import, res$export)

    # compute and plot the transport plan
    res <- mongeTransport(C, p, q, 0.2, 0.2)
    plot1DTransport(res$transportPlan, supplyList, demandList, X)
    plotUOTP(res$transportPlan, res$import, res$export)

}
