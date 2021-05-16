#' The Monge cost matrix example
#'
#' Using the same values as for the scaling algorithm example given in [paper]
#'
#' @export
ExampleMonge <- function(){

    print("This function uses the example distributions used in the scaling paper and the quadratic cost matrix.")


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
    C <-  costMatrix(X,Y, exp = 2)
    checkMongeProperty(C)

    print("The cost for mass creation and destruction is set to 1.1")

    # compute and plot the transport plan
    res <- mongeAlgorithm(C, p, q, 1.1, 1.1)
    plot1DTransport(res$transportPlan, supplyList, demandList)
    plotUOTP(res$transportPlan, res$import, res$export)


    print("The cost for mass creation and destruction is set to 0.2")

    # compute and plot the transport plan
    res <- mongeAlgorithm(C, p, q, 0.1, 0.1)
    plot1DTransport(res$transportPlan, supplyList, demandList)
    plotUOTP(res$transportPlan, res$import, res$export)



}
