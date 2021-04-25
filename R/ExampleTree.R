#' The tree metric transport example
#'
#' The tree metric transport problem from paper
#'
#' @export
ExampleTree <- function(){


    tree <- list(1, c(1, 2, 1), c(2, 3, 1), c(3, 4, 1), c(3, 5, 1),
                 c(2, 6, 1), c(1, 7, 1), c(7, 8, 1), c(7, 9, 1),
                 c(9, 10, 1), c(9, 11, 1), c(11, 12, 1), c(11, 13, 1))

    plotTree(tree)

    constructionCost <- rep(2,13)
    destructionCost <- rep(1,13)

    supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
    demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)


    transport <- treeAlgorithm(tree, supply, demand , constructionCost,
                               destructionCost, output = "list")

    # Plotting the results
    plotTree(tree, tList = transport$transportList, supply = supply, demand = demand)

}
