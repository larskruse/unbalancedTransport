#' The tree metric transport example
#'
#' The tree metric transport problem from paper
#'
#' @export
ExampleTree1 <- function(){


    tree <- list(1, c(1, 2, 1), c(2, 3, 1), c(3, 4, 1), c(3, 5, 1),
                 c(2, 6, 1), c(1, 7, 1), c(7, 8, 1), c(7, 9, 1),
                 c(9, 10, 1), c(9, 11, 1), c(11, 12, 1), c(11, 13, 1))


    constructionCost <- rep(2,13)
    destructionCost <- rep(1,13)


    cat("Mass creation cost: ", constructionCost, "\n")
    cat("Mass destruction cost: ", destructionCost, "\n")

    supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
    demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)

    cat("Supply: ", supply, "\n")
    cat("Demand: ", demand, "\n")

    print("All edge weights are 1.")


    transport <- treeAlgorithm(tree, supply, demand , constructionCost,
                               destructionCost, output = "list")


    # Plotting the results
    plotTree(tree, tList = transport$transportList, supply = supply, demand = demand)

}

#' A second tree metric transport example
#'
#'
#'
#' @export
ExampleTree2 <- function(){


    tree <- list(1, c(1,2,1), c(2,3,1), c(2,4,1), c(1,5,1), c(5,6,1), c(5,7,1), c(1,8,10),
                 c(8,9,1), c(8,10,1))


    constructionCost <- rep(2,13)
    destructionCost <- rep(2,13)

    cat("Mass creation cost: ", constructionCost, "\n")
    cat("Mass destruction cost: ", destructionCost, "\n")

    supply <- c(0,0,2,0,0,4,0,0,2,0)
    demand <- c(0,0,0,0,7,0,0,0,0,1)

    cat("Supply: ", supply, "\n")
    cat("Demand: ", demand, "\n")

    print("The edge weight from 1 to 8 is 10. Every other edge weight is 1.")


    transport <- treeAlgorithm(tree, supply, demand , constructionCost,
                               destructionCost, output = "list")

    # Plotting the results
    plotTree(tree, tList = transport$transportList, supply = supply, demand = demand)

}

