#' @title Tree Cost Matrix
#'
#' Calculating the cost matrix for a given tree.
#'
#' @param tree A tree in list format. Every element in the list is a vector
#' containing the parent node, the child node and the edge weight.
#'
#' @return The cost matrix for the given tree.
#' @noRd
getCostMatrix <- function(tree){
    treeDF <- as.data.frame(do.call(rbind, tree[-1]))
    colnames(treeDF) <- c("parent", "child", "weight")
    treeDF <- treeDF[order(treeDF$parent),]
    rootNode <- tree[1]

    uniqueNodes <- sort(unique(c(treeDF$parent, treeDF$child)))
    costM <- matrix(rep(0, length(uniqueNodes)*length(uniqueNodes)), nrow = length(uniqueNodes))

    # Calulate the path between every pair of nodes
    for(i in uniqueNodes){
        for(j in uniqueNodes[uniqueNodes < i]){

            # Finding the paths between the root node and the two current nodes
            pathFrom <- unlist(findPath(rootNode,i,treeDF))
            pathTo <- unlist(findPath(rootNode,j, treeDF))

            # Removing every but the last node that is in both paths from the
            # paths
            if(length(pathFrom) > 1 & length(pathTo) > 1){
                while(length(pathFrom) > 1 & length(pathTo) > 1 & pathFrom[2] == pathTo[2]){
                    pathTo <- pathTo[-1]
                    pathFrom <- pathFrom[-1]
                }
            }

            # Calculating the cost of both paths by adding up the edge weights
            cost <- 0
            if(length(pathTo) > 1){
                for(k in 1:(length(pathTo)-1)){
                    cost <- cost + treeDF[(treeDF$parent == pathTo[k] & treeDF$child == pathTo[k+1]),]$weight
                }
            }

            if(length(pathFrom) > 1){
                for(k in 1:(length(pathFrom)-1)){
                    cost <- cost + treeDF[(treeDF$parent == pathFrom[k] & treeDF$child == pathFrom[k+1]),]$weight

                }
            }


            # Since the transport matrix is symmetric, the cost is assigned to
            # both matrix entries.
            costM[which(uniqueNodes == i), which(uniqueNodes == j)] <- cost
            costM[which(uniqueNodes == j), which(uniqueNodes == i)] <- cost

        }

    }


    return(costM)

}





#' Transport on Trees
#' 
#' Solving unbalanced optimal transport problems that use a tree metric as cost function.
#' 
#' If the cost matrix is derived from a tree metric, the unbalanced optimal transport problem  \eqn{\min_r <C,r> + \sum_i p_i(\alpha_i - \sum_j r_{ij}) + \sum_j q_j(\beta-\sum_i r{ij})}{
#' min_r <C,r> + sum_i(p_i (a-sum_j r_ij)) + sum_j(q_j (b-sum_i r_ij))} with supply and demand measure \eqn{\alpha}{a} and \eqn{\beta}{b}
#' , transport plan \eqn{r}, construction and destruction costs \eqn{p} and \eqn{q} and cost matrix \eqn{C} can be solved efficiently in
#' O(n log²n) time. The C++ implementation of the algorithm used in this package can be found at https://github.com/joisino/treegkr. It was updated in order to make it accessible 
#' from R and to be able to compute the transport map. 
#'
#' @param tree A tree in list format. The first element is the index of the
#' root node. Every other element holds a vector describing the tree's edges: The parent node, the child node and the weight.
#' @param supply A non negative numeric supply vector. 
#' @param demand A non negative numeric demand vector
#' @param creationCost A non negative vector giving the cost to create mass at each node.
#' @param destructionCost A non negative vector giving the cost to destruct mass at each node.
#' @param output Determines the format of the output:
#' \itemize{
#' \item "transportPlan" give a standard transport plan matrix where the value at \eqn{i,j} gives the amount of mass transproted from 
#' note \eqn{i} to \eqn{j}.
#' \item "list" gives a list of vectors. Each vector holds three entries: The first gives the supply node, the second the target node of the transport
#' and the third the amount of mass transported between those nodes.
#' \item "cost" gives only the transport cost
#' }
#'
#' @return A list containing the transport cost, the mass import and export vectors, and the transport list or plan if specified in 'output'.
#' @export
#'
#' @examples
#'
#' tree <- list(1, c(1,2,1), c(2,3,1), c(2,4,1), c(1,5,1), c(5,6,1),
#'                 c(5,7,1), c(1,8,10), c(8,9,1), c(8,10,1))
#'
#' constructionCost <- rep(2,13)
#' destructionCost <- rep(2,13)
#'
#' supply <- c(0,0,2,0,0,4,0,0,2,0)
#' demand <- c(0,0,0,0,7,0,0,0,0,1)
#'
#' transport <- treeAlgorithm(tree, supply, demand,
#'                            constructionCost, destructionCost,
#'                            output = "list")
#'
#' plotTree(tree, tList = transport$transportList,
#'          supply = supply, demand = demand)
#'
#'
#'
#'
#'tree <- list(1, c(1,2,1), c(2,3,1), c(3,4,1), c(3,5,1), c(2,6,1),
#'                c(1,7,1), c(7,8,1), c(7,9,1), c(9,10,1), c(9,11,1),
#'                c(11,12,1), c(11,13,1))
#'
#'
#' constructionCost <- rep(2,13)
#' destructionCost <- rep(1,13)
#'
#'
#' supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
#' demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)
#'
#' transport <- treeAlgorithm(tree, supply, demand,
#'                            constructionCost, destructionCost,
#'                            output = "list")
#'
#' plotTree(tree, tList = transport$transportList,
#'          supply = supply, demand = demand)
#'
treeAlgorithm <- function(tree, supply, demand, creationCost, destructionCost, output = "transportPlan"){


    if(output != "list" & output != "transportPlan"  & output != "cost"){

        print("Output set to 'cost'")
        output <- "cost"

    }

    treegkrOut <- treegkr_Rcpp(tree[-1], supply, demand, creationCost, destructionCost)

    importVec = treegkrOut$import

    # importVec hat positive and negative entries. The positive entries indicate
    # imported mass, the negative entries exported mass.
    # Splitting the vector in import and export
    import <- rep(0,length(supply))
    export <- rep(0,length(supply))
    for(i in 1:length(importVec)){
        if(importVec[i] != 0){
            if(importVec[i] > 0){
                import[i] <- import[i] + importVec[i]
            }else{
                export[i] <- export[i] - importVec[i]
            }

        }

    }

    # Compute the cost matrix from the tree.
    costMatrix <- getCostMatrix(tree)

    # Computing the transport plan using the revised simplex algorithm.
    res <- transport::transport((supply-export), (demand-import), costMatrix,method = "revsimplex")


    # Computing the transport plan in list form. Each entry consists of the node the mass
    # comes from, the node the mass is transported to and the amount of mass.
    transportList <- list()
    if(length(res$from) > 0){

        for(i in 1:length(res$from)){
            if(res$mass[i] != 0){
                transportList[[length(transportList)+1]] <- c(res[[1]][i], res[[2]][i], res[[3]][i])
            }

        }
    }else{

    }

    # return a list according to the "output" variable
    if(output == "list"){

        result <- list(treegkrOut$cost, transportList, import, export)
        names(result) <- c("cost", "transportList", "import", "export")
        return(result)
    }else if (output == "transportPlan"){

        if(length(res$from) > 0){

            tPlan <- matrix(0,length(supply),length(demand))
            tPlan[cbind(res$from,res$to)] <- res$mass

        }else{
            tPlan <- matrix(0,length(supply),length(demand))

        }
        result <- list(treegkrOut$cost, tPlan, import, export)
        names(result) <- c("cost", "transportPlan", "import", "export")

    }else{
        result <- list(treegkrOut$cost, import, export)
        names(result) <- c("cost", "import", "export")
        return(result)
    }
}

