#' @importFrom transport transport




#' @title getCostMatrix
#'
#' @param tree a tree in list format
#'
#' @return the cost matrix for the given tree
#'
getCostMatrix <- function(tree){
    treeDF <- as.data.frame(do.call(rbind, tree[-1]))
    colnames(treeDF) <- c("parent", "child", "weight")
    treeDF <- treeDF[order(treeDF$parent),]
    rootNode <- tree[1]

    uniqueNodes <- sort(unique(c(treeDF$parent, treeDF$child)))
    costM <- matrix(rep(0, length(uniqueNodes)*length(uniqueNodes)), nrow = length(uniqueNodes))


    for(i in uniqueNodes){
        for(j in uniqueNodes[uniqueNodes < i]){


            pathFrom <- unlist(findPath(rootNode,i,treeDF))
            pathTo <- unlist(findPath(rootNode,j, treeDF))


            if(length(pathFrom) > 1 & length(pathTo) > 1){
                while(length(pathFrom) > 1 & length(pathTo) > 1 & pathFrom[2] == pathTo[2]){
                    pathTo <- pathTo[-1]
                    pathFrom <- pathFrom[-1]
                }
            }


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



            costM[which(uniqueNodes == i), which(uniqueNodes == j)] <- cost
            costM[which(uniqueNodes == j), which(uniqueNodes == i)] <- cost

        }

    }


    return(costM)



}





#' Title tree Algorithm
#'
#' @param tree tree in list format
#' @param supply supply vector
#' @param demand demand vector
#' @param creationCost creation cost vector
#' @param destructionCost destruction cost vector
#' @param output determines the output format
#' @param plot whether or not the results are plotted
#'
#' @return the transport cost, import and export vector, transport list or plan if specified in 'output'
#' @export
#'
treeAlgorithm <- function(tree, supply, demand, creationCost, destructionCost, output = "cost", plot = "FALSE"){

    if(output != "list" & output != "transportPlan"  & output != "both"  & output != "cost"){

        print("Output set to 'cost'")
        output <- "cost"

    }

    treegkrOut <- treegkr_Rcpp(tree[-1], supply, demand, creationCost, destructionCost)

    importVec = treegkrOut$import

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

    costMatrix <- getCostMatrix(tree)
    res <- transport::transport((supply-export), (demand-import), costMatrix,method = "revsimplex")


    transportList <- list()

    if(length(res$from) > 0){

        for(i in 1:length(res$from)){
            if(res$mass[i] != 0){
                transportList[[length(transportList)+1]] <- c(res[[1]][i], res[[2]][i], res[[3]][i])

            }

        }
    }else{

    }

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

    }else if (output == "both"){

        if(length(res$from) > 0){

            tPlan <- matrix(0,length(supply),length(demand))
            tPlan[cbind(res$from,res$to)] <- res$mass

        }else{
            tPlan <- matrix(0,length(supply),length(demand))

        }

        result <- list(treegkrOut$cost, tPlan, transportList, import, export)
        names(result) <- c("cost", "transportPlan", "transportList", "import", "export")
        return(result)

    }else{
        result <- list(treegkrOut$cost, import, export)
        names(result) <- c("cost", "import", "export")
        return(result)
    }
}

