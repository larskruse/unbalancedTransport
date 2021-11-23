#' Plotting transport between single points.
#'
#' Visualizing mass transport between point clouds. 
#' 
#' 
#'
#' @param transportPlan A non negative numeric matrix that indicates where the mass is
#' transported. The value at point \eqn{\[i,j\]} is the amount of mass transported from
#' supply point \eqn{i} to demand point \eqn{j}.
#' @param supplyList A list containing the information about the supply measure. The first element hast to be the mass
#' distribution followed by
#' a vector specifying the cost for destruction of mass at each supply point. Lastly, the positions of the supply points. 
#' Has to be provided. This can either be a vector or a matrix where each row gives the coordinates for one point.
#' @param demandList A list similar to the supplyList but holding the information about the demand distribution.
#' 
#' \if{html}{\figure{plotTransportPoints.png}}
#' \if{latex}{\figure{plotTransportPoints.png}{options: width=0.5in}}
#' 
#' @examples 
#' 
#' supplyPoints <- matrix(c(0,0,0,1), ncol = 2)
#' demandPoints <- matrix(c(3.5,4,1,0.2), ncol = 2)
#' 
#' creationCost <- rep(1.8,2)
#' destructionCost <- rep(2,2)
#' 
#' 
#' supplyMass <- c(1,2)
#' demandMass <- c(2.5,1)
#' 
#' supplyList <- list(supplyMass, destructionCost,supplyPoints)
#' demandList <- list(demandMass, creationCost, demandPoints)
#' 
#' res <- umtp(supplyList, demandList)
#' transportPlan <- res$transportPlan
#' 
#' plotTransportPoints(transportPlan, supplyList, demandList)
#'
#'
#' @export

plotTransportPoints <- function(transportPlan, supplyList, demandList){
    
    
    if(length(supplyList) == 2 ){
        
        supply <- supplyList[[1]]
        costDestruct <- rep(0, length(supplyList[[1]]))
        supplyCoordinates <- supplyList[[2]]
        
        
    }else{
        
        supply <- supplyList[[1]]
        costDestruct <- supplyList[[2]]
        supplyCoordinates <- supplyList[[3]]
        
        
    }
    
    if(length(demandList) == 2){
        
        demand <- demandList[[1]]
        costCreate <- rep(0,length(demandList[[1]]))
        demandCoordinates <- demandList[[2]]
        
    }else{
        
        demand <- demandList[[1]]
        costCreate <- demandList[[2]]
        demandCoordinates <- demandList[[3]]
        
    }
    
    
    if(is.null(costCreate)){
        costCreate <- rep(0, ncol(transportPlan))
    }
    if(is.null(costDestruct)){
        costDestruct <- rep(0, nrow(transportPlan))
        
    }
    
    if(is.null(dim(supplyCoordinates))){
        
        
        
        
        xS <- supplyCoordinates
        yS <- rep(0, length(supplyCoordinates))
        
        xD <- demandCoordinates
        yD <- rep(0,length(demandCoordinates))
        
        print(xS)
        print(yS)
        print(xD)
        print(yD)
        
    }else{
        xS <- supplyCoordinates[,1]
        yS <- supplyCoordinates[,2]
        
        xD <- demandCoordinates[,1]
        yD <- demandCoordinates[,2]
    }
    
    if(length(xS)+length(xD) > 10){
        plot(1, type = "n", xlab = "", ylab = "",
             xlim = c(min(c(xS,xD))-0.1*abs(max(c(xS,xD))-min(c(xS,xD))),max(c(xS,xD))+0.1*abs(max(c(xS,xD))-min(c(xS,xD)))),
             ylim =c(min(c(yS,yD))-0.1*abs(max(c(yS,yD))-min(c(yS,yD))),max(c(yS,yD))+0.1*abs(max(c(yS,yD))-min(c(yS,yD)))),
             asp = 1)
        if(!is.null(transportPlan)){
            for(i in 1:length(xS)){
                
                for(j in 1:length(xD)){
                    
                    if(transportPlan[i,j] > 0){
                        curvedarrow(c(xS[i],yS[i]),c(xD[j],yD[j]), lwd = 0.3,
                                    arr.pos = 0.5, arr.adj = 0.5, arr.type = "none",
                                    curve = 0, lcol = "red", arr.col = "red")
                        
                    }
                    
                }
                
            }
            
        }
        
        points(c(xS,xD), c(yS,yD) , pch = 1)
        points(xS[which(supply > 0 )],yS[which(supply > 0)],
               pch = 19,  col = "chartreuse3")
        points(xD[which(demand > 0 )],yD[which(demand > 0)], pch = 19,
                 col = "dodgerblue3")
        
        
    }else{
        
        
        
        # create an empty plot
        plot(1, type = "n", xlab = "", ylab = "",
             xlim = c(min(c(xS,xD))-0.3*abs(max(c(xS,xD))-min(c(xS,xD))),max(c(xS,xD))+0.3*abs(max(c(xS,xD))-min(c(xS,xD)))),
             ylim =c(min(c(yS,yD))-0.3*abs(max(c(yS,yD))-min(c(yS,yD))),max(c(yS,yD))+0.3*abs(max(c(yS,yD))-min(c(yS,yD)))),
             asp = 1)
        
        
        if(!is.null(transportPlan)){
            for(i in 1:length(xS)){
                
                for(j in 1:length(xD)){
                    
                    if(transportPlan[i,j] > 0){
                        curvedarrow(c(xS[i],yS[i]),c(xD[j],yD[j]), lwd = 0.5+transportPlan[[i,j]],
                                    arr.pos = 0.5, arr.adj = 0.5, arr.type = "triangle",
                                    curve = 0.0, lcol = "red", arr.col = "red")
                        
                    }
                    
                }
                
            }
            
        }
        

        points(c(xS,xD), c(yS,yD) , pch = 1)
        points(xS[which(supply > 0 )],yS[which(supply > 0)],
               pch = 19, cex = supply[which(supply > 0)],  col = "chartreuse3")
        points(xD[which(demand > 0 )],yD[which(demand > 0)], pch = 19,
               cex = abs(demand[which(demand > 0)]),  col = "dodgerblue3")
        
        
        # if a transport plan is given, add arrows to indicate mass transport
        if(!is.null(costDestruct) & !is.null(costCreate)){
            
            # plot circles as indicators for the creation/destruction cost
            theta = seq(0, 2 * pi, length = 200)
            
            for(i in 1:length(xS)){
                
                lines(x = costDestruct[i]*cos(theta) + xS[i], y = costDestruct[i]*sin(theta) + yS[i], col = 'chartreuse3')
                
            }
            
            for(i in 1:length(xD)){
                lines(x = costCreate[i]*cos(theta) + xD[i], y = costCreate[i]*sin(theta) + yD[i], col = 'dodgerblue3')
            }
            
        }
        
    }
    
    
}


#' A grid plot of the transport plan
#'
#' Creates a grid plot of a transport plan. 
#' 
#' This function creates a grid plot of the transport plan. Import and export
#' vectors can be given as additional arguments. These will be plotted along the sides of the grid 
#' plot in order to indicate where mass is added or removed. Therefore, the mass shown in
#' each row or column is equal to the supply or demand vector. 
#'
#'
#' @param transportPlan A non negative numeric matrix that indicates the mass transport.
#' The value at \eqn{\[i,j\]} is the amount of mass transported from supply point \eqn{i} to demand point \eqn{j}.
#' @param import (optional) A non negative numeric vector that give the amount of mass created at each
#' demand point. It length has to be equal to the number of columns in the transport matrix.  
#' @param export (optional) A non negative numeric vector that give the amount of mass destroyed  at each
#' supply point. It length has to be equal to the number of rows in the transport matrix.  
#'
#'
#' \if{html}{\figure{Grid.png}}
#' \if{latex}{\figure{Grid.png}{options: width=0.5in}}
#'
#'
#' @examples
#' I <- 200
#' J <- 200
#' X <- seq(0,1,length.out = I)
#' Y <- seq(0,1,length.out = J)
#' p <- supplyExample[seq(0,1000,by = 1000/I)]
#' q <- demandExample[seq(0,1000, by = 1000/J)]
#' 
#' supply <- list(p,X)
#' demand <- list(q,Y)
#' 
#' maxIter <- 200
#' eps <- 1e-3
#' 
#' 
#' suppyDiv <- list("KL", 0.04)
#' demandDiv <- list("KL", 0.04)
#' res <- regularizedTransport(supply, demand, suppyDiv, demandDiv,
#'                             maxIteration = maxIter, epsVector = eps, p = 2)
#' plotGridTransport(res$TransportPlan)
#'
#' @export
plotGridTransport <- function(transportPlan, import = NULL, export =  NULL){

    # If no import or export vector is given, only the transport plan is plotted
    if(is.null(import) | is.null(export)){


        transportPlan <- t(transportPlan[(nrow(transportPlan)):1,])

        image(transportPlan, asp = 1, axes = FALSE, ylab = "Supply", xlab = "Demand",
              col=hcl.colors(20, palette = "YlOrRd", alpha = NULL, rev = TRUE, fixup = TRUE))


    }else{

        # Distance between import/export vectors and transport plan
        emptyColRow <- max(1, ceiling(max(nrow(transportPlan), ncol(transportPlan))/20))

        leftCol <- export
        leftCol[(length(leftCol)+1):(length(leftCol)+1+emptyColRow)] <- NaN

        # Binding the import vector to the transport plan
        addRows <- matrix(rep(NaN, ncol(transportPlan)*emptyColRow), ncol = ncol(transportPlan))
        printObj <- rbind(transportPlan, addRows ,import)

        # Binding the export vector to the transport plan
        addCols <- matrix(rep(NaN, nrow(printObj)*emptyColRow), nrow = nrow(printObj))
        printObj <- cbind(leftCol, addCols, printObj)
        printObj <- t(printObj[nrow(printObj):1,])

        # Plotting the transport plan with import and export vectors
        image(printObj, asp = 1, axes = FALSE, ylab = "Supply", xlab = "Demand",
              col=hcl.colors(20, palette = "viridis", alpha = NULL, rev = TRUE, fixup = TRUE))

        # Adding the labels and axes.
        att2 <- ((emptyColRow+1):(nrow(printObj)-1))/(nrow(printObj)-1)
        att2 <- att2[seq(1,length(att2), length.out = 10)]
        lab2 <- nrow(transportPlan):1
        lab2 <- lab2[seq(1, length(lab2), length.out = 10)]
        att1 <- ((emptyColRow+1):(ncol(printObj)-1))/(ncol(printObj)-1)
        att1 <- att1[seq(1,length(att1), length.out = 10)]
        lab1 <- 1:ncol(transportPlan)
        lab1 <- lab1[seq(1, length(lab1), length.out = 10)]
        Axis(side = 2, at = att2, labels = lab2)
        Axis(side = 1, at = att1, labels = lab1)

    }

}



#' Plotting dense 1D transport
#' 
#' A function to plot optimal transport between discretizations of continuous supply and demand measures in 1D. 
#' 
#' 
#' @param transportPlan A non negative numeric matrix giving the mass transport. The value at \eqn{(i,j)}
#' gives the mass transported from supply point \eqn{i} to demand point \eqn{j}.
#' @param supplyList A list containing the non negative supply measure and the underlying discretization as vectors.
#' @param demandList A list containing the non negative demand measure and the underlying discretization as vectors.
#'
#'
#' \if{html}{\figure{1D.png}}
#' \if{latex}{\figure{1D.png}{options: width=0.5in}}
#'
#'
#' @examples 
#'  
#' I <- 1000
#' J <- 1000
#' X <- seq(0,1,length.out = I)
#' Y <- seq(0,1,length.out = J)
#' p <- supplyExample
#' q <- demandExample
#' 
#' supply <- list(p,X)
#' demand <- list(q,Y)
#' 
#' maxIter <- 200
#' eps <- 1e-3
#' 
#' 
#' suppyDiv <- list("KL", 0.04)
#' demandDiv <- list("KL", 0.04)
#' res <- regularizedTransport(supply, demand, suppyDiv, demandDiv,
#'                             maxIteration = maxIter, epsVector = eps, p = 2)
#' plot1DTransport(res$TransportPlan, supply, demand)
#' 
#'
#' @export
plot1DTransport <- function(transportPlan, supplyList, demandList){

    # Number of color intervals
    numIntervals <- 50

    X <- supplyList[[2]]

    x1Measure <- rep(1, length(supplyList[[1]]))
    y1Measure <- rep(1, length(demandList[[1]]))

    colors <- rainbow(numIntervals, s = 1, v = 1, start = 0, end = max(1, numIntervals - 1)/numIntervals, alpha = 1, rev = FALSE)

    ymax <- max(supplyList[[1]], demandList[[1]], t(transportPlan) %*% x1Measure, transportPlan %*% y1Measure)

    # Plotting the supply and demand measures
    # plot(1, type = "n", xlab = "", ylab = "", xlim = c(min(X), max(X)), ylim= c(0, ymax))

    plot(supplyList[[2]], supplyList[[1]], type = "l", lty = 3 , col = "blue", ylim= c(0, ymax), xlab = "Position", ylab = "Mass")
    

    firstSupp <- supplyList[[2]][1]
    firstDem <- demandList[[2]][1]

    # Plot the intervals for different colors
    for( i in 1:numIntervals){

        # Calculating the interval of X
        colSuppInter <- rep(1, length(X))
        colSuppInter[(i-1)/numIntervals > X | i/numIntervals < X] <- 0

        # Restricting the transport map on the interval
        subK <- transportPlan * colSuppInter

        # Adding the color intervals
        #lines(demandList[[2]], t(subK) %*% x1Measure, type = "h", col = colors[i])
        #lines(supplyList[[2]], subK %*% y1Measure, type = "h", col = colors[i])

        # polygon(c(firstDem,demandList[[2]]),c(0,t(subK) %*% x1Measure), col = colors[i])
        # polygon(c(firstSupp,supplyList[[2]]), c(0,subK %*% y1Measure), col = colors[i])
        
        polygon(c(firstDem,demandList[[2]]),c(0,colSums(subK)), col = colors[i])
        polygon(c(firstSupp,supplyList[[2]]),c(0,rowSums(subK)) , col = colors[i])

        #lines(demandList[[2]], t(subK) %*% x1Measure, type = "l", col = "black")
        #lines(supplyList[[2]], subK %*% y1Measure, type = "l", col = "black")

    }

    
    
    lines(demandList[[2]], t(transportPlan) %*% x1Measure, type = "l", col = "green")
    lines(supplyList[[2]], transportPlan %*% y1Measure, type = "l", col = "blue")
    
    lines(supplyList[[2]], rep(0, length(supplyList[[2]])), type = "l", col = "black")

    lines(supplyList[[2]], supplyList[[1]], type = "l", lty = 3 , col = "blue")
    lines(demandList[[2]], demandList[[1]], type = "l", lty = 3, col = "green")

}


#' Next Layer
#'
#' This calculates the coordinates of the tree nodes in the next layer.
#'
#' @param treeDF A tree in data.frame format
#' @param coordinates The coordinates of the tree nodes as data.frame
#' @param node The current node index
#' @param layer The current layer
#'
#' @return data.frame with coordinates and information for all tree nodes
#'
#' @noRd
#'
nextLayer <- function(treeDF, coordinates ,node, layer){

    children <- treeDF[treeDF$parent ==node,]$child
    numChildren <- length(children)

    if(numChildren == 0){
        return(coordinates)
    }
    # maximum and minimum x coordinate for the next layer
    maxX <- coordinates[coordinates$node == node,]$maxX
    minX <- coordinates[coordinates$node == node,]$minX

    # distance between each child nodes
    distance = (maxX-minX)/numChildren

    # computing the x-cooridantes of the child nodes.
    for(i in 0:(numChildren-1)){
        coordinates[nrow(coordinates)+1,] <- c(children[i+1], maxX-(i+0.5)*distance,
                                               layer, maxX-(i*distance), maxX - (i+1)*distance, node )
        # Compute the next layer for each child.
        coordinates = nextLayer(treeDF, coordinates, children[i+1], layer-1)
    }

    return(coordinates)


}


# #' findPath
# #'
# #' Finding the path between a node and a node in its subtree.
# #'
# #' @param from node at which the path starts
# #' @param to node at which the path ends. Must be in the subtree of 'from'
# #' @param treeDF tree in data.frame format
# #'
# #' @return A list of nodes from 'from' to 'to'
# #' @noRd
# findPath <- function(from, to, treeDF){
# 
#     # if the node has children, check if one of them is the wanted node
#     if(length(treeDF[treeDF$parent == from,]$child) > 0){
#         children <- treeDF[treeDF$parent == from,]$child
# 
#         # if one child is the wanted node, return the current node and the child
#         if(to %in% children){
#             return(c(from,to))
#         }
# 
#         # if non of the children is the wanted node, search in all child nodes
#         for(i in 1:length(children)){
#             path <- findPath(children[i], to, treeDF)
#             # if the node is found in a child node, add the current node to the return list
#             if(!is.null(path)){
#                 return(c(from,path))
#             }
#         }
# 
#         return(NULL)
#     }else{
#         return(NULL)
#     }
# 
# }

#' findPath
#'
#' Finding the path between a node and a node in its subtree.
#'
#' @param from node at which the path starts
#' @param to node at which the path ends. Must be in the subtree of 'from'
#' @param treeDF tree in data.frame format
#'
#' @return A list of nodes from 'from' to 'to'
#' @noRd
findPath <- function(from, to, treeDF){
    
    if(from == to){
        return(from)
    }
    
    
    parent = treeDF[treeDF$child == to,]$parent
    if(identical(parent, numeric(0))){
        return(NULL)
    }
    
    if(parent == from){
        
        return(c(from,to))
        
    }else{
        path <- findPath(from, parent, treeDF)
        if(is.null(path)){
            return(NULL)
        }else{
            return( c(path, to))    
        }
        
    }

}







#' Plotting transport on trees
#'
#'
#' This function visualizes the transport of mass on a tree. 
#'
#'
#' @param tree A tree structure in list format:
#' The first element is the index of the root node.
#' The other elements are vectors defining the edges of the tree. Each of these vectors has to be 
#'  of the form \eqn{(parent_node_index, child_node_index, edge_weight)}.
#' @param tList (optional) The mass transport as list. Each element is a vector of the form
#'  \eqn{(source_node_index, target_node_index, mass)}.
#' @param supply (optional) A non negative numeric vector giving the mass supply at each tree node. The value at 
#' position \eqn{i} give the supply at node \eqn{i}.
#' @param demand (optional) A non negative numeric vector giving the mass demand at each tree node. The value at 
#' position \eqn{i} give the supply at node \eqn{i}.
#'
#'
#' \if{html}{\figure{tree.png}}
#' \if{latex}{\figure{tree.png}{options: width=0.5in}}
#'
#'
#' @examples
#' tree <- list(1, c(1, 2, 1), c(2, 3, 1), c(3, 4, 1), c(3, 5, 1),
#' c(2, 6, 1), c(1, 7, 1), c(7, 8, 1), c(7, 9, 1),
#' c(9, 10, 1), c(9, 11, 1), c(11, 12, 1), c(11, 13, 1))
#'
#' plotTree(tree)
#'
#'
#' supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
#' demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)
#'
#' plotTree(tree, supply = supply, demand = demand)
#'
#'
#' tList = list(c(3,6,1), c(4,8,1))
#' plotTree(tree, tList, supply, demand)
#'
#' @export

plotTree <- function(tree, tList = NULL , supply = NULL, demand = NULL){

    if(length(tList) == 0){
        tList <- NULL
    }

    # create a tree data frame
    treeDF <- as.data.frame(do.call(rbind, tree[-1]))
    colnames(treeDF) <- c("parent", "child", "weight")
    treeDF <- treeDF[order(treeDF$parent),]

    rootNode <- tree[1]



    # initiate a data frame to hold the coordianates for the tree nodes
    coordinates <- data.frame(c(rootNode, 0, 0, 100,-100,-1))
    colnames(coordinates) <- c("node", "x", "layer", "maxX", "minX", "parent")


    # compute all coordinates
    coordinates <- nextLayer(treeDF, coordinates, rootNode, -1)

    maxLayer <- min(coordinates$layer)

    coordinates$layer <- coordinates$layer * 100/(-maxLayer)
    maxLayer <- min(coordinates$layer)


    # add the supply to the coordiantes
    if(!is.null(supply) & !is.null(demand)){
        supDem <- supply-demand
        coordinates$supply <- supDem[coordinates$node]
    }


    # create an empty plot
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(-110, 110), ylim = c(maxLayer-1,1), axes = FALSE)



    # plot the edges
    for(i in 1:nrow(treeDF)){
        segments(coordinates[coordinates$node == treeDF[i,]$parent,]$x,
                 coordinates[coordinates$node == treeDF[i,]$parent,]$layer,
                 coordinates[coordinates$node == treeDF[i,]$child,]$x,
                 coordinates[coordinates$node == treeDF[i,]$child,]$layer)
    }





    # If a transport plan is given, plot the transport paths as arrows
    if(!is.null(tList)){

        treeDF$tMass <- rep(0, nrow(treeDF))

        for(i in 1:length(tList)){
            # Finding a path for each transport list entry
            pathTo <- unlist(findPath(tList[[i]][1], tList[[i]][2], treeDF))
            pathFrom <- NULL

            # If the path was not found, the receiving node is not a child of the
            # origin node
            if(is.null(pathTo)){

                # compute the paths from the trees root node to each of the two nodes.
                # The path between the two nodes can be computed from these two paths.
                pathTo <- unlist(findPath(rootNode, tList[[i]][2], treeDF))
                pathFrom <- unlist(findPath(rootNode, tList[[i]][1], treeDF))

                while(length(pathFrom) > 1 & length(pathTo) > 1 & pathFrom[2] == pathTo[2]){
                    pathTo <- pathTo[-1]
                    pathFrom <- pathFrom[-1]

                }


            }

            if(length(pathTo) > 1){

                for(j in 1:(length(pathTo)-1)){
                    treeDF[(treeDF$parent == pathTo[j] & treeDF$child == pathTo[j+1]),]$tMass =
                        tList[[i]][3] + treeDF[(treeDF$parent == pathTo[j] & treeDF$child == pathTo[j+1]),]$tMass
                }

            }
            if(length(pathFrom) > 1){
                for(j in 1:(length(pathFrom)-1)){
                    treeDF[(treeDF$parent == pathFrom[j] & treeDF$child == pathFrom[j+1]),]$tMass =
                        - tList[[i]][3] + treeDF[(treeDF$parent == pathFrom[j] & treeDF$child == pathFrom[j+1]),]$tMass

                }

            }


        }


        arrowsDF <- treeDF[treeDF$tMass != 0, ]
        # Plotting the arrows.
        # The line width indicates the amount of mass moves along that edge.
        for(i in 1:nrow(arrowsDF)){

            fromX <- coordinates[coordinates$node == arrowsDF[i,]$parent,]$x
            fromY <- coordinates[coordinates$node == arrowsDF[i,]$parent,]$layer
            toX <- coordinates[coordinates$node == arrowsDF[i,]$child,]$x
            toY <- coordinates[coordinates$node == arrowsDF[i,]$child,]$layer

            if(arrowsDF[i,]$tMass > 0){
                curvedarrow(c(fromX, fromY), c(toX, toY), arr.adj = 1, arr.pos = 0.5,arr.type = "triangle", curve = 0.1,
                            lwd = abs(arrowsDF[i,]$tMass), lcol = "red", arr.col = "red")

            }else{

                curvedarrow(c(toX, toY),c(fromX, fromY), arr.adj = 1, arr.pos = 0.5,arr.type = "triangle", curve = 0.1,
                            lwd = abs(arrowsDF[i,]$tMass), lcol = "red", arr.col = "red")

            }



        }
    }


    # If the supply and demand are not given, plot all nodes in black.
    if(is.null(supply) | is.null(demand)){
        points(coordinates$x, coordinates$layer,pch = 19 )

        # Otherwise plot supply nodes in green and demand nodes in blue.
        # The radius of the node indicates the amount of supply / demand: The bigger the node
        # the more mass is supplied or demanded.
        # Nodes without supply or demand are plotted in black.
    }else{
        points(coordinates[coordinates$supply  == 0, ]$x, coordinates[coordinates$supply  == 0, ]$layer ) #,  pch = 19 )

        points(coordinates[coordinates$supply  > 0, ]$x, coordinates[coordinates$supply  > 0, ]$layer,
               pch = 19, cex = abs(coordinates[coordinates$supply > 0, ]$supply),  col = "chartreuse3")

        points(coordinates[coordinates$supply  < 0, ]$x, coordinates[coordinates$supply  < 0, ]$layer,
               pch = 19, cex = abs(coordinates[coordinates$supply  < 0, ]$supply),  col = "dodgerblue3")


    }



    # Adding the keys to the plot.
    text(coordinates$x, coordinates$layer, labels = coordinates$node, pos = 4)

}


