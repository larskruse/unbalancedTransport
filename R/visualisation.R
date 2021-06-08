#' Plotting transport between single points.
#'
#' Visualizing mass transport between point clouds. The coordinates can either be given directly or calculated from a given
#' distance matrix. 
#'
#' @param transportPlan A non negative numeric matrix that indicates where the mass is
#' transported. The value at point \eqn{\[i,j\]} is the amount of mass transported from
#' supply point \eqn{i} to demand point \eqn{j}.
#' @param supply (optional) A non negative numeric vector giving the mass supply at each point. The value at 
#' position \eqn{i} give the supply at point \eqn{i}.
#' @param demand (optional) A non negative numeric vector giving the mass demand at each point. The value at 
#' position \eqn{i} give the supply at point \eqn{i}.
#' @param creationDestructionCost (optional) A numeric vector giving the cost for mass creation and destruction.
#' @param distanceMatrix (optional) A non negative matrix containing the distances between all points. It is used
#' to calcuated the positions of the points if they are not explicitly given.
#'
#' @export
plotTransportByCost <- function(transportPlan, supply = NULL, demand = NULL, distanceMatrix = NULL,  creationDestructionCost = rep(0, length(x))){

    # calculate the positions
        
    positions <- positionsByEigenvalue(distanceMatrix)    
        
    
    
    positions <- positionsByEigenvalue(distanceMatrix)
    print("Coordinates:")
    print(positions)

    # if the positions are in 1 dimension, add a second by using a 0L vector as
    # y coordinates
    if(is.null(dim(positions))){
        x <- positions
        y <- rep(0, length(positions))

    # otherwise, use the projection on the first two dimensions
    }else{
        x <- positions[,1]
        y <- positions[,2]
    }
    dim <- length(x)

    if(!is.null(supply) & !is.null(demand)){
        supDem <- supply-demand
    }

    # create an empty plot
    plot(1, type = "n", xlab = "", ylab = "",
         xlim = c(min(x)-0.3*abs(max(x)-min(x)),max(x)+0.3*abs(max(x)-min(x))),
         ylim =c(min(y)-0.3*abs(max(y)-min(y)),max(y)+0.3*abs(max(y)-min(y))),
         asp = 1)


    # if a transport plan is given, add arrows to indicate mass transport
    if(!is.null(transportPlan)){

        # plot circles as indicators for the creation/destruction cost
        theta = seq(0, 2 * pi, length = 200)
        for (i in 1:dim){
            lines(x = creationDestructionCost[i] * cos(theta) + x[i], y = creationDestructionCost[i] * sin(theta) + y[i], col = "grey")

            for (j in 1:dim){
                if (transportPlan[[i,j]]>0 && i != j){
                    curvedarrow(c(x[i],y[i]),c(x[j],y[j]), lwd = 0.5+transportPlan[[i,j]],
                                         arr.pos = 0.5, arr.adj = 0.5, arr.type = "triangle",
                                         curve = 0.2)
                }
            }
        }
    }

    # add the points
    if(is.null(supply) | is.null(demand)){
        points(x, y ,pch = 19 )


    }else{

        # if supply and demand values are given, the color indicates the supply
        # and the size indicates the mass at each point

        points(x[which(supDem == 0 )],y[which(supDem == 0)],  pch = 19 )

        points(x[which(supDem > 0 )],y[which(supDem > 0)],
               pch = 19, cex = supDem[which(supDem > 0)],  col = "chartreuse3")

        points(x[which(supDem < 0 )],y[which(supDem < 0)],
               pch = 19, cex = abs(supDem[which(supDem < 0)]),  col = "dodgerblue3")


    }

}


#' A grid plot of the transport plan
#'
#' This function creates a grid plot of a transport plan. Import and export vectors can
#' be given as additional arguments. These will be plotted along the sides of the grid 
#' plot in order to indicate where mass is added or removed. Hence, the mass shown in
#' each row or column is equal to the supply or demand vectors. 
#'
#' @param transportPlan A non negative numeric matrix that indicates where the mass is
#' transported. The value at point \eqn{\[i,j\]} is the amount of mass transported from
#' supply point \eqn{i} to demand point \eqn{j}.
#' @param import (optional) A non negative numeric vector that give the amount of mass created at each
#' demand point. It length has to be equal to the number of columns in the transport matrix.  
#' @param export (optional) A non negative numeric vector that give the amount of mass destroyed  at each
#' supply point. It length has to be equal to the number of rows in the transport matrix.  
#'
#' @examples
#' 
#' transport <- matrix(runif(9), nrow = 3)
#' import <- runif(3)
#' export <- runif(3)
#' 
#' gridPlotTransport(transport)
#' gridPlotTransport(transport, import, export)
#'
#' @export
gridPlotTransport <- function(transportPlan, import = NULL, export =  NULL){

    # If no import or export vector is given, only the transport plan is plotted
    if(is.null(import) | is.null(export)){


        transportPlan <- t(transportPlan[(nrow(transportPlan)):1,])

        image(transportPlan, asp = 1, axes = FALSE, ylab = "Supply", xlab = "Demand",
              col=hcl.colors(20, palette = "viridis", alpha = NULL, rev = TRUE, fixup = TRUE))


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
#' @param transportPlan A non negative numeric matrix giving the mass transport. The value at \eqn{(i,j)}
#' gives the mass transported from supply point \eqn{i} to demand point \eqn{j}.
#' @param supplyList A list containing the non negative supply measure and the underlying discretization as vectors.
#' @param demandList A list containing the non negative demand measure and the underlying discretization as vectors.
#'
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
    lines(demandList[[2]], demandList[[1]], type = "l", lty = 3, col = "red")

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

        polygon(c(firstDem,demandList[[2]]),c(0,t(subK) %*% x1Measure), col = colors[i])
        polygon(c(firstSupp,supplyList[[2]]), c(0,subK %*% y1Measure), col = colors[i])

        #lines(demandList[[2]], t(subK) %*% x1Measure, type = "l", col = "black")
        #lines(supplyList[[2]], subK %*% y1Measure, type = "l", col = "black")

    }

    lines(demandList[[2]], t(transportPlan) %*% x1Measure, type = "l", col = "red")
    lines(supplyList[[2]], transportPlan %*% y1Measure, type = "l", col = "blue")
    # 
    lines(supplyList[[2]], rep(0, length(supplyList[[2]])), type = "l", col = "black")


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

    # if the node has children, check if one of them is the wanted node
    if(length(treeDF[treeDF$parent == from,]$child) > 0){
        children <- treeDF[treeDF$parent == from,]$child

        # if one child is the wanted node, return the current node and the child
        if(to %in% children){
            return(c(from,to))
        }

        # if non of the children is the wanted node, search in all child nodes
        for(i in 1:length(children)){
            path <- findPath(children[i], to, treeDF)
            # if the node is found in a child node, add the current node to the return list
            if(!is.null(path)){
                return(c(from,path))
            }
        }

        return(NULL)
    }else{
        return(NULL)
    }

}


#' Plotting transport on trees
#'
#' This function visualizes the transport of mass on a tree. 
#'
#' @param tree A tree structure in list format: The first element is the index of the root node
#'  followed by multiple vectors defining the edges of the tree. Each of these vectors has to be 
#'  of the form \eqn{(parent, child, weight)}.
#' @param tList (optional) The mass transport as list. Each element is a vector of the form
#'  \eqn{(source, target, mass)}.
#' @param supply (optional) A non negative numeric vector giving the mass supply at each tree node. The value at 
#' position \eqn{i} give the supply at node \eqn{i}.
#' @param demand (optional) A non negative numeric vector giving the mass demand at each tree node. The value at 
#' position \eqn{i} give the supply at node \eqn{i}.
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


