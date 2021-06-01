

#' @title legendre Entropy
#' @param lambda num value
#' @param x num value
#' @param DivFun div fun value
#' @param param1 num val
#' @param param2 num val
#' @noRd
# #' @export
#'
legendre_entropy <- function(lambda, x, DivFun, param1 = 0, param2 = 0){
    
    if(DivFun == "KL"){
        
        return(lambda*(exp(x/lambda) - 1))
        
    }else if(DivFun == "Power"){
        print("le power")
        print(x)
        
        if(param1 == 0){
            return( -lambda *log(1-(x/lambda)))
        }else{
            return(lambda*(1- 1/param1) * ((1+x/(lambda *(param1 - 1)))**param1 - 1) )
        }
        
    }else if(DivFun == "RG"){
        return(pmax(param1 * x, param2*x))
        
    }else if(DivFun == "TV"){
        return(x)
    }
    
    
}

#' @title legendre grad
#' @param lambda num value
#' @param x num value
#' @param DivFun div fun value
#' @param param1 num val
#' @param param2 num val
#'
#' @noRd
grad_legrende <- function(lambda, x, DivFun, param1 = 0, param2 = 0){
    
    if(DivFun == "KL"){
        
        return(exp(x/(lambda)))
        
    }else if(DivFun == "Power"){
        
        return((1-(x/(lambda * (1-param1))))**(param1-1))
        
    }else if(DivFun == "RG"){
        return(pmax(-param1*sign(x), param1 * sign(x)))
        
    }else if(DivFun == "TV"){
        return(rep(1, length(x)))
        
        
    }
    
}

#' @title exp exponent
#' @param f num vec
#' @param g num  vec 
#' @param C cost mat
#'
#' @noRd
expC <- function(f, g, C){
    
    res = matrix(0, length(f), length(g))
    
    for(i in 1:length(f)){
        for(j in 1:length(g)){
            
            res[i,j] <- f[i]+g[j]-C[i,j] 
            
        }
        
    }
    
    return(res)
    
}

#' @title Sink div
#' @param supplyList sup list
#' @param demandList dem list
#' @param eps num val
#' @param iterMax num val
#' @param tol num val
#' @param method cost method
#' @param exp cost exponent
#' @param p cost parameter
#' @param wfr wfr param
#' @param Cxy cost matrix
#' @param Cxx cost matrix
#' @param Cyy cost matrix
#'
#' @export
sinkhorn_divergence <- function(supplyList, demandList, eps, iterMax = 1000, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){

    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    X <- supplyList[[lenSup]]
    Y <- demandList[[lenDem]]
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(X, Y , method, exp, wfr, p)
    }
    
    if(is.null(Cxx)){
        Cxx <- costMatrix(X, X, method, exp, wfr, p)
    }
    
    if(is.null(Cyy)){
        Cyy <- costMatrix(Y, Y , method, exp, wfr, p)
    }
    
    print(Cxy)
    print(Cxx)
    print(Cyy)
    
    if(supplyList[[1]] != demandList[[1]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithmFromCost(Cxy, supplyList, demandList, iterMax, eps, tol)
    res_x = sinkhornAlgorithmFromCost(Cxx, supplyList, supplyList, iterMax, eps, tol)
    res_y = sinkhornAlgorithmFromCost(Cyy, demandList, demandList, iterMax, eps, tol)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
    f_x1 = res_x$dual_f
    f_x2 = res_x$dual_g
    
    g_y1 = res_y$dual_f
    g_y2 = res_y$dual_g
    
    
    print("f_xy ...")
    print(f_xy)
    print(g_xy)
    print(f_x1)
    print(f_x2)
    print(g_y1)
    print(g_y2)
    
    cat("\n\n")
    
    if(supplyList[[1]] == "TV"){
        
        print("TV")
        
        
        print("func")
        print(sum(supplyList[[2]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)))
        print(sum(demandList[[2]] * (g_xy - 0.5*g_y1 -0.5*g_y2)))
        
        
        
        func = sum(supplyList[[2]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) + sum(demandList[[2]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
        
        print(func)
        
        supxdem = supplyList[[2]] %*% t(demandList[[2]])
        supxsup = supplyList[[2]] %*% t(supplyList[[2]])
        demxdem = demandList[[2]] %*% t(demandList[[2]])
        

        func = func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        print((1-exp(expC(f_xy,g_xy,Cxy)/eps)))
        
        print(sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        print(0.5*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))))
        print(0.5*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps))))
        print(eps * sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        print(0.5*eps*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))))
        print(0.5*eps*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps))))
        
        
        print(func)
        
        return(func)
        
        
    }else if(supplyList[[1]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        
        print("outs")
        print(-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2))
        print(0.5 * (-legendre_entropy(0, -f_x1, supplyList[[1]], param1, param2)))
        print(0.5 * (-legendre_entropy(0, -f_x2, supplyList[[1]], param1, param2)))
        print(-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2))
        print(0.5 * (-legendre_entropy(0, -g_y1, supplyList[[1]], param1, param2)))
        print(0.5 * (-legendre_entropy(0, -g_y2, supplyList[[1]], param1, param2)))
        
        
        
        supxdem = supplyList[[2]] %*% t(demandList[[2]])
        supxsup = supplyList[[2]] %*% t(supplyList[[2]])
        demxdem = demandList[[2]] %*% t(demandList[[2]])
        
        
        print("exps")
        
        print(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))
        
        
        print("funcs")
        func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)-
                                          0.5 * (-legendre_entropy(0, -f_x1, supplyList[[1]], param1, param2))-
                                          0.5 * (-legendre_entropy(0, -f_x2, supplyList[[1]], param1, param2)))) +
                                        sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2) -
                                            0.5 * (-legendre_entropy(0, -g_y1, supplyList[[1]], param1, param2)) -
                                            0.5 * (-legendre_entropy(0, -g_y2, supplyList[[1]], param1, param2))))
                  
        print(sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)
                                     - 0.5 * (-legendre_entropy(0, -f_x1, supplyList[[1]], param1, param2))
                                     - 0.5 * (-legendre_entropy(0, -f_x2, supplyList[[1]], param1, param2)))))
        
        print(sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)
                                     - 0.5 * (-legendre_entropy(0, -g_y1, supplyList[[1]], param1, param2))
                                     - 0.5 * (-legendre_entropy(0, -g_y2, supplyList[[1]], param1, param2)))))
          
        print(func)
        
       
    
        
        print(sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        print(0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))))
        print(0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps))))
        
        func = func + sum(eps * supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))) -
            0.5*sum(eps * supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) -
            0.5*sum(eps * demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        print(func)
        return(func)

        
    }else{
        
        param1 <- 0
        
        if(supplyList[[1]] == "Berg"){
            supplyList[[1]] <- "Power"
            demandList[[1]] <- "Power"
            
            param1 <- 0
            
            
        }else if(supplyList[[1]] == "Hellinger"){
            supplyList[[1]] <- "Power"
            demandList[[1]] <- "Power"
            
            param1 <- -1
            
        }else if(supplyList[[1]] == "Power"){
            
            param1 <- supplyList[[4]]
        }
        
        
        
        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[1]], param1)
        outf_xx <- -legendre_entropy(supplyList[[3]], -f_x1, supplyList[[1]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -f_x1, supplyList[[1]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1)
        outg_yy <- -legendre_entropy(supplyList[[3]], -g_y1, supplyList[[1]], param1) - 0.5*eps*grad_legrende(supplyList[[3]], -g_y1, supplyList[[1]], param1)
        
        print("outs")
        print(outf_xy)
        print(outf_xx)
        print(outf_xy-outf_xx)
        print(sum(supplyList[[2]]*(outf_xy-outf_xx)))
        
        
        print(outg_xy)
        print(outg_yy)
        print(outg_xy-outg_yy)
        print(sum(demandList[[2]]*(outg_xy-outg_yy)))
        
        
        out <- sum(supplyList[[2]]*(outf_xy-outf_xx)) + sum(demandList[[2]] * (outg_xy - outg_yy))
        
        return(out)
    }
    
    
    
}





#' @title regularized output
#' @param supplyList sup list
#' @param demandList dem list
#' @param eps num val
#' @param iterMax num val
#' @param tol num val
#' @param method cost method
#' @param exp cost exponent
#' @param p cost parameter
#' @param wfr wfr param
#' @param Cxy cost matrix
#' @export
#' 
regularizedf_ot <- function(supplyList, demandList, eps, iterMax = 100, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL){
    
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
    }
    
    print(Cxy[1:10, 1:10])

    if(supplyList[[1]] != demandList[[1]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithmFromCost(Cxy, supplyList, demandList,
                                       iterMax, eps, tol)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
    
    print("solved")
    print("f_x, g_y: ")
    print(f_xy)
    print(g_xy)
    
    if(supplyList[[1]] == "TV"){
        
        print("funcs")
        print(sum(supplyList[[2]] * f_xy))
        print( sum(demandList[[2]] * g_xy))
        

        func = sum(supplyList[[2]] * f_xy) + sum(demandList[[2]] * g_xy)
        
        print(func)
        
      
        expFun <- supplyList[[2]] %*% t(demandList[[2]]) * (1-exp(expC(f_xy,g_xy,Cxy)/eps))
        
        
        print((1-exp(expC(f_xy,g_xy,Cxy)/eps)))
        print(sum(eps * expFun))
        
        
        func = func + sum(eps * expFun)
        return(func)
        
       
        
    }else if(supplyList[[1]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        supxdem <- supplyList[[2]] %*% t(demandList[[2]])
        
        print("funcs")
        print(sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2))))
        
        print( sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2))))
        
        
        func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)))+ 
            sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)))
        print(func)
        
        print(sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))))
        
        func = func + sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        
        print(func)
        
        return(func)
        

    
    }else{
        
        param1 <- 0
        
        if(supplyList[[1]] == "Berg"){
            supplyList[[1]] <- "Power"
            demandList[[1]] <- "Power"
            
            param1 <- 0
            
            
        }else if(supplyList[[1]] == "Hellinger"){
            supplyList[[1]] <- "Power"
            demandList[[1]] <- "Power"
            
            param1 <- -1
            
        }else if(supplyList[[1]] == "Power"){
            
            param1 <- supplyList[[4]]
        }
        

        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[1]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1)
        
        print("pre outs")
        print(legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, 0))
        print(grad_legrende(supplyList[[3]], -f_xy, supplyList[[1]], param1))
        
        print(legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, 0))
        print(grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1))
        
        print("outs")
        print(outf_xy)
        print(outg_xy)

        print("funcs")
        print(sum(supplyList[[2]] * outf_xy))
        print(sum(supplyList[[2]] * outf_xy)+ sum(demandList[[2]] * outg_xy))
        print(sum(supplyList[[2]] * outf_xy)+ sum(demandList[[2]] * outg_xy)+ eps*(sum(supplyList[[2]]) * sum(demandList[[2]])))
        print(eps*(sum(supplyList[[2]]) + sum(demandList[[2]])))
        res = sum(supplyList[[2]] * outf_xy) + sum(demandList[[2]] * outg_xy) + eps*(sum(supplyList[[2]]) * sum(demandList[[2]]))
        
        return(res)
        
    }

}



#' 
#' #' @title hausdorff div
#' #' @param supplyList sup list
#' #' @param demandList dem list
#' #' @param eps num val
#' #' @param iterMax num val
#' #' @param tol num val
#' #' @param method cost method
#' #' @param exp cost exponent
#' #' @param p cost parameter
#' #' @param wfr wfr param
#' #' @param Cxx cost matrix
#' #' @param Cyy cost matrix
#' #' @param Cxy cost matrix 
#' #' @export
#' #'
#' hausdorff_divergence <- function(supplyList, demandList, eps, iterMax = 1000, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
#'                                 Cxy = NULL, Cxx = NULL, Cyy = NULL){
#' 
#'     lenSup <- length(supplyList)
#'     lenDem <- length(demandList)
#' 
#'     if(is.null(Cxy)){
#'         Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
#'     }
#' 
#'     if(is.null(Cxx)){
#'         Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], method, exp, wfr, p)
#'     }
#' 
#'     if(is.null(Cyy)){
#'         Cyy <- costMatrix(demandList[[lenDem]], demandList[[lenDem]], method, exp, wfr, p)
#'     }
#' 
#' 
#' 
#'     if(supplyList[[1]] != demandList[[1]]){
#' 
#'         print("Please chose the same entropy for supply and demand.")
#' 
#'     }
#' 
#'     
#'     reg <- 0
#'     param1 <- 0
#'     param2 <- 0
#'     
#'     
#' 
#'     
#'     supply <- supplyList[[2]]
#'     demand <- demandList[[2]]
#'     
#' 
#'     res_x = sinkhornAlgorithmFromCost(Cxx, supplyList, supplyList, iterMax, eps, tol)
#'     res_y = sinkhornAlgorithmFromCost(Cyy, demandList, demandList, iterMax, eps, tol)
#' 
#'     f_x = res_x$dual_f
#'     g_y = res_y$dual_g
#' 
#'     if(supplyList[[1]] == "KL"){
#'         div <- 1
#'         reg <- supplyList[[3]]
#'         
#'     }else if(supplyList[[1]] == "TV"){
#'         div <- 2
#'         reg <- supplyList[[3]]
#'     }else if(supplyList[[1]] == "RG"){
#'         div <- 3
#'         param1 <- supplyList[[3]]
#'         param2 <- supplyList[[4]]
#'         
#'         if(param1 < 0 || param2 < param1){
#'             stop("0 <= Alpha <= Beta")
#'         }
#'         
#'     }else if(supplyList[[1]] == "Berg"){
#'         supplyList[[1]] <- "Power"
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- 0
#'     }else if(supplyList[[1]] == "Hellinger"){
#'         supplyList[[1]] <- "Power"
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- -1
#'     }else if(supplyList[[1]] == "Power"){
#'         div <- 4
#'         reg <- supplyList[[3]]
#'         param1 <- supplyList[[4]]
#'     }else{
#'         stop("Please supply a divergence")
#'     }
#'     
#'     
#'     g_xy = Hausdorff_Vec_Rcpp(Cxy, supply, f_x, reg, param1, param2, div, eps)
#'     f_xy = Hausdorff_Vec_Rcpp(t(Cxy), demand, g_y, reg, param1, param2, div, eps)
#'     
#'     
#'     
#'     print("f_x, g_y, ...")
#'     print(f_x)
#'     print(g_y)
#'     print(g_xy)
#'     print(f_xy)
#'     
#'     
#'     print("part outs")
#'     print(legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     print(grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     
#'     print("outs")
#'     
#'     print(legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2))
#'     print(legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_xy,supplyList[[1]], param1, param2))
#'     
#'     print(legendre_entropy(supplyList[[3]], -g_y, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_y,supplyList[[1]], param1, param2))
#'     print(legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2))
#' 
#'     res <- sum(supplyList[[2]] * (legendre_entropy(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -f_x, supplyList[[1]], param1, param2) -
#'                                       legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, param2) - eps * grad_legrende(supplyList[[3]], -f_xy,supplyList[[1]], param1, param2))) +
#'         sum(demandList[[2]] *(legendre_entropy(supplyList[[3]], -g_y, supplyList[[1]], param1, param2) + eps * grad_legrende(supplyList[[3]], -g_y,supplyList[[1]], param1, param2)-
#'                                   legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2) - eps * grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1, param2)))
#' 
#'     
#'     
#'     
#'     return(res)
#' }
