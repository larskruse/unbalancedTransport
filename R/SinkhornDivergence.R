legendre_entropy <- function(lambda, x, DivFun, param1 = 0, param2 = 0){
    
    if(DivFun == "KL"){
        
        return(lambda*exp(x/lambda) - 1)
        
    }else if(DivFun == "Power"){
        
        if(param1 == 0){
            return( -lambda *(1-log(x/lambda)))
        }else{
            return(lambda*(1- 1/param1) * ((1+x/(lambda *(param1 - 1)))**param1 - 1) )
        }
        
    }else if(DivFun == "RG"){
        return(pmax(param1 * x, param2*x))
        
    }else if(DivFun == "TV"){
        return(x)
    }
    
    
}


grad_legrende <- function(lambda, x, DivFun, param1 = 0){
    
    if(DivFun == "KL"){
        
        return(exp(x/(lambda)))
        
    }else{
        
        return((1-(x/(lambda * (1-param1))))**(param1-1))
        
    }
    
}


expC <- function(f, g, C){
    
    res = matrix(0, length(f), length(g))
    
    for(i in 1:length(f)){
        for(j in 1:length(g)){
            
            res[i,j] <- f[i]+g[j]-C[i,j] 
            
        }
        
    }
    
    return(res)
    
}


sinkhorn_divergence <- function(supplyList, demandList, eps, iterMax = 100, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){

    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
    }
    
    if(is.null(Cxx)){
        Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], method, exp, wfr, p)
    }
    
    if(is.null(Cyy)){
        Cyy <- costMatrix(demandList[[lendem]], supplyList[[lendem]], method, exp, wfr, p)
    }
    
    
    
    if(supplyList[[1]] != demandList[[1]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithmFromCost(Cxy, supplyList, demandList,
                                       iterMax, eps, tol)
    res_x = sinkhornAlgorithmFromCost(Cxx, supplyList, supplyList,
                                      iterMax, eps, tol)
    res_y = sinkhornAlgorithmFromCost(Cyy, demandList, demandList,
                                      iterMax, eps, tol)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
    f_x1 = res_x$dual_f
    f_x2 = res_x$dual_g
    
    g_y1 = res_y$dual_f
    g_y2 = res_y$dual_g
    
    
    if(supplyList[[1]] == "TV"){
        
        func = sum(supplyList[[2]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) + sum(demandList[[2]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
        
        supxdem = supplyList[[2]] %*% demandList[[2]]
        supxsup = supplyList[[2]] %*% supplyList[[2]] 
        demxdem = demandList[[2]] %*% demandList[[2]]
        

        func = func + sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))
                - 0.5*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) 
                - 0.5*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)
        
        
    }else if(supplyList[[1]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        -legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)
        

        func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)
                    - 0.5 * (-legendre_entropy(0, -f_x1, supplyList[[1]], param1, param2))
                    - 0.5 * (-legendre_entropy(0, -f_x2, supplyList[[1]], param1, param2))))
                + sum(demadnList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)
                    - 0.5 * (-legendre_entropy(0, -g_y1, supplyList[[1]], param1, param2))
                    - 0.5 * (-legendre_entropy(0, -g_y2, supplyList[[1]], param1, param2))))
                    
                    
        
        supxdem = supplyList[[2]] %*% demandList[[2]]
        supxsup = supplyList[[2]] %*% supplyList[[2]] 
        demxdem = demandList[[2]] %*% demandList[[2]]
        
        
        func = func + sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))
        - 0.5*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) 
        - 0.5*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
        
        return(func)

        
    }else{
        
        outf_xy <- -legendre_entropy(-f_xy) - 0.5*eps*grad_legrende(-f_xy)
        outf_xx <- -legendre_entropy(-f_x1) - 0.5*eps*grad_legrende(-f_x)
        outg_xy <- -legendre_entropy(-g_xy) - 0.5*eps*grad_legrende(-g_xy)
        outg_yy <- -legendre_entropy(-g_y1) - 0.5*eps*grad_legrende(-g_y)
        
        out <- sum(supplyList[[2]]*(outf_xx-outf_xy)) + sum(demandList[[2]] * (outg_yy - outg_xy))
        
        return(out)
    }
    
    
    
}






regularizedf_ot <- function(supplyList, demandList, eps, iterMax = 100, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
                            Cxy = NULL, Cxx = NULL, Cyy = NULL){
    
    lenSup <- length(supplyList)
    lenDem <- length(demandList)
    
    if(is.null(Cxy)){
        Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
    }

    if(supplyList[[1]] != demandList[[1]]){
        
        print("Please chose the same entropy for supply and demand.")
        
    }
    
    res_xy = sinkhornAlgorithmFromCost(Cxy, supplyList, demandList,
                                       iterMax, eps, tol)
    
    f_xy = res_xy$dual_f
    g_xy = res_xy$dual_g
    
   
    
    if(supplyList[[1]] == "TV"){

        func = sum(supplyList[[2]], f_xy) + sum(demandList[[2]], g_xy)
        
        expFun <- supplyList[[2]] %*% demandList[[2]] * (1-exp(expC(f_xy,g_xy,Cxy)/eps))
        
        func = func + sum(eps * expFun)
        return(func)
        
       
        
    }else if(supplyList[[1]] == "RG" ){
        
        param1 <- supplyList[[3]]
        param2 <- supplyList[[4]]
        
        supxdem <- supplyList[[2]] %*% demandList[[2]]
        
        -legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)     
        

        func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)))
                   + sum(demandList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)))

        func = func + sum(eps * (supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps))))
        
        return(func)
        

    
    }else{
        
        if(supplyList[[1]] == "Power"){
            
            param1 = supplyList[[4]]
            
        }else{
            param1 = 0
        }
        
    
        outf_xy <- -legendre_entropy(supplyList[[3]], -f_xy, supplyList[[1]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -f_xy, supplyList[[1]], param1)
        outg_xy <- -legendre_entropy(supplyList[[3]], -g_xy, supplyList[[1]], param1, 0) - 0.5*eps*grad_legrende(supplyList[[3]], -g_xy, supplyList[[1]], param1)
        
        res = sum(supplyList[[2]] * outf_xy) + sum(demandList[[2]] * outg_xy) + eps*(sum(supplyList[[2]]) + sum(demandList[[2]]))
        
        return(res)
        
    }

}





# hausdorff_divergence <- function(supplyList, demandList, eps, iterMax = 100, tol = 1e-3, method = "euclidean", exp = 1, p = 2,  wfr = FALSE,
#                                 Cxy = NULL, Cxx = NULL, Cyy = NULL){
#     
#     lenSup <- length(supplyList)
#     lenDem <- length(demandList)
#     
#     if(is.null(Cxy)){
#         Cxy <- costMatrix(supplyList[[lenSup]], demandList[[lenDem]], method, exp, wfr, p)
#     }
#     
#     if(is.null(Cxx)){
#         Cxx <- costMatrix(supplyList[[lenSup]], supplyList[[lenSup]], method, exp, wfr, p)
#     }
#     
#     if(is.null(Cyy)){
#         Cyy <- costMatrix(demandList[[lendem]], supplyList[[lendem]], method, exp, wfr, p)
#     }
#     
#     
#     
#     if(supplyList[[1]] != demandList[[1]]){
#         
#         print("Please chose the same entropy for supply and demand.")
#         
#     }
#     
#     res_xy = sinkhornAlgorithmFromCost(Cxy, supplyList, demandList,
#                                        iterMax, eps, tol)
#     res_x = sinkhornAlgorithmFromCost(Cxx, supplyList, supplyList,
#                                       iterMax, eps, tol)
#     res_y = sinkhornAlgorithmFromCost(Cyy, demandList, demandList,
#                                       iterMax, eps, tol)
#     
#     f_xy = res_xy$dual_f
#     g_xy = res_xy$dual_g
#     
#     f_x1 = res_x$dual_f
#     f_x2 = res_x$dual_g
#     
#     g_y1 = res_y$dual_f
#     g_y2 = res_y$dual_g
#     
#     
#     if(supplyList[[1]] == "TV"){
#         
#         func = sum(supplyList[[2]] * (f_xy - 0.5*f_x1 - 0.5*f_x2)) + sum(demandList[[2]] * (g_xy - 0.5*g_y1 -0.5*g_y2))
#         
#         supxdem = supplyList[[2]] %*% demandList[[2]]
#         supxsup = supplyList[[2]] %*% supplyList[[2]] 
#         demxdem = demandList[[2]] %*% demandList[[2]]
#         
#         
#         func = func + sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))
#         - 0.5*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) 
#         - 0.5*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
#         
#         return(func)
#         
#         
#     }else if(supplyList[[1]] == "RG" ){
#         
#         param1 <- supplyList[[3]]
#         param2 <- supplyList[[4]]
#         
#         -legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)
#         
#         
#         func = sum(supplyList[[2]] * (-legendre_entropy(0, -f_xy, supplyList[[1]], param1, param2)
#                                       - 0.5 * (-legendre_entropy(0, -f_x1, supplyList[[1]], param1, param2))
#                                       - 0.5 * (-legendre_entropy(0, -f_x2, supplyList[[1]], param1, param2))))
#         + sum(demadnList[[2]] * (-legendre_entropy(0, -g_xy, supplyList[[1]], param1, param2)
#                                  - 0.5 * (-legendre_entropy(0, -g_y1, supplyList[[1]], param1, param2))
#                                  - 0.5 * (-legendre_entropy(0, -g_y2, supplyList[[1]], param1, param2))))
#         
#         
#         
#         supxdem = supplyList[[2]] %*% demandList[[2]]
#         supxsup = supplyList[[2]] %*% supplyList[[2]] 
#         demxdem = demandList[[2]] %*% demandList[[2]]
#         
#         
#         func = func + sum(supxdem * (1-exp(expC(f_xy,g_xy,Cxy)/eps)))
#         - 0.5*sum(supxsup * (1-exp(expC(f_x1,f_x2,Cxx)/eps))) 
#         - 0.5*sum(demxdem * (1-exp(expC(g_y1,g_y2,Cyy)/eps)))
#         
#         return(func)
#         
#         
#     }else{
#         
#         outf_xy <- -legendre_entropy(-f_xy) - 0.5*eps*grad_legrende(-f_xy)
#         outf_xx <- -legendre_entropy(-f_x1) - 0.5*eps*grad_legrende(-f_x)
#         outg_xy <- -legendre_entropy(-g_xy) - 0.5*eps*grad_legrende(-g_xy)
#         outg_yy <- -legendre_entropy(-g_y1) - 0.5*eps*grad_legrende(-g_y)
#         
#         out <- sum(supplyList[[2]]*(outf_xx-outf_xy)) + sum(demandList[[2]] * (outg_yy - outg_xy))
#         
#         return(out)
#     }
#     
#     
#     
# }
