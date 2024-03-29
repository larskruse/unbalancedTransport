% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umtp.R
\name{umtp}
\alias{umtp}
\title{Solver for the unbalanced mass transportation problem}
\usage{
umtp(supplyList, demandList, exp = 1, p = 2, wfr = FALSE, costMatrix = NULL)
}
\arguments{
\item{supplyList}{A list containing the information about the supply. The first element hast to be the distribution followed by
a vector specifying the cost for destruction of mass at each supply point. If the cost matrix
is not provided, the third element has to be the positions of the supply points.
This can either be a vector or a matrix where each row gives the coordinates for one point.}

\item{demandList}{A list similar to the supplyList but holding the information about the demand distribution.}

\item{exp}{(optional) The exponent applied to the cost function. Default value is 1.}

\item{p}{(optional) The parameter for calculating the L_p cost. Can be a positive real
number or Inf. Default value is 2 calculating the euclidean distance}

\item{wfr}{(optional) Computes the cost matrix needed for the Wasserstein-Fisher-Rao
 distance \eqn{c(x,y) = -\log(\cos^2_+(d(x,y)))}{c(x,y) = -log(cos_+(d(x,y)²))}.
The default value is "false".}

\item{costMatrix}{(optional) Instead of having the algorithm compute the
cost matrix, a custom cost matrix can be passed to the algorithm.}
}
\description{
Solving the unbalanced mass transportation problem by extending to a balanced
problem and applying the networkflow algorithm.
}
\details{
This function solves unbalanced optimal transport problems by extending them
to a balanced problem and applying the networkflow algorithm.
The unbalanced mass transportation problem \eqn{<C,r> + \sum_i p_i(\alpha_i - \sum_j r_{ij}) + \sum_j q_j(\beta-\sum_i r{ij})}{
<C,r> + sum_i(p_i (a_i-sum_j r_ij)) + sum_j(q_j (b_j-sum_i r_ij))} with supply
and demand measure \eqn{\alpha}{a} and \eqn{\beta}{b}, cost matrix C, 
transport plan \eqn{r} and construction and destruction costs \eqn{p} and 
\eqn{q} is transformed to a balanced transport problem by extending the cost
 matrix as well as supply and demand vectors. 

The resulting transport problem is then solved using the networkflow algorithm provided by the \code{\link[transport]{transport}} package.
}
\examples{
supplyPoints <- matrix(c(0,0,0,0,1,2), ncol = 2)
demandPoints <- matrix(c(3,3.5,4,1.7,1,0.2), ncol = 2)


p <- c(2,1,2)
q <- c(1,2,3)

costCreate <- rep(10,3)
costDestruct <- rep(10,3)

supplyList <- list(p,costDestruct, supplyPoints)
demandList <- list(q,costCreate, demandPoints)


umtp(supplyList, demandList, exp = 2)


}
\references{
\insertRef{Guittet2002}{unbalancedTransport}
}
