#' @title Thematic map accuracy for fuzzy classification/reference data under simple random sampling.
#'
#' @description Implements the estimators described in Binaghi \emph{et al}. (1999) for overall accuracy, producer's accuracy, user's accuracy, and index of fuzziness. Includes \emph{no} precision estimates.
#'
#' @param r dataframe or matrix. Reference data (rows relate to pixels/points and columns to class membership degrees). The columns must be named (see details).
#' @param m dataframe or matrix. Membership degrees of a soft classification (rows relate to pixels/points and columns to class membership degrees). The columns must be named (see details).
#' @inheritParams olofsson
#'
#' @details The columns of both arguments must be named to explicitly and clearly identify the class that each column refers to. The columns may be sorted differently, and the order of \code{r} will be used for displaying the results.
#'
#' @return A list with the estimates and matrices, and a csv if defined.
#' \item{OA}{overall accuracy}
#' \item{UA}{user's accuracy}
#' \item{PA}{producer's accuracy}
#' \item{matrix}{fuzzy confusion error. Rows and columns represent map and reference class labels, respectively.}
#' \item{numeric}{index of fuzziness of \code{r}}
#' \item{numeric}{index of fuzziness of \code{m}}
#'
#' @references
#' Binaghi, E.; Brivio, P. A.; Ghezzi, P.; Rampini, A. (1999). \href{https://doi.org/10.1016/S0167-8655(99)00061-6}{A fuzzy set-based accuracy assessment of soft classification}. \emph{Pattern Recognit. Lett.}, 20, 935-948.
#'
#' @author Hugo Costa
#'
#' @seealso \code{\link{stehman2007}}
#'
#' @examples
#' ## Example of table 2 in Binaghi et al. (1999), p. 940.
#' # case a (perfect matching):
#' r<-matrix(c(0.4,0.4,0.4), ncol = 3)
#' m<-matrix(c(0.4,0.4,0.4), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' binaghi(r,m) # compare to paper
#'
#' # case b (under-estimation)
#' r<-matrix(c(0.4,0.4,0.4), ncol = 3)
#' m<-matrix(c(0.2,0.4,0.4), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' binaghi(r,m) # compare to paper
#'
#' # case c (over-estimation):
#' r<-matrix(c(0.4,0.4,0.4), ncol = 3)
#' m<-matrix(c(0.6,0.4,0.4), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' binaghi(r,m) # compare to paper
#'
#' ## Example of table 3 in Binaghi et al. (1999), p. 941.
#' # case a (perfect matching):
#' r<-matrix(c(0.7,0.2,0.1), ncol = 3)
#' m<-matrix(c(0.7,0.2,0.1), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' binaghi(r,m) # compare to paper
#'
#' # case b (under-estimation)
#' r<-matrix(c(0.7,0.2,0.1), ncol = 3)
#' m<-matrix(c(0.6,0.3,0.1), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' binaghi(r,m) # compare to paper
#' # (note the error in the paper: UA and PA are swiped)
#'
#' # case c (over-estimation):
#' r<-matrix(c(0.7,0.2,0.1), ncol = 3)
#' m<-matrix(c(0.8,0.1,0.1), ncol = 3)
#' colnames(r) <- colnames(m) <- c("q1", "q2", "q3")
#' e<-binaghi(r,m) # compare to paper
#' # (note the error in the paper: UA and PA are swiped)
#'
#' ## additional example
#' #m<-readRDS("D:/mypapers/11-mapaccuracy/data/fcla.rds")
#' #r<-readRDS("D:/mypapers/11-mapaccuracy/data/fref.rds")
#' data(fcla)
#' data(fref)
#' binaghi(fref,fcla)
#' @export
binaghi<-function(r, m, margins=TRUE){

  # check arguments
  r<-as.matrix(r)
  m<-as.matrix(m)
  if(is.null(colnames(r))) stop("The columns of r must be named.", call. = FALSE)
  if(is.null(colnames(m))) stop("The columns of m must be named.", call. = FALSE)
  if(length(colnames(r))>length(unique(colnames(r)))) stop("Repeated names detected in argument r.", call. = FALSE)
  if(length(colnames(m))>length(unique(colnames(m)))) stop("Repeated names detected in argument m.", call. = FALSE)
  if(nrow(r)!=nrow(m)){stop("Number of rows in 'r' and 'm' differ.", call. = FALSE)}

  virt<-matrix(data=1:ncol(r), ncol=ncol(r))
  colnames(virt)<-colnames(r)
  m<-gtools::smartbind(virt, m, fill=0)
  m<-m[-1,]
  m<-as.matrix(m)

  # # sort m and add missing columns if any
  # if(ncol(r)>ncol(m)){
  #   virt<-matrix(data=1:ncol(r), ncol=ncol(r))
  #   colnames(virt)<-colnames(r)
  #   m<-gtools::smartbind(virt, m, fill=0)
  #   m<-m[-1,]
  # }else if(nrow(m)>1){
  #   m<-m[,colnames(r)]  # This line destroys the matrix if nrow=1
  # }



  #fuzzy matrix ("min" fuzzy set operator used)
  minop<-NULL
  for (i in 1:nrow(m)){
    for(j in 1:ncol(m)){
      for(k in 1:ncol(m)){
        minop<-rbind(minop,c(colnames(m)[j], colnames(r)[k], min(as.numeric(m[i,j]),as.numeric(r[i,k]))))
      }
    }
  }
  minop<-as.data.frame(minop,stringsAsFactors=F)
  minop[,3]<-as.numeric(minop[,3])
  minop<-reshape::cast(minop, V1~V2, fun.aggregate=sum, value="V3", margins=FALSE)
  rownames(minop)<-minop[,1]
  minop<-minop[,-1]

  # if(!missin(weights)){
  #   minop<-minop*weights
  # }

  #Margins of the fuzzy matrix
  mref<-apply(r,2,sum)
  mcla<-apply(m,2,sum)

  #Accuracy measures
  O<-sum(diag(as.matrix(minop)))/sum(mref); names(O)<-""
  U<-diag(as.matrix(minop))/mcla
  P<-diag(as.matrix(minop))/mref

  # U<-t(as.data.frame(U))
  # P<-t(as.data.frame(P))

  #Add margins to confusion matrix
  if(margins){
    minop<-cbind(minop,mcla)
    minop<-rbind(minop,c(mref,sum(mref)))
    names(minop)[length(names(minop))]<-"sum"; rownames(minop)[length(names(minop))]<-"sum"
  }


  # index of fuzziness
  indexf<-function(x){
    y<-x>0.5
    mean(apply(abs(x-y),2,sum)/apply(x,2,sum))
  }
  IFr<-indexf(r)
  IFc<-indexf(m)

  # return
  list(OA=O,
       UA=U,
       PA=P,
       matrix=minop,
       IFr=IFr,
       IFc=IFc)

}
