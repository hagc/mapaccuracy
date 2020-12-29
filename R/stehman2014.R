#' Thematic 2 map accuracy and area under stratified random sampling when the strata differ from the map classes.
#'
#' @description Implements the estimators described in Stehman (2014) for overall accuracy, producer's accuracy, user's accuracy, and area. Includes precision estimates.
#'
#' @param s character vector. Strata class labels. The object will be coarsed to factor.
#' @inheritParams olofsson
#' @param Nh_strata numeric vector. Area of the classes in the strata. It must be named (see details).
#'
#' @details Argument \code{Nh_strata} must be named to explicitly and clearly identify the class that each area refers to.
#' The order of \code{Nh} will be used for displaying the results.
#'
#' The labels included in \code{s}, \code{r} and \code{m} must match among them and the names of \code{Nh_strata}.
#'
#' The estimates of area are given in the same unit of the mapped areas (if \code{Nh_strata} express hectares, area estimates express hectares as well).
#'
#' @return A list with the estimates and error matrix.
#' \item{OA}{overall accuracy}
#' \item{UA}{user's accuracy}
#' \item{PA}{producer's accuracy}
#' \item{area}{area}
#' \item{SEoa}{standard error of OA}
#' \item{SEua}{standard error of UA}
#' \item{SEpa}{standard error of PA}
#' \item{SEa}{standard error of area}
#' \item{matrix}{confusion error (area proportion). Rows and columns represent map and reference class labels, respectively}
#'
#' @references
#' Stehman, S. V. (2014).\href{https://doi.org/10.1080/01431161.2014.930207}{Estimating area and map accuracy for stratified random sampling when the strata are different from the map classes}. \emph{Int. J. Remote Sens.}, 35, 4923-4939.
#'
#' @author Hugo Costa
#'
#' @examples
#' # Numerical example in Stheman (2014)
#' s<-c(rep("A",10), rep("B",10), rep("C",10), rep("D",10))
#' m<-c(rep("A",7), rep("B",3), "A", rep("B",11), rep("C",6), "B", "B", rep("D",10))
#' r<-c(rep("A",5), "C", "B", "A", "B", "C", "A", rep("B",5),
#'               "A", "A", "B", "B", rep("C",5), "D", "D", "B", "B", "A", rep("D",7), "C", "C", "B")
#' Nh_strata<-c("A"=40000, "B"=30000, "C"=20000, "D"=10000)
#' e<-stehman2014(s, r, m, Nh_strata)
#' (e$area/sum(Nh_strata))[1] # Proportion of area of class A (compare with paper in p. 4932)
#' (e$area/sum(Nh_strata))[3] # Proportion of area of class C (p. 4932)
#' e$OA                       # Overall accuracy (p. 4932)
#' e$UA[2]                    # User's accuracy of class B (compare with paper in p. 4934)
#' e$PA[2]                    # Producer's accuracy of class B (p. 4934)
#' e$matrix[2,3]              # Cell (2, 3) of the error matrix (p. 4935)
#' e$SEa[1]/Nh_strata[1]      # Standard error (SE) for proportion of area of class A (p. 4935)
#' e$SEa[3]/Nh_strata[3]      # Standard error (SE) for proportion of area of class C (p. 4935)
#' e$SEoa                     # SE for overall accuracy (p. 4936)
#' e$SEua[2]                  # SE for user's accuracy of class B (p. 4936)
#' e$SEpa[2]                  # SE for producer's accuracy of class B (p. 4936)
#'
#' # change class order
#' stehman2014(s, m, r, Nh_strata[c(4,2,1,3)], margins=FALSE)
#' @export
stehman2014<-function(s, r, m, Nh_strata, margins=TRUE){

  # check arguments
  s<-unlist(s)
  r<-unlist(r)
  m<-unlist(m)
  Nh_strata<-unlist(Nh_strata)
  if(is.null(names(Nh_strata))) stop("Nh_strata must be named.", call. = FALSE)
  if(length(names(Nh_strata))>length(unique(names(Nh_strata)))) stop("Repeated names detected in Nh_strata.", call. = FALSE)
  # mapaccuracy:::.check_labels(names(Nh_strata), s, r, m)
  mapaccuracy:::.check_labels(s, names(Nh_strata))
  # mapaccuracy:::.check_labels(r, names(Nh_strata))
  mapaccuracy:::.check_length(s, r, m)

  # covert arguments
  s<-factor(s, levels = names(Nh_strata))
  r<-factor(r, levels = names(Nh_strata))
  m<-factor(m,levels = names(Nh_strata))


  ## PREPARE TABLE 2 OF STEHMAN 2014

  sample<-data.frame(s=s, m=m, r=r)
  classes1<-levels(s)
  classes2<-levels(m)
  classes3<-levels(r)


  # proportion of area of the classes (equation 14)
  q<-length(Nh_strata)
  area_class<-list()
  for (i in 1:q){
    area_class[[i]]<-as.numeric(sample$r==classes1[i])
  }
  A<-Reduce(sum,area_class)

  # comparison between map and reference classes (equation 12)
  O<-as.numeric(sample$m==sample$r)

  # user's accuracy (equations 18 and 19)
  U_class_y<-list()
  U_class_x<-list()
  for (i in 1:q){
    U_class_y[[i]]<-as.numeric(O&sample$r==classes1[i])
    U_class_x[[i]]<-as.numeric(sample$m==classes1[i])
  }

  # producer's accuracy (equation 22 and 23)
  P_class_y<-list()
  P_class_x<-list()
  for (i in 1:q){
    P_class_y[[i]]<-as.numeric(O&sample$r==classes1[i])
    P_class_x[[i]]<-as.numeric(sample$r==classes1[i])
  }

  # "table 2"
  table2<-data.frame(s=sample$s)
  for (i in 1:q){
    table2<-cbind(table2,area_class[[i]])
    names(table2)[length(table2)]<-paste("Area",classes1[i])
  }
  for (i in 1:q){
    table2<-cbind(table2,U_class_y[[i]],U_class_x[[i]])
    n<-length(table2)
    names(table2)[(n-1):n]<-c(paste("U_y",classes1[i],sep="_"),paste("U_x",classes1[i],sep="_"))
  }
  for (i in 1:q){
    table2<-cbind(table2,P_class_y[[i]],P_class_x[[i]])
    n<-length(table2)
    names(table2)[(n-1):n]<-c(paste("P_y",classes1[i],sep="_"),paste("P_x",classes1[i],sep="_"))
  }
  table2<-cbind(table2,O)

  for(i in 1:q){
    for(j in 1:q){
      table2<-cbind(table2,as.numeric(sample$m==classes1[i]&sample$r==classes1[j]))
      n<-length(table2)
      names(table2)[n]<-paste0("p",i,j)
    }
  }


  ## PREPARE TABLE 3  OF STEHMAN 2014
  tmean<-stats::aggregate(.~sample$s, data=table2, mean)
  tvar<-stats::aggregate(.~sample$s, data=table2, stats::var)

  strata<-list()
  covP<-list()
  covU<-list()
  for(i in 1:q){
    a<-which(table2$s==classes1[i])
    strata[[i]]<-table2[a,]
    covP[[i]]<-list()
    covU[[i]]<-list()
    for(j in 1:q){
      covU[[i]][[j]]<-stats::cov(strata[[i]][,paste("U_y",classes1[j],sep="_")],strata[[i]][,paste("U_x",classes1[j],sep="_")])
      covP[[i]][[j]]<-stats::cov(strata[[i]][,paste("P_y",classes1[j],sep="_")],strata[[i]][,paste("P_x",classes1[j],sep="_")])
    }
  }

  # simplify a bit the estruture of the lists covU and covP
  for(i in 1:q){
    covU[[i]]<-unlist(covU[[i]])
    covP[[i]]<-unlist(covP[[i]])
  }
  temp1<-as.data.frame(covU)
  temp2<-as.data.frame(covP)
  temp1<-t(temp1)
  temp2<-t(temp2)
  rownames(temp1)<-NULL
  rownames(temp2)<-NULL


  ## ACCURACY MEASURES

  # Proportion and area of the classes
  Y<-NULL
  for(i in 1:q){
    Y<-c(Y,sum(Nh_strata*tmean[,paste("Area",classes1[i])])/sum(Nh_strata))
  }
  area<-Y*sum(Nh_strata)
  names(area)<-names(Nh_strata)

  # Overall Accuracy
  O<-sum(Nh_strata*tmean[,"O"])/sum(Nh_strata)

  # User's accuracy
  U <- R1 <- R2 <- NULL
  for(i in 1:q){
    R1<-c(R1,sum(Nh_strata*tmean[,paste("U_y",classes1[i],sep="_")]))
    R2<-c(R2,sum(Nh_strata*tmean[,paste("U_x",classes1[i],sep="_")]))
    U<-c(U,R1[i]/R2[i])
  }
  names(U)<-names(Nh_strata)

  # Producer's accuracy
  P <- R3 <- R4 <- NULL
  for(i in 1:q){
    R3<-c(R3,sum(Nh_strata*tmean[,paste("P_y",classes1[i],sep="_")]))
    R4<-c(R4,sum(Nh_strata*tmean[,paste("P_x",classes1[i],sep="_")]))
    P<-c(P,R3[i]/R4[i])
  }
  names(P)<-names(Nh_strata)

  # Cell i,j of the error matrix
  matrix<-matrix(nrow=q,ncol=q)
  for(i in 1:q){
    for(j in 1:q){
      matrix[i,j]<-sum(Nh_strata*tmean[,paste0("p",i,j)])/sum(Nh_strata)
    }
  }

  # Standard error (SE) for proportion
  Va<-NULL
  for(i in 1:q){
    Va<-c(Va,1/sum(Nh_strata)^2*sum(Nh_strata^2*tvar[,paste("Area",classes1[i])]/table(sample$s)[i])) # This was 10, but I changed to table(sample$s)[i]
  }
  SEa<-sqrt(Va)
  SEa<-Nh_strata*SEa
  names(SEa)<-names(Nh_strata)

  # Standard error (SE) for overall accuracy
  Vo<-1/sum(Nh_strata)^2*sum(Nh_strata^2*tvar[,"O"]/table(sample$s)) # This was 10, but I changed to sample$s)
  SEo<-sqrt(Vo)

  # Standard error (SE) for user's accuracy
  Vu<-NULL
  for(i in 1:q){
    Vu<-c(Vu,1/R2[i]^2*sum(Nh_strata^2*(tvar[,paste("U_y",classes1[i],sep="_")]+U[i]^2*tvar[,paste("U_x",classes1[i],sep="_")]-2*U[i]*temp1[,i])/table(sample$s)[i])) # This was 10, but I changed to npixels
  }
  SEu<-sqrt(Vu)
  names(SEu)<-names(Nh_strata)

  # Standard error (SE) for producer's accuracy
  Vp<-NULL
  for(i in 1:q){
    Vp<-c(Vp,1/R4[i]^2*sum(Nh_strata^2*(tvar[,paste("P_y",classes1[i],sep="_")]+P[i]^2*tvar[,paste("P_x",classes1[i],sep="_")]-2*P[i]*temp2[,i])/table(sample$s)[i])) # This was 10, but I changed to npixels
  }
  SEp<-sqrt(Vp)
  names(SEp)<-names(Nh_strata)

  # Add margins to confusion matrix
  colnames(matrix)<-classes1; rownames(matrix)<-classes1
  if(margins){
    matrix<-cbind(matrix,apply(matrix,1,sum))
    matrix<-rbind(matrix,c(apply(matrix,2,sum)))
    colnames(matrix)[length(colnames(matrix))]<-"sum"; rownames(matrix)[length(rownames(matrix))]<-"sum"
  }


  # return
  list(OA=O,
       UA=U,
       PA=P,
       area=area,
       SEoa=SEo,
       SEua=SEu,
       SEpa=SEp,
       SEa=SEa,
       matrix=matrix)
}
