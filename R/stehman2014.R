#' Thematic map accuracy and area under stratified random sampling when the strata differ from 
#' the map classes.
#'
#' Implements the estimators described in Stehman (2014) for overall accuracy, producer's 
#' accuracy, user's accuracy, and area. Includes precision estimates.
#' 
#' Argument \code{Nh_strata} must be named to explicitly and clearly identify the stratum that 
#' each entry refers to.
#' 
#' In the error matrix returned, the entries corresponding to no observed cases will present
#' \code{NA} rather than 0. This is to emphasize the difference between the absence of cases
#' and the presence of some (few) cases that represent a very small proportion of area (almost
#' zero) and thus possibly rounded to zero. However, \code{NA} means zero proportion of area.
#'
#' @param s character vector. Strata class labels. The object will be coerced to factor.
#' @inheritParams olofsson
#' @param Nh_strata numeric vector. Number of pixels forming each stratum. It must be named (see details).
#' @param order character vector. Order of the classes to be displayed in the results.
#'
#' @return A list with the estimates and error matrix.
#' \item{OA}{overall accuracy}
#' \item{UA}{user's accuracy}
#' \item{PA}{producer's accuracy}
#' \item{area}{area proportion}
#' \item{SEoa}{standard error of OA}
#' \item{SEua}{standard error of UA}
#' \item{SEpa}{standard error of PA}
#' \item{SEa}{standard error of area proportion}
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
#' e$area[1]                  # Proportion of area of class A (compare with paper in p. 4932)
#' e$area[3]                  # Proportion of area of class C (p. 4932)
#' e$OA                       # Overall accuracy (p. 4932)
#' e$UA[2]                    # User's accuracy of class B (compare with paper in p. 4934)
#' e$PA[2]                    # Producer's accuracy of class B (p. 4934)
#' e$matrix[2,3]              # Cell (2, 3) of the error matrix (p. 4935)
#' e$SEa[1]                   # Standard error (SE) for proportion of area of class A (p. 4935)
#' e$SEa[3]                   # Standard error (SE) for proportion of area of class C (p. 4935)
#' e$SEoa                     # SE for overall accuracy (p. 4936)
#' e$SEua[2]                  # SE for user's accuracy of class B (p. 4936)
#' e$SEpa[2]                  # SE for producer's accuracy of class B (p. 4936)
#' 
#'
#' # change class order
#' stehman2014(s, m, r, Nh_strata, order=c("D","C","B","A"))
#' 
#' 
#' # When the number of strata differs from the number of map classes
#' s<-c(rep("A",5), rep("AA",5), rep("B",10), rep("C",10), rep("D",10))
#' m<-c(rep("A",4), "B", rep("A",4), "B", "A", rep("B",11), rep("C",6), "B", "B", rep("D",10))
#' r<-c(rep("A",4), "C", rep("A",4), "C", "A", rep("B",5),
#'      "A", "A", "B", "B", rep("C",5), "D", "D", "B", "B", "A", rep("D",7), "C", "C", "B")
#' Nh_strata<-c("A"=20000, "AA"=20000, "B"=30000, "C"=20000, "D"=10000)
#' stehman2014(s, r, m, Nh_strata)
#' 
#' 
#' @export
stehman2014<-function(s, r, m, Nh_strata, margins=TRUE, order = sort(unique(r))){

  # check arguments
  s<-unlist(s)
  r<-unlist(r)
  m<-unlist(m)
  Nh_strata<-unlist(Nh_strata)
  if(is.null(names(Nh_strata))) stop("Nh_strata must be named.", call. = FALSE)
  if(length(names(Nh_strata))>length(unique(names(Nh_strata)))) stop("Repeated names detected in Nh_strata.", call. = FALSE)
  
  # compare arguments
  mapaccuracy:::.check_labels(names(Nh_strata), s)
  Nh_strata<-Nh_strata[names(Nh_strata)%in%s] # keep only the strata found in s
  mapaccuracy:::.check_labels(r, m)
  mapaccuracy:::.check_length(s, r, m)

  # covert arguments
  s<-factor(s, levels = names(Nh_strata))
  r<-factor(r, levels = order)
  m<-factor(m, levels = order)
  
  
  # equation 2 of Stehman 2014
  eq2<-function(col_name){
    sum(sample[,col_name]/sample$incl)/sum(Nh_strata)
  }
  
  # equation 25 of Stehman 2014
  eq25<-function(Nh, nh, var){
    sum(Nh^2*(1-nh/Nh)*var/nh)*(1/sum(Nh_strata)^2)
  }
  
  # equation 28 of Stehman 2014
  eq28<-function(R, X, Nh, nh, vary, varx, cov){
    sum(Nh^2*(1-nh/Nh)*(vary+R^2*varx-2*R*cov)/nh)*(1/X^2)
  }
  
  
  ## INCLUSION PROBABILITIES
  nh<-table(s)
  if(min(nh)<2){
    warning(sprintf("The following strata include only one observation: %s", paste(names(nh[nh<2]), collapse = "; ")),
            call. = FALSE)
  }
  nh<-nh[sapply(paste0("^",names(nh),"$"), grep, names(Nh_strata))] # match order of Nh_strata
  incl<-nh/Nh_strata
  incl<-sapply(s, function(x){incl[x]})
  

  ## CASE 1: when the objective is to estimate proportion of area of reference class k, 
  ## overall accuracy, and proportion of area in cell (i, j) of the error matrix

  sample<-data.frame(s=s, m=m, r=r, incl=incl)
  classes1<-levels(s)
  classes2<-levels(m)
  classes3<-levels(r)
  
  
  # yu for overall accuracy (eq 12 of Stehman 2014) 
  sample$yu_O<-as.numeric(m==r)
  
  # yu for proportion of area of reference classes (eq 14 of Stehman 2014)
  for (i in 1:length(classes3)){
    col_name<-paste0("yu_r",i)
    sample[, col_name]<-as.numeric(sample$r==classes3[i])
  }
  
  # yu for proportion of area with map class i and reference class j (eq 16 of Stehman 2014) 
  for (i in 1:length(classes2)){
    for (j in 1:length(classes3)){
      col_name<-paste0("yu_m",i,"r",j)
      sample[, col_name]<-as.numeric(sample$m==classes2[i] & sample$r==classes3[j])
    }
  }
  
  # Overall accuracy
  O<-eq2("yu_O")
  
  # proportion of area of reference classes
  A<-NULL
  for (i in 1:length(classes3)){
    col_name<-paste0("yu_r",i)
    A<-c(A, eq2(col_name))
  }
  names(A)<-order
  
  # proportion of area with map class i and reference class j
  M<-matrix(nrow = length(classes2), ncol = length(classes3))
  for (i in 1:length(classes2)){
    for (j in 1:length(classes3)){
      col_name<-paste0("yu_m",i,"r",j)
      M[i, j]<-ifelse(eq2(col_name)==0, NA, eq2(col_name))
    }
  }
  
  # Standard error (SE) for overall accuracy
  tmean<-stats::aggregate(.~sample$s, data=sample, mean)
  tvar<-stats::aggregate(.~sample$s, data=sample, stats::var)
  tvar[is.na(tvar)]<-0 # when sample$s includes strata with one point only, var is NA, which is replaced by zero
  
  Vo<-eq25(Nh=Nh_strata, nh=nh, var=tvar$yu_O)
  SEo<-sqrt(Vo)
  
  # Standard error (SE) for proportion of area of reference classes
  Va<-NULL
  for(i in 1:length(classes3)){
    Va<-c(Va, eq25(Nh=Nh_strata, nh=nh, var=tvar[,paste0("yu_r",i)])) 
  }
  SEa<-sqrt(Va)
  names(SEa)<-classes3
  
  # Standard error (SE) for proportion of area with map class i and reference class j ###### isto não aparece no exemplo numérico do artigo!
  Vm<-matrix(data=0, nrow = length(classes2), ncol = length(classes3))
  for (i in 1:length(classes2)){
    for (j in 1:length(classes3)){
      col_name<-paste0("yu_m",i,"r",j)
      Vm[i, j]<-Vm[i, j]+eq25(Nh=Nh_strata, nh=nh, var=tvar[,col_name])
    }
  }
  SEm<-sqrt(Vm)
  
  
  
  ## CASE 2: when the objective is to estimating a ratio R in the case 
  ## of the user’s accuracy, and producer’s accuracy,
  
  # yu for user's accuracy
  for (i in 1:length(classes2)){
    col_name<-paste0("yu_U",i)
    sample[, col_name]<-as.numeric(sample$m==classes2[i] & sample$m==sample$r)
  }
  for (i in 1:length(classes2)){
    col_name<-paste0("xu_U",i)
    sample[, col_name]<-as.numeric(sample$m==classes2[i])
  }
  
  # yu for producer's accuracy
  for (i in 1:length(classes3)){
    col_name<-paste0("yu_P",i)
    sample[, col_name]<-as.numeric(sample$r==classes3[i] & sample$r==sample$m)
  }
  for (i in 1:length(classes3)){
    col_name<-paste0("xu_P",i)
    sample[, col_name]<-as.numeric(sample$r==classes3[i])
  }
  
  # user's accuracy
  tmean<-stats::aggregate(.~sample$s, data=sample, mean)      # repeat this
  tvar<-stats::aggregate(.~sample$s, data=sample, stats::var) # repeat this
  tvar[is.na(tvar)]<-0 # when sample$s includes strata with one point only, var is NA, which is replaced by zero
  
  R_U_numerator<-NULL
  R_U_dominator<-NULL
  for (i in 1:length(classes2)){
    col_name<-paste0("yu_U",i)
    R_U_numerator<-c(R_U_numerator, sum(tmean[,col_name]*Nh_strata))
    col_name<-paste0("xu_U",i)
    R_U_dominator<-c(R_U_dominator, sum(tmean[,col_name]*Nh_strata))
  }
  U<-R_U_numerator/R_U_dominator
  names(U)<-classes2
  
  # producer's accuracy
  R_P_numerator<-NULL
  R_P_dominator<-NULL
  for (i in 1:length(classes3)){
    col_name<-paste0("yu_P",i)
    R_P_numerator<-c(R_P_numerator, sum(tmean[,col_name]*Nh_strata))
    col_name<-paste0("xu_P",i)
    R_P_dominator<-c(R_P_dominator, sum(tmean[,col_name]*Nh_strata))
  }
  P<-R_P_numerator/R_P_dominator
  names(P)<-classes2

  
  # Standard error (SE) for user's accuracy
  strata<-list()
  covP<-list()
  covU<-list()
  for(i in 1:length(classes1)){
    # a<-which(sample$m==classes1[i])
    strata[[i]]<-sample[sample$s==classes1[i],]
    covP[[i]]<-list()
    covU[[i]]<-list()
    for(j in 1:length(classes3)){
      covU[[i]][[j]]<-stats::cov(strata[[i]][,paste0("yu_U",j)],strata[[i]][,paste0("xu_U",j)])
      covP[[i]][[j]]<-stats::cov(strata[[i]][,paste0("yu_P",j)],strata[[i]][,paste0("xu_P",j)])
    }
  }
  
  # (simplify a bit the structure of the lists covU and covP)
  covU<-lapply(covU, unlist)
  covP<-lapply(covP, unlist)
  temp1<-as.data.frame(covU)
  temp2<-as.data.frame(covP)
  temp1<-t(temp1)
  temp2<-t(temp2)
  rownames(temp1)<-NULL
  rownames(temp2)<-NULL
  temp1[is.na(temp1)]<-0 # when sample$s includes strata with one point only
  temp2[is.na(temp2)]<-0 # when sample$s includes strata with one point only
  
  Vu<-NULL
  for(i in 1:length(classes2)){
    Vu<-c(Vu, eq28(R=U[i], X=R_U_dominator[i], Nh=Nh_strata, nh=nh[i], 
                   vary=tvar[,paste0("yu_U",i)], varx=tvar[,paste0("xu_U",i)], cov=temp1[,i]))
  }
  SEu<-sqrt(Vu)
  names(SEu)<-classes2
    
  # Standard error (SE) for producer's accuracy
  Vp<-NULL
  for(i in 1:length(classes3)){
    Vp<-c(Vp, eq28(R=P[i], X=R_P_dominator[i], Nh=Nh_strata, nh=nh[i], 
                   vary=tvar[,paste0("yu_P",i)], varx=tvar[,paste0("xu_P",i)], cov=temp2[,i]))
  }
  SEp<-sqrt(Vp)
  names(SEp)<-classes3
  
  
  # Add margins to confusion matrix
  colnames(M)<-classes3; rownames(M)<-classes2
  if(margins){
    M<-cbind(M,apply(M,1,sum, na.rm=TRUE))
    M<-rbind(M,c(apply(M,2,sum, na.rm=TRUE)))
    colnames(M)[length(colnames(M))]<-"sum"; rownames(M)[length(rownames(M))]<-"sum"
  }


  # return
  list(OA=O,
       UA=U,
       PA=P,
       area=A,
       SEoa=SEo,
       SEua=SEu,
       SEpa=SEp,
       SEa=SEa,
       matrix=M)
}
