#' Thematic map accuracy and area.
#'
#' Implements the estimators described in Olofsson \emph{et al}. (2013, 2014) for overall accuracy,
#' producer's accuracy, user's accuracy, and area. Includes precision estimates.
#'
#' Argument \code{Nh} must be named to explicitly and clearly identify the class that each area refers to.
#' The order of \code{Nh} will be used for displaying the results.
#'
#' The labels included in \code{r} and \code{m} must match among them and the names of \code{Nh}.
#' 
#' In the error matrix returned, the entries corresponding to no observed cases will present
#' \code{NA} rather than \code{0}. This is to emphasize the difference between the absence of cases
#' and the presence of some (few) cases that represent a very small proportion of area (almost
#' zero) and thus possibly rounded to zero. However, \code{NA} means zero proportion of area.
#'
#' @param r character vector. Reference class labels. The object will be coerced to factor.
#' @param m character vector. Map class labels. The object will be coerced to factor.
#' @param Nh numeric vector. Area, number of pixels, or proportion of the classes in the map. 
#'           It must be named (see details).
#' @param margins logical. If FALSE, the error matrix produced includes no margins (sum of the 
#'                rows and columns).
#' @param csv_filepath (optional) character. Write results to a csv file (not implemented)
#'
#' @return A list with the estimates and error matrix, and a csv if defined.
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
#' Olofsson, P.; Foody, G. M.; Stehman, S. V.; Woodcock, C. E. (2013). \href{https://doi.org/10.1016/j.rse.2012.10.031}{Making better use of accuracy data in land change studies: Estimating accuracy and area and quantifying uncertainty using stratified estimation}. \emph{Remote Sens. Environ.}, 129, 122-131.
#'
#' Olofsson, P.; Foody, G. M.; Herold, M.; Stehman, S. V.; Woodcock, C. E.; Wulder, M. A. (2014). \href{https://doi.org/10.1016/j.rse.2014.02.015}{Good practices for estimating area and assessing accuracy of land change}. \emph{Remote Sens. Environ.}, 148, 42-57.
#'
#' @author Hugo Costa
#'
#' @examples
#' ## Example 1 in Olofsson et al. (2013)
#' r<-c(rep("1",102),rep("2",280),rep("3",118))
#' m<-c(rep("1",97) ,rep("2",3), rep("3",2),rep("2",279),
#'               "3",rep("1",3),rep("2",18),rep("3",97))
#' Nh<-c(22353, 1122543, 610228)
#' names(Nh)<-c("1", "2", "3")
#' a<-olofsson(r, m, Nh)
#'
#' # compare to paper:
#' a$area[1]                                       # eq. 9
#' a$area[1]*sum(Nh)                               # eq. 10
#' a$SEa[1]*sum(Nh)                                # eq. 12
#' a$area[1]*sum(Nh)-qnorm(0.975)*a$SEa[1]*sum(Nh) # 95% CI lower bound (note typo in the paper)
#' a$area[1]*sum(Nh)+qnorm(0.975)*a$SEa[1]*sum(Nh) # 95% CI upper bound
#' a$UA[1]                                         # eq. 14
#' a$PA[1]                                         # eq. 15
#' a$OA                                            # eq. 16
#' a$UA                                            # table 4
#' qnorm(0.975)*a$SEua                             # table 4
#' a$PA                                            # table 4
#' qnorm(0.975)*a$SEpa                             # table 4
#' a$matrix                                        # table 4
#'
#'
#' ## Example 2 in Olofsson et al. (2013)
#' r<-c(rep("1", 129), rep("2", 403), rep("3", 611))
#' m<-c(rep("1", 127), "2", "2", rep("1", 66), rep("2", 322), rep("3", 15), rep("1", 54),
#'      rep("2", 17), rep("3", 540))
#' Nh<-c(0.007, 0.295, 0.698)
#' names(Nh)<-c("1", "2", "3")
#' b<-olofsson(r, m, Nh)
#' 
#' # compare to paper (table 6):
#' b$OA
#' qnorm(0.975)*b$SEoa
#' b$UA
#' qnorm(0.975)*b$SEua
#' b$PA
#' qnorm(0.975)*b$SEpa
#' 
#'
#' ## Example of table 8 in Olofsson et al. (2014)
#' r<-c(rep(1,69),rep(2,56),rep(3,175),rep(4,340))
#' m<-c(rep(1,66), 3, rep(4,2), rep(2,55), 4, rep(1,5), rep(2,8),
#'               rep(3,153),rep(4,9),rep(1,4),rep(2,12),rep(3,11),rep(4,313))
#' r[r==1] <- m[m==1] <- "Deforestation"
#' r[r==2] <- m[m==2] <- "Forest gain"
#' r[r==3] <- m[m==3] <- "Stable forest"
#' r[r==4] <- m[m==4] <- "Stable non-forest"
#' Nh<-c("Deforestation"=200000, "Forest gain"=150000,
#'       "Stable forest"=3200000, "Stable non-forest"=6450000) * 30^2 # Landsat pixel area = 30^2
#' e<-olofsson(r, m, Nh)
#'
#' # compare to paper, left-hand of p. 54:
#' e$UA                            # User's accuracy
#' qnorm(0.975)*e$SEua             # 95% CI width
#' e$PA                            # Producer's accuracy
#' qnorm(0.975)*e$SEpa             # 95% CI width
#' e$OA                            # Overall accuracy
#' qnorm(0.975)*e$SEoa             # 95% CI width
#'
#' # compare to paper, right-hand of p. 54:
#' e$area[1]*sum(Nh)/10000         # deforestation in hectares
#' qnorm(0.975)*e$SEa[1]*sum(Nh)/10000     # 95% CI width in hectares
#' e$area[2]*sum(Nh)/10000                 # forest gain in hectares
#' qnorm(0.975)*e$SEa[2]*sum(Nh)/10000     # 95% CI width in hectares
#' e$area[3]*sum(Nh)/10000                 # stable forest in hectares
#' qnorm(0.975)*e$SEa[3]*sum(Nh)/10000     # 95% CI width in hectares
#' e$area[4]*sum(Nh)/10000                 # stable non-forest in hectares
#' qnorm(0.975)*e$SEa[4]*sum(Nh)/10000     # 95% CI width in hectares
#'
#' # change class order
#' order<-c(4,2,1,3)
#' olofsson(r, m, Nh[c(4,2,1,3)])
#' @export
olofsson<-function(r, m, Nh, margins=TRUE, csv_filepath){

  # check arguments
  r<-unlist(r)
  m<-unlist(m)
  Nh<-unlist(Nh)
  if(is.null(names(Nh))) stop("Nh must be named.", call. = FALSE)
  if(length(names(Nh))>length(unique(names(Nh)))) stop("Repeated names detected in Nh.", call. = FALSE)
  .check_labels(names(Nh), r, m)
  .check_labels(r, names(Nh))
  .check_length(r,m)

  # convert arguments
  r<-factor(r, levels = names(Nh))
  m<-factor(m, levels = names(Nh))


  # error matrix
  matrix<-table(m, r)

  # Sampling settings
  q <-length(Nh)       # number of classes
  A <-sum(Nh)          # total area
  nh<-rowSums(matrix)  # stratum sample size
  Wh<-Nh/A             # map class proportions

  # confusion matrix (estimated area proportions)
  props<-prop.table(matrix,1)
  for(j in 1:q){
    props[j,]<-Wh[j]*props[j,]
  }

  # Accuracy and area estimates
  OA<-sum(diag(props),na.rm=T)
  UA<-diag(props)/rowSums(props)
  PA<-diag(props)/colSums(props)
  # tarea<-A*colSums(props)
  tarea<-colSums(props)


  # standard error of OA
  VAR_OA<-0
  for(j in 1:q){
    temp<-Wh[j]^2*UA[j]*(1-UA[j])/(nh[j]-1)
    if(is.na(temp)){temp<-0}
    if(is.infinite(temp)){temp<-0}
    VAR_OA<-VAR_OA + temp
  }
  SEoa<-sqrt(VAR_OA)
  names(SEoa)<-NULL


  # standard error of UA
  VAR_UA<-NULL
  for(j in 1:q){
    temp<-UA[j]*(1-UA[j])/(nh[j]-1)
    if(is.na(temp)){temp<-0}
    if(is.infinite(temp)){temp<-0}
    VAR_UA<-c(VAR_UA, temp)
  }
  SEua<-sqrt(VAR_UA)


  # standard error of PA
  N.j<-NULL
  for(j in 1:q){
    temp<-0
    for(i in 1:q){
      temp<-sum(temp,Nh[i]/nh[i]*matrix[i,j],na.rm=T)
    }
    N.j<-c(N.j, temp)
  }

  VAR_PA<-NULL
  for(j in 1:q){
    # temp1<-N.j[j]^2*(1-PA[j])^2*UA[j]*(1-UA[j])/(nh[j]-1)
    temp1<-Nh[j]^2*(1-PA[j])^2*UA[j]*(1-UA[j])/(nh[j]-1)
    if(is.na(temp1)){temp1<-0}
    if(is.infinite(temp1)){temp1<-0}
    temp2<-0
    seque<-1:q
    seque<-seque[-j]
    for(i in seque){
      temp2<-sum(temp2, Nh[i]^2*matrix[i,j]/nh[i]*(1-matrix[i,j]/nh[i])/(nh[i]-1), na.rm=T)
    }
    if(is.na(temp2)){temp2<-0}
    if(is.infinite(temp2)){temp2<-0}

    VAR_PA<-c(VAR_PA, (1/N.j[j]^2)*(temp1+PA[j]^2*temp2))
  }
  SEpa<-sqrt(VAR_PA)


  # standard error of the area
  VAR_A<-list()
  for(j in 1:q){
    a<-Wh^2
    b<-matrix[,j]/nh
    c<-(1-matrix[,j]/nh)
    d<-nh-1
    VAR_A[[j]]<-sum(a*(b*c/d))
  }
  VAR_A<-unlist(VAR_A)
  SEa<-sqrt(VAR_A)
  # SEa<-A*SEa
  names(SEa)<-names(Nh)

  # gather calculations together
  if(!missing(csv_filepath)){
    # export as a table to csv
  }

  # Add margins to confusion matrix (potentially useful for exporting only)
  props[matrix==0]<-NA
  if(margins){
    props<-cbind(props,rowSums(props, na.rm = TRUE))
    props<-rbind(props,c(colSums(props, na.rm = TRUE)))
    colnames(props)[ncol(props)]<-"sum"; rownames(props)[nrow(props)]<-"sum"
  }

  # return
  list(OA=OA,
       UA=UA,
       PA=PA,
       area=tarea,
       SEoa=SEoa,
       SEua=SEua,
       SEpa=SEpa,
       SEa =SEa,
       matrix=props)
}
