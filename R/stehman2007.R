#' @title Thematic map accuracy for fuzzy classification/reference data under stratified random sampling.
#'
#' @description Implements the estimators described in Stehman \emph{et al}. (2007) for overall accuracy, producer's accuracy, and user's accuracy. Includes precision estimates.
#'
#' @inheritParams binaghi
#' @param Nh numeric vector. Mapped area (or pixel counts) of the classes. It must be named (see details).
#' @param old use the old code?
#'
#' @details The columns of both \code{r} and \code{m} must be named to explicitly and clearly identify the class that each column refers to.
#' Columns may be sorted differently.
#' Argument \code{Nh} must also be named to identify the class that each area refers to.
#' The order of \code{Nh} will be used for displaying the results.
#' Strata area (or pixel counts) is based on the map's spatial unit (e.g. pixel) maximum membership value.
#'
#' @return A list with the estimates and error matrix.
#' \item{OA}{overall accuracy}
#' \item{UA}{user's accuracy}
#' \item{PA}{producer's accuracy}
#' \item{SEoa}{standard error of OA}
#' \item{SEua}{standard error of UA}
#' \item{SEpa}{standard error of PA}
#' \item{matrix}{confusion error (proportion). Rows and columns represent map and reference class labels, respectively}
#'
#' @references
#' Stehman, S. V.; Arora, M. K.; Kasetkasem, T.; Varshney, P. K (2007). \href{https://doi.org/10.14358/PERS.73.2.165}{Estimation of fuzzy error matrix accuracy measures under stratified random sampling}. \emph{Photogramm. Eng. Remote Sensing}, 73, 165–173.
#'
#' @author Hugo Costa
#'
#' @examples
#' # Table A1, page 170 of Stehman et al. (2007)
#' r<-c(0.5, 0.2, 0.3, # stratum 1, pixel 1
#'         0.4, 0.2, 0.4, # stratum 1, pixel 2
#'         0.5, 0.2, 0.3, # stratum 1, pixel 3
#'         0.8, 0.2, 0.0, # stratum 1, pixel 4
#'         0.4, 0.6, 0.0, # stratum 2, pixel 1
#'         0.2, 0.6, 0.2, # stratum 2, pixel 2
#'         0.6, 0.1, 0.3, # stratum 2, pixel 3
#'         0.5, 0.5, 0.0, # stratum 2, pixel 4
#'         0.2, 0.1, 0.7, # stratum 3, pixel 1
#'         0.6, 0.4, 0.0, # stratum 3, pixel 2
#'         0.0, 0.3, 0.7, # stratum 3, pixel 3
#'         0.9, 0.0, 0.1) # stratum 3, pixel 4
#' r<-matrix(r,ncol=3,byrow = TRUE)
#'
#' m<-c(0.6, 0.3, 0.1, # stratum 1, pixel 1
#'         0.4, 0.3, 0.3, # stratum 1, pixel 2
#'         0.6, 0.0, 0.4, # stratum 1, pixel 3
#'         1.0, 0.0, 0.0, # stratum 1, pixel 4
#'         0.1, 0.7, 0.2, # stratum 2, pixel 1
#'         0.3, 0.4, 0.3, # stratum 2, pixel 2
#'         0.3, 0.5, 0.2, # stratum 2, pixel 3
#'         0.0, 1.0, 0.0, # stratum 2, pixel 4
#'         0.1, 0.1, 0.8, # stratum 3, pixel 1
#'         0.3, 0.3, 0.4, # stratum 3, pixel 2
#'         0.0, 0.3, 0.7, # stratum 3, pixel 3
#'         0.2, 0.3, 0.5) # stratum 3, pixel 4
#' m<-matrix(m,ncol = 3, byrow = TRUE)
#' colnames(r) <- colnames(m) <- c("N1", "N2", "N3")
#' Nh<-c("N1"=500, "N2"=300, "N3"=200)
#' e<-stehman2007(r, m, Nh)
#'
#' # compare to paper (differences due to number rounding):
#' e$OA          # eq. 170
#' e$SEoa^2      # eq. 170
#' e$SEoa        # eq. 171
#' e$UA[2]       # eq. 171
#' e$SEua[2]^2   # eq. 172 ######################## ERRO no meu codigo: should be 0.0826
#' e$SEua[2]     # eq. 172 ######################## ERRO no meu codigo: should be 0.287
#' e$PA[3]       # eq. 172
#' e$SEpa[3]^2   # eq. 172
#' e$SEpa[3]     # eq. 173
#' @export
stehman2007<-function(r, m, Nh, old=TRUE){
  # check arguments
  r<-as.matrix(r)
  m<-as.matrix(m)
  Nh<-unlist(Nh)
  if(is.null(colnames(r))) stop("The columns of r must be named.", call. = FALSE)
  if(is.null(colnames(m))) stop("The columns of m must be named.", call. = FALSE)
  if(is.null(names(Nh))) stop("Nh must be named.")
  if(length(colnames(r))>length(unique(colnames(r)))) stop("Repeated names detected in r.", call. = FALSE)
  if(length(colnames(m))>length(unique(colnames(m)))) stop("Repeated names detected in m.", call. = FALSE)
  if(length(names(Nh))>length(unique(names(Nh)))) stop("Repeated names detected in Nh.", call. = FALSE)
  if(nrow(r)!=nrow(m)){stop("Number of rows of 'r' and 'm' differ.", call. = FALSE)}

  # sort m and add missing columns if any
  virt<-matrix(data=1:ncol(r), ncol=ncol(r))
  colnames(virt)<-colnames(r)
  m<-gtools::smartbind(virt, m, fill=0)
  m<-m[-1,]
  m<-as.matrix(m)

  # stats
  q<-length(Nh)
  N<-sum(Nh)

  if(old){ ########################################################################################################## old code

    nh<-c(4,4,4)

    ###Divide data by stratum
    S<-list()
    for(s in 1:q){
      S[[s]]<-list()
      S[[s]][[1]]<-m[which(m[,s]==apply(m,1,max)),]
      S[[s]][[2]]<-r[which(m[,s]==apply(m,1,max)),]
      S[[s]][[3]]<-matrix((0*q*q),nrow=q,ncol=q) # matrices of zeros
    }

    ###Fuzzy mean error matrices for each stratum ("min" fuzzy set operator used)
    gu<-list()
    for (s in 1:q){
      gu[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        gu[[s]]<-0
        next
      }
      for (k in 1:nrow(S[[s]][[1]])){
        gu[[s]][[k]]<-0
        for(i in 1:q){
          for(j in 1:q){
            Eij<-min(S[[s]][[1]][k,i],S[[s]][[2]][k,j])
            S[[s]][[3]][i,j]<-S[[s]][[3]][i,j]+Eij/nrow(S[[s]][[1]])
            if(i==j){
              gu[[s]][[k]]<-c(gu[[s]][[k]],Eij)
            }
          }
        }
      }
    }

    ### Overall accuracy and variance
    diag<-0
    for(s in 1:q){
      diag<-diag+(sum(diag(S[[s]][[3]]))*Nh[s])
    }
    OA_Stehman<-diag/N

    for (s in 1:q){
      gu[[s]]<-sapply(gu[[s]],sum)
    }
    gu2<-sapply(gu,mean)
    gu3<-sapply(gu,stats::sd)^2
    gu3<-stats::na.omit(gu3)
    VAR_OA_Stehman<-(1/N^2)*sum(Nh[which(Nh>0)]^2/nh[which(nh>0)]*(1-nh[which(nh>0)]/Nh[which(Nh>0)])*gu3)
    SE_OA_Stehman<-sqrt(VAR_OA_Stehman)


    ### User's accuracy and Producer's accuracy
    M<-list()
    for (s in 1:q){
      M[[s]]<-apply(S[[s]][[1]],2,mean)
    }
    R<-list()
    for (s in 1:q){
      R[[s]]<-apply(S[[s]][[2]],2,mean)
    }

    UA_Stehman<-NULL
    PA_Stehman<-NULL
    numerator<-list()
    denominatorM<-list()
    denominatorR<-list()
    for (i in 1:q){
      numerator[[i]]<-0
      denominatorM[[i]]<-0
      denominatorR[[i]]<-0
      for(s in 1:q){
        numerator[[i]]<- sum(numerator[[i]], S[[s]][[3]][i,i]*Nh[s], na.rm =T)
        denominatorM[[i]]<- sum(denominatorM[[i]], M[[s]][i]*Nh[s], na.rm =T)
        denominatorR[[i]]<- sum(denominatorR[[i]], R[[s]][i]*Nh[s], na.rm =T)
      }
      UA_Stehman<-c(UA_Stehman,numerator[[i]]/denominatorM[[i]])
      PA_Stehman<-c(PA_Stehman,numerator[[i]]/denominatorR[[i]])
    }
    UA_Stehman[is.na(UA_Stehman)] <- 0
    PA_Stehman[is.na(PA_Stehman)] <- 0
    names(UA_Stehman)<-colnames(r)
    names(PA_Stehman)<-colnames(r)

    ### variance UA_Stehman and PA_Stehman
    ev<-list()
    for (s in 1:q){
      ev[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        ev[[s]]<-0
        next
      }
      for(i in 1:q){
        temp<-NULL
        for (k in 1:nrow(S[[s]][[1]])){
          Eij<-min(S[[s]][[1]][k,i],S[[s]][[2]][k,i])
          temp<-c(temp,Eij)
          ev[[s]][[i]]<-temp
        }
      }
    }
    mv<-list()
    rv<-list()
    for (s in 1:q){
      mv[[s]]<-list()
      rv[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        mv[[s]]<-0
        rv[[s]]<-0
        next
      }
      for(i in 1:q){
        temp<-NULL
        tempr<-NULL
        for (k in 1:nrow(S[[s]][[1]])){
          Ei<-S[[s]][[1]][k,i]
          Ri<-S[[s]][[2]][k,i]
          temp<-c(temp,Ei)
          tempr<-c(tempr,Ri)
          mv[[s]][[i]]<-temp
          rv[[s]][[i]]<-tempr
        }
      }
    }
    emv<-list()
    erv<-list()
    for (s in 1:q){
      emv[[s]]<-list()
      erv[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        emv[[s]]<-0
        erv[[s]]<-0
        next
      }
      for(i in 1:q){
        emv[[s]][[i]]<-ev[[s]][[i]]*mv[[s]][[i]]
        erv[[s]][[i]]<-ev[[s]][[i]]*rv[[s]][[i]]
      }
    }
    ev2<-list()
    ev3<-list()
    mv2<-list()
    mv3<-list()
    emv2<-list()
    rv2<-list()
    rv3<-list()
    erv2<-list()
    for (s in 1:q){
      ev2[[s]]<-sapply(ev[[s]],mean)
      ev3[[s]]<-sapply(ev[[s]],stats::sd)^2
      mv2[[s]]<-sapply(mv[[s]],mean)
      mv3[[s]]<-sapply(mv[[s]],stats::sd)^2
      emv2[[s]]<-sapply(emv[[s]],sum)
      rv2[[s]]<-sapply(rv[[s]],mean)
      rv3[[s]]<-sapply(rv[[s]],stats::sd)^2
      erv2[[s]]<-sapply(erv[[s]],sum)
    }
    #just change the way calculations are stored
    e1<-list()
    e2<-list()
    m1<-list()
    m2<-list()
    em<-list()
    er1<-list()
    er2<-list()
    r1<-list()
    r2<-list()
    er<-list()
    for (i in 1:q){
      tempe1<-NULL
      tempe2<-NULL
      tempm1<-NULL
      tempm2<-NULL
      tempem<-NULL
      temper1<-NULL
      temper2<-NULL
      tempr1<-NULL
      tempr2<-NULL
      temper<-NULL
      for(s in 1:q){
        if(length(ev2[[s]])==1){ #evitar erro no c?digo
          tempe1<-c(tempe1,0)
          tempe2<-c(tempe2,0)
          tempm1<-c(tempm1,0)
          tempm2<-c(tempm2,0)
          tempem<-c(tempem,0)
          temper1<-c(temper1,0)
          temper2<-c(temper2,0)
          tempr1<-c(tempr1,0)
          tempr2<-c(tempr2,0)
          temper<-c(temper,0)
          next
        }else{
          tempe1<-c(tempe1,ev2[[s]][[i]])
          tempe2<-c(tempe2,ev3[[s]][[i]])
          tempm1<-c(tempm1,mv2[[s]][[i]])
          tempm2<-c(tempm2,mv3[[s]][[i]])
          tempem<-c(tempem,emv2[[s]][[i]])
          temper1<-c(temper1,ev2[[s]][[i]])
          temper2<-c(temper2,ev3[[s]][[i]])
          tempr1<-c(tempr1,rv2[[s]][[i]])
          tempr2<-c(tempr2,rv3[[s]][[i]])
          temper<-c(temper,erv2[[s]][[i]])
        }
      }
      e1[[i]]<-tempe1
      e2[[i]]<-tempe2
      m1[[i]]<-tempm1
      m2[[i]]<-tempm2
      em[[i]]<-tempem
      er1[[i]]<-temper1
      er2[[i]]<-temper2
      r1[[i]]<-tempr1
      r2[[i]]<-tempr2
      er[[i]]<-temper
    }
    ss<-list()
    ssr<-list()
    for (s in 1:q){
      temp<-NULL
      tempr<-NULL
      for (i in 1:q){
        temp<-c(temp,(em[[s]][[i]]-nh[s]*e1[[s]][[i]]*m1[[s]][[i]])/(nh[s]-1))
        tempr<-c(tempr,(er[[s]][[i]]-nh[s]*er1[[s]][[i]]*r1[[s]][[i]])/(nh[s]-1))
      }
      ss[[s]]<-temp
      ssr[[s]]<-tempr
    }

    VAR_UA_Stehman<-NULL
    VAR_PA_Stehman<-NULL
    for (i in 1:q){
      temp<-(e2[[i]]+UA_Stehman[i]^2*m2[[i]]-2*UA_Stehman[i]*ss[[i]])
      VAR_UA_Stehman<-c(VAR_UA_Stehman,(1/unlist(denominatorM)[i]^2)*sum(Nh^2*(1-nh/Nh)/nh*temp,na.rm=T))
      tempr<-(er2[[i]]+PA_Stehman[i]^2*r2[[i]]-2*PA_Stehman[i]*ssr[[i]])
      VAR_PA_Stehman<-c(VAR_PA_Stehman,(1/unlist(denominatorR)[i]^2)*sum(Nh^2*(1-nh/Nh)/nh*tempr,na.rm=T))
    }
    VAR_UA_Stehman[is.na(VAR_UA_Stehman)]<-0
    VAR_PA_Stehman[is.na(VAR_PA_Stehman)]<-0
    names(VAR_UA_Stehman)<-colnames(r)
    names(VAR_PA_Stehman)<-colnames(r)
    SE_UA_Stehman<-sqrt(VAR_UA_Stehman)
    SE_PA_Stehman<-sqrt(VAR_PA_Stehman)



  }else{
    ########################################################################################################### New code
    ########################################################################################################### New code
    ########################################################################################################### New code
    ########################################################################################################### New code
    ########################################################################################################### New code
    ########################################################################################################### New code
    ########################################################################################################### New code

    ###Divide data by stratum (tables A2 and A3)
    S<-list()
    S2<-list()           #A2
    S3<-list()           #A2
    M2<-list()
    R2<-list()
    nh<-NULL
    for(s in 1:q){
      S[[s]]<-list()
      S2[[s]]<-list()

      # isolate data by stratum
      temp1<-r[which(m[,s]==apply(m,1,max)),]
      temp2<-m[which(m[,s]==apply(m,1,max)),]
      nh<-c(nh, nrow(temp2))
      S[[s]][[1]]<-temp2
      S[[s]][[2]]<-temp1
      S3[[s]]<-binaghi(temp1,temp2, margins=FALSE)$matrix/nrow(temp2) # Table A3
      R2[[s]]<-apply(temp1,2,mean)
      M2[[s]]<-apply(temp2,2,mean)

      for(i in 1:nrow(temp1)){
        row1<-matrix(temp1[i,],ncol=q)
        row2<-matrix(temp2[i,],ncol=q)
        colnames(row1)<-colnames(row2)<-colnames(r)
        S2[[s]][[i]]<-binaghi(row1,row2, margins=FALSE)$matrix
      }
    }

    # table A4
    S4<-list()
    for(i in 1:length(S2)){
      S4[[i]]<-list()
      S2[[i]]<-lapply(S2[[i]],as.matrix)
      S2[[i]]<-lapply(S2[[i]],diag)
      sapply(S2[[i]],sum)
      S4[[i]][[1]]<-mean(sapply(S2[[i]],sum)) # table A4
      S4[[i]][[2]]<-stats::sd(sapply(S2[[i]],sum))^2 # table A4
    }
    S3<-lapply(S3,as.matrix)
    ( S3diag<-lapply(S3,diag) )
    S3diag<-sapply(S3diag,sum)*Nh
    OA_Stehman2<-sum(S3diag)/sum(Nh)

    #( gu22<-unlist(sapply(S4,'[',1)) )
    ( gu32<-unlist(sapply(S4,'[',2)) )

    VAR_OA_Stehman<-(1/sum(Nh)^2)*sum(Nh[which(Nh>0)]^2/nh[which(nh>0)]*(1-nh[which(nh>0)]/Nh[which(Nh>0)])*gu32)
    SE_OA_Stehman<-sqrt(VAR_OA_Stehman)


    ### User's accuracy and Producer's accuracy
    UA_Stehman2<-NULL
    PA_Stehman2<-NULL
    numerator<-list()
    denominatorM<-list()
    denominatorR<-list()
    for (i in 1:q){
      numerator[[i]]<-0
      denominatorM[[i]]<-0
      denominatorR[[i]]<-0
      for(s in 1:q){
        numerator[[i]]<- sum(numerator[[i]], S3[[s]][i,i]*Nh[s], na.rm =T)
        denominatorM[[i]]<- sum(denominatorM[[i]], M2[[s]][i]*Nh[s], na.rm =T)
        denominatorR[[i]]<- sum(denominatorR[[i]], R2[[s]][i]*Nh[s], na.rm =T)
      }
      UA_Stehman2<-c(UA_Stehman2,numerator[[i]]/denominatorM[[i]])
      PA_Stehman2<-c(PA_Stehman2,numerator[[i]]/denominatorR[[i]])
    }
    UA_Stehman2[is.na(UA_Stehman2)] <- 0
    PA_Stehman2[is.na(PA_Stehman2)] <- 0
    names(UA_Stehman2)<-colnames(r)
    names(PA_Stehman2)<-colnames(r)

    ### variance UA_Stehman2 and PA_Stehman2
    ev<-list()
    for (s in 1:q){
      ev[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        ev[[s]]<-0
        next
      }
      for(i in 1:q){
        temp<-NULL
        for (k in 1:nrow(S[[s]][[1]])){
          Eij<-min(S[[s]][[1]][k,i],S[[s]][[2]][k,i])
          temp<-c(temp,Eij)
          ev[[s]][[i]]<-temp
        }
      }
    }
    mv<-list()
    rv<-list()
    for (s in 1:q){
      mv[[s]]<-list()
      rv[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        mv[[s]]<-0
        rv[[s]]<-0
        next
      }
      for(i in 1:q){
        temp<-NULL
        tempr<-NULL
        for (k in 1:nrow(S[[s]][[1]])){
          Ei<-S[[s]][[1]][k,i]
          Ri<-S[[s]][[2]][k,i]
          temp<-c(temp,Ei)
          tempr<-c(tempr,Ri)
          mv[[s]][[i]]<-temp
          rv[[s]][[i]]<-tempr
        }
      }
    }
    emv<-list()
    erv<-list()
    for (s in 1:q){
      emv[[s]]<-list()
      erv[[s]]<-list()
      if(nrow(S[[s]][[1]])==0){
        emv[[s]]<-0
        erv[[s]]<-0
        next
      }
      for(i in 1:q){
        emv[[s]][[i]]<-ev[[s]][[i]]*mv[[s]][[i]]
        erv[[s]][[i]]<-ev[[s]][[i]]*rv[[s]][[i]]
      }
    }
    ev2<-list()
    ev3<-list()
    mv2<-list()
    mv3<-list()
    emv2<-list()
    rv2<-list()
    rv3<-list()
    erv2<-list()
    for (s in 1:q){
      ev2[[s]]<-sapply(ev[[s]],mean)
      ev3[[s]]<-sapply(ev[[s]],stats::sd)^2
      mv2[[s]]<-sapply(mv[[s]],mean)
      mv3[[s]]<-sapply(mv[[s]],stats::sd)^2
      emv2[[s]]<-sapply(emv[[s]],sum)
      rv2[[s]]<-sapply(rv[[s]],mean)
      rv3[[s]]<-sapply(rv[[s]],stats::sd)^2
      erv2[[s]]<-sapply(erv[[s]],sum)
    }
    #just change the way calculations are stored
    e1<-list()
    e2<-list()
    m1<-list()
    m2<-list()
    em<-list()
    er1<-list()
    er2<-list()
    r1<-list()
    r2<-list()
    er<-list()
    for (i in 1:q){
      tempe1<-NULL
      tempe2<-NULL
      tempm1<-NULL
      tempm2<-NULL
      tempem<-NULL
      temper1<-NULL
      temper2<-NULL
      tempr1<-NULL
      tempr2<-NULL
      temper<-NULL
      for(s in 1:q){
        if(length(ev2[[s]])==1){ #evitar erro no c?digo
          tempe1<-c(tempe1,0)
          tempe2<-c(tempe2,0)
          tempm1<-c(tempm1,0)
          tempm2<-c(tempm2,0)
          tempem<-c(tempem,0)
          temper1<-c(temper1,0)
          temper2<-c(temper2,0)
          tempr1<-c(tempr1,0)
          tempr2<-c(tempr2,0)
          temper<-c(temper,0)
          next
        }else{
          tempe1<-c(tempe1,ev2[[s]][[i]])
          tempe2<-c(tempe2,ev3[[s]][[i]])
          tempm1<-c(tempm1,mv2[[s]][[i]])
          tempm2<-c(tempm2,mv3[[s]][[i]])
          tempem<-c(tempem,emv2[[s]][[i]])
          temper1<-c(temper1,ev2[[s]][[i]])
          temper2<-c(temper2,ev3[[s]][[i]])
          tempr1<-c(tempr1,rv2[[s]][[i]])
          tempr2<-c(tempr2,rv3[[s]][[i]])
          temper<-c(temper,erv2[[s]][[i]])
        }
      }
      e1[[i]]<-tempe1
      e2[[i]]<-tempe2
      m1[[i]]<-tempm1
      m2[[i]]<-tempm2
      em[[i]]<-tempem
      er1[[i]]<-temper1
      er2[[i]]<-temper2
      r1[[i]]<-tempr1
      r2[[i]]<-tempr2
      er[[i]]<-temper
    }
    ss<-list()
    ssr<-list()
    for (s in 1:q){
      temp<-NULL
      tempr<-NULL
      for (i in 1:q){
        temp<-c(temp,(em[[s]][[i]]-nh[s]*e1[[s]][[i]]*m1[[s]][[i]])/(nh[s]-1))
        tempr<-c(tempr,(er[[s]][[i]]-nh[s]*er1[[s]][[i]]*r1[[s]][[i]])/(nh[s]-1))
      }
      ss[[s]]<-temp
      ssr[[s]]<-tempr
    }

    VAR_UA_Stehman<-NULL
    VAR_PA_Stehman<-NULL
    for (i in 1:q){
      temp<-(e2[[i]]+UA_Stehman2[i]^2*m2[[i]]-2*UA_Stehman2[i]*ss[[i]])
      VAR_UA_Stehman<-c(VAR_UA_Stehman,(1/unlist(denominatorM)[i]^2)*sum(Nh^2*(1-nh/Nh)/nh*temp,na.rm=T))
      tempr<-(er2[[i]]+PA_Stehman2[i]^2*r2[[i]]-2*PA_Stehman2[i]*ssr[[i]])
      VAR_PA_Stehman<-c(VAR_PA_Stehman,(1/unlist(denominatorR)[i]^2)*sum(Nh^2*(1-nh/Nh)/nh*tempr,na.rm=T))
    }
    VAR_UA_Stehman[is.na(VAR_UA_Stehman)]<-0
    VAR_PA_Stehman[is.na(VAR_PA_Stehman)]<-0
    names(VAR_UA_Stehman)<-colnames(r)
    names(VAR_PA_Stehman)<-colnames(r)
    SE_UA_Stehman<-sqrt(VAR_UA_Stehman)
    SE_PA_Stehman<-sqrt(VAR_PA_Stehman)

    # fuzzy matrix: EXPERIMENTAL
    matrix<-list()
    for (s in 1:q){
      matrix[[s]]<-S3[[s]]*Nh[s]/N
    }
    matrix<-Reduce('+',matrix)

    OA_Stehman<-OA_Stehman2
    UA_Stehman<-UA_Stehman2
    PA_Stehman<-PA_Stehman2



  } #end else old/new



  # return
  list(OA=OA_Stehman,
       UA=UA_Stehman,
       PA=PA_Stehman,
       # area=tarea,
       SEoa=SE_OA_Stehman,
       SEua=SE_UA_Stehman,
       SEpa=SE_PA_Stehman)#,
       #SEa =SEa,
       #matrix=matrix)

}
