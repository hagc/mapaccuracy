#' Internal functions.
#'
#' Internal functions to be called from within another functions.
#'
#' \code{.check_classes} checks (tests) if the arguments passed to functions comply with their expected classes.
#'
#' \code{.check_labels} checks (tests) if the character arguments passed to functions comply with their expected content.
#'
#' \code{.check_length} checks (tests) if the length of objects match.
#'
#' @keywords internal
#'
#' @param x expected objects (character vector)
#' @param ...  arguments to be checked.
#'
#' @return Error message if the test fails. Nothing if the test passes.
#' 
#' @name internal.mapaccuracy
NULL

#' @rdname internal.mapaccuracy
.check_classes<-function(x, ...){

  # check argument classes
  classes.expected<-as.character(x)
  allclasses<-c("character","factor","numeric","data.frame","matrix","list","logical")
  if(!all(sapply(classes.expected,'%in%',allclasses))) stop("Argument 'classes' includes one or more classes not recognized. Use only the following: ", paste(allclasses,collapse = ", "), call. = FALSE)


  # check if the arguments provided in '...' match at least one of the classes defined in 'classes')
  arguments<-list(...)
  classes.objects<-sapply(arguments,class)
  if(!any(sapply(classes.objects,'%in%',classes.expected))) stop("One or more arguments passed to the function are not of the expected class.", call. = FALSE)

}

#' @rdname internal.mapaccuracy
.check_factors<-function(x, ...){

  # check argument factors
  factors.expected<-factor(x)
  allfactors<-c("character","factor","numeric","data.frame","matrix","list","logical")
  if(!all(sapply(factors.expected,'%in%',allfactors))) stop("Argument 'factors' includes one or more factors not recognized. Use only the following: ",paste(allfactors,collapse = ", "))


  # check if the arguments provided in '...' match at least one of the factors defined in 'factors')
  arguments<-list(...)
  factors.objects<-sapply(arguments,class)
  if(!any(sapply(factors.objects,'%in%',factors.expected))) stop("One or more arguments passed to the function are not of the expected class.", call. = FALSE)

}

#' @rdname internal.mapaccuracy
.check_labels<-function(x, ...){
  x<-unique(x)
  arguments=sapply(list(...),unique)
  if(!all(sapply(arguments,'%in%',x))){
    stop("Arguments should include only: ", paste(x,collapse = ", "), call. = FALSE)
  }
}

#' @rdname internal.mapaccuracy
.check_length<-function(...){
  arguments<-list(...)
  if(length(arguments)>1){
    length<-lapply(arguments,length)
    if(!do.call(assertthat::are_equal,length)) stop("Length of arguments differ.", call. = FALSE)
  }
}
