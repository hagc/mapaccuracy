#' Internal functions.
#'
#' Internal functions to be called from within another functions.
#'
#' \code{.check_classes} checks (tests) if the arguments passed to functions comply with their expected classes.
#'
#' \code{.check_factors} delete?.
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
#' @examples
#' mapaccuracy2:::.check_classes("character", "a")                 # no error found
#' mapaccuracy2:::.check_classes(c("character","list"), "a", "a")  # no error found
#' mapaccuracy2:::.check_classes(c("character","list"), list("a")) # no error found
#' \dontrun{
#'   mapaccuracy2:::.check_classes("character", list("a"))         # error (not a character)
#' }
#' fun<-function(x){
#'   mapaccuracy2:::.check_classes("numeric",x)
#'   print(x+1)
#' }
#' fun(1)       # no error found
#' \dontrun{
#'   fun("a")   # error as argument is not numeric
#' }
#'
#' mapaccuracy2:::.check_factors("character", "a")                 # no error found
#' mapaccuracy2:::.check_factors(c("character","list"), "a", "a")  # no error found
#' mapaccuracy2:::.check_factors(c("character","list"), list("a")) # no error found
#' \dontrun{
#'   mapaccuracy2:::.check_factors("character", list("a"))         # error (not a character)
#' }
#' fun<-function(x){
#'   mapaccuracy2:::.check_factors("numeric",x)
#'   print(x+1)
#' }
#' fun(1)       # no error found
#' \dontrun{
#'   fun("a")   # error as argument is not numeric
#' }
#'
#' a<-c("a","b","c")
#' b<-c("a","b","d")
#' mapaccuracy2:::.check_labels(a,a)
#' mapaccuracy2:::.check_labels(a,a[1:2])
#' mapaccuracy2:::.check_labels(a,a,a)
#' \dontrun{
#' mapaccuracy2:::.check_labels(a,b)
#' }
#'
#' mapaccuracy2:::.check_length("a", "b")
#' \dontrun{
#'   mapaccuracy2:::.check_length("a", c("b","b"))
#' }
#'
#' @name .mapaccuracy-internal
NULL

#' @rdname .mapaccuracy-internal
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

#' @rdname .mapaccuracy-internal
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

#' @rdname .mapaccuracy-internal
.check_labels<-function(x, ...){
  x<-unique(x)
  arguments=sapply(list(...),unique)
  if(!all(sapply(arguments,'%in%',x))){
    stop("Arguments should include only: ", paste(x,collapse = ", "), call. = FALSE)
  }
}

#' @rdname .mapaccuracy-internal
.check_length<-function(...){
  arguments<-list(...)
  if(length(arguments)>1){
    length<-lapply(arguments,length)
    if(!do.call(assertthat::are_equal,length)) stop("Length of arguments differ.", call. = FALSE)
  }
}
