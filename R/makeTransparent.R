#' Function  to make transparent versions of R colors.

#' @param someColor takes a vector of color as \code{col} in \code{col2rgb}
#' @param alpha determines the transparency on a 0-255 scale
#'
#' @export

makeTransparent <- function(someColor, alpha = 100) {
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){
    rgb(red=curcoldata[1], green=curcoldata[2],
        blue=curcoldata[3], alpha=alpha, maxColorValue=255)
  })
}
