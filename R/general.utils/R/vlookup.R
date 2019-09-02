#' vlookup
#'
#' Vlookup functionality from excel
#' @param this df key value
#' @keywords 
#' @export
#' @examples
#' vlookup(this, df, key, value)

require(tidyverse)

vlookup <- function(this, df, key, value) {
  m <- match(this, df[[key]])
  df[[value]][m]
}