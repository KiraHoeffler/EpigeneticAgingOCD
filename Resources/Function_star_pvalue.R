
star_p_value <- function(pvalue){
  if (pvalue < 0.05){
    text_star <- "*"
  } 
  if (pvalue < 0.01) {
    text_star <- "**"
  }
  if (pvalue < 0.001) {
    text_star <- "***"
  } 
  if (pvalue < 0.0001) {
    text_star <- "****"
  } 
  if (pvalue >= 0.05) {
    text_star <- "ns"
  } 
  return(text_star)
}
