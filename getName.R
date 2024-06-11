get_name <- function(str1, number)
{
  res <- c()
  for (i in (1 : number)){
    this_str <- sprintf("%s%s", str1, as.character(i))
    res = c(res, this_str)
  }
  return(res)
}