doc_date <- function() 
{
  format(Sys.Date(), "%d %B %Y")
}

pkg_ver <- function(name) 
{
  paste(name, packageVersion(name))
}