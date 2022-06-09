require(broom)

convertCoxPHToTextOut = function(resCoxOS, FILENAME, OUTPUTPATH = "~/Dropbox/PhD/public/papers/CNSigs2 v6/_Excels_with_figure_data/") {
  
  # Convert to tidy
  x = tidy(resCoxOS)
  
  # Add estimate in exponential
  x$expEst = exp(x$estimate)
  
  # Add lower and upper boundaries
  x$lowerEst = exp(x$estimate-2*x$std.error)
  x$upperEst = exp(x$estimate+2*x$std.error)
  
  write.table(x, file.path(OUTPUTPATH, FILENAME), 
              sep='\t',row.names = FALSE, col.names = TRUE, quote=FALSE)
  
  return(x)
}