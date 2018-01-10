#' @title qc_outlierCheck
#' @description This function takes a dataframe, and does a simple check for 
#' outliers based on the boxplot.stats() routine. Specifically, it identifies 
#' data which is 1.5, 5, or 10 times the length of the box away from the box.  
#' These are identified as "Possible", "Likely", or "Extreme" in the resultant 
#' "OL_<field>" column, respectively. Other data is identified as "OK", which 
#' indicates that it is within a reasonable range of the mean.
#' @param df a dataframe to be analyzed.
#' @param field a fieldname within the df to be analyzed.
#' @return the df, with an additional field - "OL_<field>".
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
qc_outlierCheck <- function(df, field) {
  this.field <- df[,field]
  #figure out some outlier groups
  outlier015 <- boxplot.stats(this.field,coef=1.5)$out
  outlier05 <- boxplot.stats(this.field,coef=5)$out
  outlier10 <- boxplot.stats(this.field,coef=10)$out
  #remove overlaps from outlier sets
  outlier05 = setdiff(outlier05, outlier10)
  outlier015 = setdiff(setdiff(outlier015, outlier05), outlier10)
  OL_field <- paste0("OL_",field)
  df[OL_field]<-"OK"
  if (length(outlier10)>0){
    for (i in 1:length(outlier10)){
      df[which(df[field]==outlier10[i]),OL_field]<-"Extreme"
    }
  }
  if (length(outlier05)>0){
    for (i in 1:length(outlier05)){
      df[which(df[field]==outlier05[i]),OL_field]<-"Likely"
    }
  }
  if (length(outlier015)>0){
    for (i in 1:length(outlier015)){
      df[which(df[field]==outlier015[i]),OL_field]<-"Possible"
    }
  }
  return(invisible(df))
}
