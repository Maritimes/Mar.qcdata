#' @title uniqueValuesPerField
#' @description This function takes a dataframe, and for each column, it finds 
#' 1) the unique values that exist, and 2) the number of occurrences of each 
#' value.  This can be useful when trying to find blunders in data entry, where 
#' an odd value was entered.
#' Additionally, the function can attempt to find instances where a certain 
#' number of records were categorized with a particular coding via the 
#' \code{target_n} parameter.  
#' @param df a dataframe to be analyzed.
#' @param cols_ignore the default is \code{NULL}. This is a vector of column
#' names that you want to ignore.  Any columns listed here will not be included
#' in the results.
#' @param target_n the default is \code{NULL}.  This is an integer that will be 
#' used to seek cases where some value is listed in a field this many times.  
#' For example, while comparing data to another source, you might notice that 
#' you are off by 26 records.  By adding \code{target_n=26}, you can return only 
#' those columns which have a parameter that had 26 entries in the data.
#' @param target_n_wiggle the default is \code{0}.  This parameter allows some 
#' wiggle-room on matching the value specified by \code{target_n}.  
#' \code{target_n_wiggle=0} is strict, and only returns exact matches, but 
#' positive values specify a percentage (that is applied to the number of rows)
#' on the supplied \code{df}.  For example, if the dataframe is 1000 rows, the
#' \code{target_n =500}, and the \code{target_n_wiggle =10}, the acceptable 
#' value range will be 400:600 (10%).  However, on a dataframe of 10 rows with 
#' \code{target_n = 5} and the same \code{target_n_wiggle =10}, the acceptable 
#' value range will be 4:6.
#' @return a list of column names from the original df, with the unique 
#' values found in each, as well as the number of records for each unique value.
#' Fields that were found to be either unique (i.e. every value different) or 
#' uniform (i.e. every value the same) will be dropped.
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
uniqueValuesPerField <-function(df,
                      cols_ignore=NULL,
                      target_n = NULL,
                      target_n_wiggle=0){

  rownames(df)<-seq.int(nrow(df))
  intcols=lapply(df, function(x)unique(x))
  uniformCols = names(intcols)[lengths(intcols) == 1]
  idCols = names(intcols)[lengths(intcols) == nrow(df)]
  intcols= intcols[!names(intcols) %in% uniformCols] 
  intcols= intcols[!names(intcols) %in% idCols] 
  if (!is.null(cols_ignore)) {
    intcols = intcols[!names(intcols) %in% cols_ignore]
    df = df[!names(df) %in% cols_ignore]
  }
  df = df[names(df) %in% names(intcols)]
  intIndex<-rownames(df)
  if(!is.null(target_n) & target_n_wiggle>0){
    target_n_wiggle_perc=(target_n_wiggle/100)*nrow(df)
    target_n_max = target_n + target_n_wiggle_perc
    target_n_min = target_n - target_n_wiggle_perc
    cat(paste(
"Wiggle-wiggle:
target_n: ",target_n,"
target_n_wiggle: ", target_n_wiggle, "
df rows: ",nrow(df), "
acceptable value range: ",target_n_min," : ",target_n_max,"\n\n"
))
  }else{
    target_n_max = target_n
    target_n_min = target_n
  }

  res=list()
  for (i in 1:length(names(intcols))){
    these = c(
      with(df, tapply(intIndex, df[names(intcols)[i]], 
                      FUN = function(x) length(unique(x))))
    )
    if(!is.null(target_n)) {
    
      theseTarget = these[these<=target_n_max & these>=target_n_min ]
      if(is.integer(theseTarget) && length(theseTarget) == 0)next()
      theseOthers = length(intcols[intIndex])-sum(theseTarget)
      names(theseOthers)="SumOthers"
      theseNA = length(intcols[intIndex])-(sum(theseOthers)+sum(theseTarget))
      names(theseNA)="SumNA"
      info = c(theseTarget,theseOthers, theseNA)
    }else{
      theseAll = sum(these)
      names(theseAll)="SumNonNA"
      theseNA = length(intcols[intIndex])-sum(these)
      names(theseNA)="SumNA"
      info = c(these,theseAll,theseNA)
    }
    res[[names(intcols)[i]]]=  info
  }
  if(length(idCols)>0) cat(paste0("The following columns were found to be unique for every row:\n",paste0(sort(idCols), collapse=", "),"\n"))
  if(length(uniformCols)>0) cat(paste0("The following columns were found to have a single value:\n",paste0(sort(uniformCols), collapse=", "),"\n"))
  if (length(res)==0 & !is.null(target_n)){
        print("No values found with your specified target")
        return()
  }
  return(res)
}