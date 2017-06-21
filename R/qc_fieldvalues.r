#' @title qc_fieldvalues
#' @description This function takes a table, and for every unque value of one
#' field, plots the values of other specified fields against eachother.  For
#' example, you might QC you catch table by species - in which case, for every
#' different species, you might plot the TOTNO vs the TOTWGT.
#' @param qc.df default is \code{NULL}. This identifies the name of the table
#' you want to plot values from.
#' @param field.master default is \code{NULL}.  This identifies the field you
#' want to generates plots by.  In the example above, this would be SPECIES.
#' param field.master.values default is \code{NULL}.  In cases where you only
#' want to QC a few different potential values for\code{field.master}, you can
#' specify them here.
#' @param field.master.values default is \code{NULL}. If you are only interested
#' in a few values of your \code{field.master}, you can enter them as a vector
#' here (e.g. a couple years, or a list of species codes, etc)
#' @param fields.qc default is \code{NULL}. This is a vector of between 2 and 5
#' fields you want to plot against eachother.  If left blank, and there are less
#' than 6 fields, they will all be plotted.
#' @param min.n default is 5.  This is the minimum number of values you need
#' before plots will be done
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom stats complete.cases
#' @importFrom stats lm
#' @importFrom stats loess
#' @importFrom stats na.omit
#' @importFrom utils capture.output
#' @importFrom utils data
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom GGally ggpairs
#' @family qc
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
qc_fieldvalues<-function(qc.df=NULL,
                         field.master = NULL,
                         field.master.values = NULL,
                         fields.qc=NULL,
                         min.n = 5){
  if (class(qc.df) == "character") {
    qc.table.nm = qc.df
    qc.table = get(qc.df)
  }else{
    qc.table.nm = "dataframe"
    qc.table = qc.df
  }
  if (is.null(fields.qc)) fields.qc = colnames(qc.table[, !names(qc.table) %in% field.master])
  if (length(fields.qc)>5) stop("Please select between 2 and 5 fields as 'fields.qc'")
  if (is.null(field.master.values)){
    field.master.values.orig = field.master.values
    field.master.values = unique(qc.table[order(qc.table[field.master]),field.master])
  }else{
    field.master.values.orig = field.master.values
  }

  #paginate tabular results
  paginator<-function(the.data, rec_Start, the.title){
    rowspage =40
    for (i in 1:ceiling((NROW(the.data)-rec_Start)/rowspage)){
      plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
      title(the.title,cex.main=1.5)
      text(0.1,0.35,adj=c(0,0), paste0(capture.output(skipped[rec_Start:min(rowspage,NROW(the.data)),]), collapse='\n'), pos=4, family="mono")
      rec_Start = rec_Start+rowspage
    }
  }

  skipped = data.frame(Value = character(), obs= numeric())
  dropped.cols = c()
  this.pdf = paste0("qc_fieldvalues_",qc.table.nm,"_by_",field.master,".pdf")
  pdf(this.pdf,width=8.5,height=11)
  op<-par(mar=c(0,0,2,0))
  plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
  par(op)

   title(paste0("QC Analysis of ",qc.table.nm),cex.main=1.5)
   text(0.1,0.9,adj=c(0,0),lab=paste0("Plotting values of ",paste(fields.qc, collapse = ", ")," for ",field.master))
   text(0.15,0.85,adj=c(0,0),lab=paste0("Where: "))
   text(0.15,0.8,adj=c(0,0),lab=paste0("There are at least ",min.n," values"))
   if (!is.null(field.master.values.orig)){
     text(0.15,0.75,adj=c(0,0),lab=paste0("And ", field.master," in the following list:"))
     text(0.15,0.7,adj=c(0,0),lab=paste(field.master.values.orig, collapse = ", "))
   }else{
     text(0.15,0.75,adj=c(0,0),lab=paste0("For all values of ",field.master))
   }
   text(0.1,0.55,adj=c(0,0),lab=paste0("qc_fieldvalues(qc.table = '",qc.table.nm,"',
                                      field.master = '",field.master,"',
                                      field.master.values = c(",paste(field.master.values.orig, collapse = ", "),"),
                                      fields.qc=c('",paste(fields.qc, collapse = "', '"),"'),
                                      min.n = ",min.n,")"))
  for (i in 1:length(field.master.values)){
    this = qc.table[qc.table[field.master] == field.master.values[i],fields.qc]
    #drop columns that are entirely NAs (capture names)
    dropped.cols = c(dropped.cols, names(this[sapply(this, function(x) all(is.na(x)))]))
    this[sapply(this, function(x) all(is.na(x)))] <- NULL
    this = this[complete.cases(this),] #can't plot data if one is an NA
    if (NROW(this)>min.n-1) # need more than n rows to do plots
    {
      col_names <- sapply(this, function(col) length(unique(col)) < 7)
      this[ , col_names] <- lapply(this[ , col_names] , factor)
      my_fn <- function(data, mapping, ...){
        p <- ggplot(data = data, mapping = mapping) +
          geom_point() +
          geom_smooth(method=loess, fill="red", color="red", ...) +
          geom_smooth(method=lm, fill="blue", color="blue", ...)
      }
      megaplot = ggpairs(this,
                         title = paste0(field.master,": ",field.master.values[i]),
                         lower = list(continuous = my_fn))
      print(megaplot)
    }else{
      skipped = rbind(skipped,data.frame(Value = field.master.values[i],
                                         obs = nrow(na.omit(qc.table[qc.table[field.master] == field.master.values[i],fields.qc])))
      )
    }
  }

# at the end of the document, identify records that were filtered out, either
# because they don't have enough records, or because they were all NA...
   if (NROW(skipped)>0){
     names(skipped)[1] = field.master
     names(skipped)[2] = "Observations (n)"
     paginator(skipped, 0, "Skipped Values")
   }
   if (NROW(as.data.frame(dropped.cols))>0){
     dropped.cols = as.data.frame(dropped.cols)
     paginator(dropped.cols, 0, "Un-assessed Columns (all values were NA)")
   }
    dev.off()
   cat(paste0("Wrote file to: ",file.path(getwd(),this.pdf)))
}

