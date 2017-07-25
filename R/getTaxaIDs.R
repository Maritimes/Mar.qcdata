#' @title getTaxaIDs
#' @description This function attempts to determine the aphiaIDs (WoRMS) and
#' TSNs (ITIS) using first the scientific name, and then the common name. If an
#' aphiaid is found, but no TSN, it will try to find the TSN using the aphiaid.
#' If a TSN is found, but no aphiaid, it will try to find the aphiaid using the
#' TSN.
#' Additionally, in addition to all of the original fields, the returned data
#' will also include the WoRMS aphiaid, the ITIS TSN, and will indicate the
#' method used to find each.  Also for each, a field will indicate whether one,
#' many, or no matches were found.  Lastly, if the species has a different
#' "accepted name" than the one provided, that will also be returned.
#' @param spec_list the dataframe containing information to be decoded to TSN and
#' aphiaIDs
#' @param sci_col the name of the column of the dataframe containing the
#' scientific names
#' @param comm_col the name of the column of the dataframe containing the
#' common names
#' @importFrom utils head
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom bio.datawrangling make_oracle_cxn
#' @importFrom ritis search_scientific
#' @importFrom ritis search_common
#' @importFrom taxize itis_acceptname
#' @importFrom taxize get_wormsid
#' @importFrom taxizesoap worms_extid
#' @importFrom jsonlite fromJSON
#' @family qc
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
getTaxaIDs <- function(spec_list=NULL, sci_col=NULL, comm_col=NULL){

# Functions ---------------------------------------------------------------
  sapply_pb <- function(X, FUN, ...)
  {
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

    wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
    }
    res <- sapply(X, wrapper, ...)
    close(pb)
    res
  }

    proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), perl = TRUE)

    WoRMS2ITIS <- function(searchterm = NULL){
    this = worms_extid(searchterm, type="tsn")
    return(this)
  }

  chk_WoRMS <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    this = get_wormsid(query=searchterm, searchtype=searchtype, ask=ask, verbose=verbose)
    this = as.data.frame(cbind(SEARCHTERM = searchterm, APHIAID=as.vector(this),MATCH_WORMS=attr(this,"match"), METHOD_WORMS = searchtype))
    if (NROW(this[this$MATCH_WORMS=="not found",])>0) this[this$MATCH_WORMS=="not found",]$METHOD_WORMS<-NA
    if (NROW(this[this$MATCH_WORMS=="NA due to ask=FALSE",])>0) this[this$MATCH_WORMS=="NA due to ask=FALSE",]$MATCH_WORMS<-"multi match"
    return(this)
  }

    chk_ITIS <- function(searchterm = NULL, searchtype='scientific', ask=FALSE, verbose=FALSE){
    this = data.frame(SEARCHTERM = NA, TSN = NA, MATCH_ITIS = NA, METHOD_ITIS = NA)
    for (i in 1:NROW(searchterm)){
      if (searchtype == 'scientific'){
        thisi = search_scientific(searchterm[i], wt = "json", raw = FALSE)
      }else{
        thisi = search_common(searchterm[i], wt = "json", raw = FALSE)
      }
      if (nrow(thisi)==0) {
        thisi=data.frame(SEARCHTERM = searchterm[i], TSN = NA, MATCH_ITIS = "not found", METHOD_ITIS = NA)
      }else if (nrow(thisi)==1) {
        thisi = cbind(SEARCHTERM= searchterm[i],TSN = thisi$tsn, MATCH_ITIS = "found", METHOD_ITIS = searchtype)
      }else{
        thisi = cbind(SEARCHTERM = searchterm[i],TSN = NA, MATCH_ITIS = "multi match", METHOD_ITIS = searchtype)
        #could also populate field with vector of matches, but this takes it out
        #of sync with the worms data and some spp get MANY matches (e.g. 100s),
        #TSN = paste0(thisi$tsn, collapse=", ")
      }
      this = rbind(this, thisi)
    }
    this = this[rowSums(is.na(this)) != ncol(this),]
    return(this)
  }

# Processing Steps --------------------------------------------------------
  spec_list[,sci_col] = proper(spec_list[,sci_col])

  # WoRMS
  # check for matches based on sci name, then common name, then combine results
  cat("Checking WoRMS for Scientific names\n")
  chk_WoRMSspec_list_sci = as.data.frame(t(sapply_pb(spec_list[,sci_col],chk_WoRMS)))
  chk_WoRMSspec_list_sci = chk_WoRMSspec_list_sci[!chk_WoRMSspec_list_sci$MATCH=="not found",]
  if (NROW(chk_WoRMSspec_list_sci) < NROW(spec_list)){
    cat("Checking WoRMS for common names\n")
    chk_WoRMSspec_list_comm = as.data.frame(t(sapply_pb(spec_list[!(spec_list[,sci_col] %in% chk_WoRMSspec_list_sci$SEARCHTERM),][,comm_col],chk_WoRMS, searchtype="common")))
    #chk_WoRMSspec_list_comm = chk_WoRMS(searchterm = spec_list[!(spec_list[,sci_col] %in% chk_WoRMSspec_list_sci$SEARCHTERM),][,comm_col], searchtype="common")
  }
  spec_list_WoRMS = merge(spec_list, chk_WoRMSspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM")
  if (exists("chk_WoRMSspec_list_comm")){
    spec_list_WoRMS = rbind(spec_list_WoRMS, merge(spec_list, chk_WoRMSspec_list_comm, by.x = comm_col, by.y = "SEARCHTERM"))
  }
  spec_list_WoRMS = data.frame("SPEC"=spec_list_WoRMS$SPEC,
                               "SPECIES"=spec_list_WoRMS$SPECIES,
                               "COMM"=spec_list_WoRMS$COMM,
                               "APHIAID"=unlist(spec_list_WoRMS$APHIAID),
                               "MATCH_WORMS"=unlist(spec_list_WoRMS$MATCH_WORMS ),
                               "METHOD_WORMS"=unlist(spec_list_WoRMS$METHOD_WORMS))
  # ITIS
  # check for matches based on sci name, then common name
  cat("Checking ITIS for Scientific names\n")
  chk_ITISspec_list_sci = as.data.frame(t(sapply_pb(spec_list[,sci_col],chk_ITIS)))
  chk_ITISspec_list_sci = chk_ITISspec_list_sci[!chk_ITISspec_list_sci$MATCH=="not found",]
  if (NROW(chk_ITISspec_list_sci) < NROW(spec_list)){
    cat("Checking ITIS for common names\n")
    chk_ITISspec_list_comm = as.data.frame(t(sapply_pb(spec_list[!spec_list[,sci_col] %in% chk_ITISspec_list_sci$SEARCHTERM,][,comm_col],chk_ITIS, searchtype="common")))
    #chk_ITISspec_list_comm = chk_ITIS(searchterm = spec_list[!spec_list[,sci_col] %in% chk_ITISspec_list_sci$SEARCHTERM,][,comm_col], searchtype="common")
  }
     spec_list_ITIS = merge(spec_list, chk_ITISspec_list_sci, by.x = sci_col, by.y = "SEARCHTERM")
  if (exists("chk_ITISspec_list_comm")){
    spec_list_ITIS = rbind(spec_list_ITIS, merge(spec_list, chk_ITISspec_list_comm, by.x = comm_col, by.y = "SEARCHTERM"))
  }
     spec_list_ITIS = data.frame("SPEC"=spec_list_ITIS$SPEC,
                                  "SPECIES"=spec_list_ITIS$SPECIES,
                                  "COMM"=spec_list_ITIS$COMM,
                                  "TSN"=unlist(spec_list_ITIS$TSN),
                                  "MATCH_ITIS"=unlist(spec_list_ITIS$MATCH_ITIS ),
                                  "METHOD_ITIS"=unlist(spec_list_ITIS$METHOD_ITIS))

  # merge WoRMS and ITIS results into single dataframe
   spec_list_ID = merge(spec_list_WoRMS, spec_list_ITIS)

   if (requireNamespace("taxizesoap", quietly = TRUE)==FALSE){
     cat("If you had installed taxizesoap from github, the script would now be
          attempting to determine ITIS codes using the aphiaids it has already
          found.  To install if for next time, you would do:\n
          \n
         require(devtools)\n
         install_github('ropensci\\taxizesoap')")
   } else {


  cat("Trying to find missing ITIS IDs using found WoRMS IDs...\n")
  potential_ITIS = spec_list_ID[!is.na(spec_list_ID$APHIAID) & is.na(spec_list_ID$TSN),]
  if (NROW(potential_ITIS)>0){
    spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID,]$TSN<-sapply_pb(potential_ITIS$APHIAID, WoRMS2ITIS)
    if (NROW(spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & spec_list_ID$TSN =="",])>0)
      spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & spec_list_ID$TSN =="",]$TSN <- NA
    if (NROW(spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & !is.na(spec_list_ID$TSN),])>0)
      spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & !is.na(spec_list_ID$TSN),]$METHOD_ITIS <- "from aphiaid"
    if (NROW(spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & !is.na(spec_list_ID$TSN) & spec_list_ID$MATCH_ITIS == "not found",])>0)
      spec_list_ID[spec_list_ID$APHIAID %in% potential_ITIS$APHIAID & !is.na(spec_list_ID$TSN) & spec_list_ID$MATCH_ITIS == "not found",]$MATCH_ITIS <-"found"
  }
   }
  cat("Check that we're using 'accepted' versions of th ITIS IDs\n")
  TSNCheck = data.frame(t(sapply_pb(spec_list_ID[!is.na(spec_list_ID$TSN),]$TSN, itis_acceptname)))
  spec_list_ID = merge(spec_list_ID, TSNCheck[,c("submittedtsn","acceptedtsn","acceptedname")], by.x = "TSN", by.y = "submittedtsn", all.x=TRUE)
  spec_list_ID$TSNFinal = NA
  spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$TSNFinal <- spec_list_ID[!is.na(spec_list_ID$acceptedtsn),]$acceptedtsn
  spec_list_ID[is.na(spec_list_ID$acceptedtsn),]$TSNFinal <- spec_list_ID[is.na(spec_list_ID$acceptedtsn),]$TSN
  spec_list_ID$TSN <- unlist(spec_list_ID$TSNFinal)
  spec_list_ID$TSNFinal <- NULL
  spec_list_ID$acceptedtsn <- NULL
  spec_list_ID$SPEC_SUGGEST <- toupper(spec_list_ID$acceptedname)
  spec_list_ID$acceptedname <- NULL

  # check if we can use any of the tsns we have to get the worms records
  cat("Trying to find missing WoRMS IDs using found ITIS IDs...\n")
  for (i in 1:NROW(spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),])){
    this <- tryCatch(
      {
        fromJSON(paste0("http://www.marinespecies.org/rest/AphiaRecordByExternalID/",spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$TSN,"?type=tsn"))$AphiaID
      },
      error=function(cond){
        #message(cond)
      }
    )
    if (!is.null(this)) {
      spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$MATCH_WORMS <- "from TSN"
      spec_list_ID[is.na(spec_list_ID$APHIAID) & !is.na(spec_list_ID$TSN),][i,]$APHIAID<-this
    }
  }
  #
  return(spec_list_ID)
}
