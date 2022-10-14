##################################################
# This function takes the text results of NUSE
# and RLE and an affybatch object. It looks for
# the presence of outliers (defined as in the
# default R boxplot command) for both the median
# and the IQR of the two statistics. It creates two
# folders, one for the 'rejected' samples and one for
# the 'doubt' ones. For both an arrayQualityMetrics
# HTML file is created in order to have a clearer
# picture of the problem.
##################################################


##################################################
# is.outlier = function that returns a logical
# vector with TRUE for the presumed outliers and
# FALSE otherwise.
# An outlier is defined as a value that is 1.5*IQR
# times greater than the 3rd quartile or 1.5*IQR
# times less than the 1st Quartile.
##################################################

is.outlier <- function(x) {
  if(!is.numeric(x)) stop("x must be numeric")
  if(any(is.na(x))) {
    warning("NAs detected, I'm removing them")
    x <- x[!is.na(x)]
  }
  m <- median(x)
  iqr <- IQR(x)

  ## m + 2*iqr is the same as 1.5*iqr + 3rd quartile
  ## the same applies for m-2*iqr
  
  return((x < m - 2 * iqr) | (x > m + 2 * iqr))
}

##################################################
# QCfilter is the main function. It uses is.outlier
# to decide whether a value is an outlier for
# nuseFile$median OR nuseFile$IQR and puts a TRUE in
# a logical vector idxNuse. It does the same for RLE
# and creates the vector idxRle. The logical OR of
# the two vectors is taken as the index of the possibly
# bad chips.
# A numeric test is performed to assess if the values
# are in the rejection or in the doubt region, and
# arrayQualityMetrics is called and used separately
# on the two categories.
##################################################
     
QCfilter <- function(plmFit = plmFit,
                     nuseMedianCutoff = 1.10,
                     nuseIQRcutoff = 0.10,
                     rleMedianCutoff = 0.2,
                     rleIQRcutoff = 1.0) {

  library(affy)
  library(affyPLM)
  library(geneplotter)

  if(class(plmFit) != "PLMset") stop("plmFit must belong to the PLMset class")
  
  print("**************************************************")
  print("I'm calculating the NUSE statistsics, be patient")
  print("**************************************************")

  nuseFile <- as.data.frame(t(NUSE(plmFit, type = "stats")))

  print("**************************************************")
  print("I'm calculating the RLE statistsics, be patient")
  print("**************************************************")
  
  rleFile <- as.data.frame(t(RLE(plmFit, type = "stats")))

  idxNuseMedian <- is.outlier(nuseFile$median)
  idxNuseIQR <- is.outlier(nuseFile$IQR)
  
  idxNuse <- idxNuseMedian | idxNuseIQR

  ## Rejection criteria for the NUSE - the thresholds are
  ## absolutely arbitrary

  idxNuseReject <- ((nuseFile$median > nuseMedianCutoff) | (nuseFile$IQR > nuseIQRcutoff))
  
  ## a doubt chip is a chip that is in idxNuse but not in idxNuseReject

  idxNuseDoubt <- idxNuse & (!idxNuseReject)

  ## same thing for RLE
  
  idxRleMedian <- is.outlier(rleFile$median)
  idxRleIQR <- is.outlier(rleFile$IQR)

  idxRle <- idxRleMedian | idxRleIQR

  ## Rejection criteria for the RLE - the thresholds are
  ## again absolutely arbitrary

  idxRleReject <- ((abs(rleFile$median) > rleMedianCutoff) | (rleFile$IQR > rleIQRcutoff))
  idxRleDoubt <- idxRle & (!idxRleReject)
  
  idxSuspects <- idxNuse | idxRle
  idxReject <- idxNuseReject | idxRleReject

  ## doubt chips are those that are 'doubt' for both and
  ## 'reject' for none. This because a chip could be 'doubt'
  ## for NUSE but reject for RLE  
  
  idxDoubt <- ((idxNuseDoubt | idxRleDoubt) & (!idxReject))

  print("**************************************************")
  print(paste("I have flagged as REJECTED", sum(idxReject), "chips"))
  print(paste("and I have flagged as DOUBT", sum(idxDoubt), "chips"))

  
  ##################################################
  ## if there are rejected or doubt chips we proceed
  ## to calling arrayQualityMetrics
  ##################################################

  any.doubt <- as.logical(sum(idxDoubt))
  any.reject <- as.logical(sum(idxReject))

  ##################################################
  ## Diagnostic Plots For The Rejected Samples
  ##################################################
  
  if(any.reject) {
    rejectedSamples <- sampleNames(plmFit)[idxReject]
    print("The diagnostic plots for the REJECTED chips are in the RejectedChips folder")    

    if(!dir.create("RejectedChips")) stop("ERROR: Cannot create the directory")
    if(!dir.create("RejectedChips/MAplots")) stop("ERROR: Cannot create the directory")
    if(!dir.create("RejectedChips/residuals")) stop("ERROR: Cannot create the directory")

    print("***************************************************************************")
    print("I'm generating the diagnostic plots for the REJECTED chips - be patient...")
    print("***************************************************************************")
    
    for(i in seq(along = rejectedSamples)) {

      ## MA plots
      png(paste("RejectedChips/MAplots/", rejectedSamples[i], ".png", sep = ""), 800, 800)
      MAplot(plmFit, which = which(idxReject)[i], plot.method = "smoothScatter")
      dev.off()
      
      ## residual plots
      png(paste("RejectedChips/residuals/", rejectedSamples[i], ".png", sep = ""), 800, 800)
      image(plmFit, which = which(idxReject)[i], type = "resids")
      dev.off()
    }  
  }
  
  ##################################################
  ## Same Plots For The Doubt Samples
  ##################################################
  
  if(any.doubt) {
    doubtSamples <- sampleNames(plmFit)[idxDoubt]

    print("The diagnostic plots for the DOUBT chips are in the DoubtChips folder")    

    if(!dir.create("DoubtChips")) stop("ERROR: Cannot create the directory")
    if(!dir.create("DoubtChips/MAplots")) stop("ERROR: Cannot create the directory")
    if(!dir.create("DoubtChips/residuals")) stop("ERROR: Cannot create the directory")

    print("************************************************************************")
    print("I'm generating the diagnostic plots for the DOUBT chips - be patient...")
    print("************************************************************************")
    
    for(i in seq(along = doubtSamples)) {

      ## MA plots
      png(paste("DoubtChips/MAplots/", doubtSamples[i], ".png", sep = ""), 800, 800)
      MAplot(plmFit, which = which(idxDoubt)[i], plot.method = "smoothScatter")
      dev.off()
      
      ## residual plots
      png(paste("DoubtChips/residuals/", doubtSamples[i], ".png", sep = ""), 800, 800)
      image(plmFit, which = which(idxDoubt)[i], type = "resids")
      dev.off()
    }
  }
}
