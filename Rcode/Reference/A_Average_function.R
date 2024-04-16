# This is a file containing the functions for the spectral analysis

#---------------------------------------------------------
# Compute the averages spectrum based on the 4 replicates 
# (script from the manual by A. Stevens and L Ramirez-Lopez)
#---------------------------------------------------------
SpecAver <- function (x, repNum, ID, thr = 0.03){
  
  # ID is a column array containing the names of the spectra in 'x'
  # Structure of ID: name_pitrRep_depth_spectraRep (pitRep = replicate number of the soil pit; spectraRep = replicate number of the spectrum)
  
  # Sanity check (|| is OR)
  # It is checked if the 'x' input is a matrix with numeric values
  if(!is.matrix(x) || !is.numeric(x)) {stop("the 'x' argument is not a matrix")}
  # It is checked if the 'repNum' input is a number with length 1 (1-9)
  if(!is.numeric(repNum) || !length(repNum)==1) {
    stop("the 'repNum' argument is wrong it must be a numeric vector of length=1")}
  # It is checked if the number of rows is a multiple of repNum
  if(!((nrow(x) %% repNum) == 0)) {
    stop("Check your matrix.
         The number of rows is not multiple of the number of repetitions")}
  # It is checked if the number of ID's is the same as the number of spectra
  if(((nrow(x) == nrow(ID))) == 0) {
    stop("Check your ID's.
         The number of rows in the ID file does is not equal to the number of spectra")}
  
  # Tasks:
  # calculate the number of samples
  n_samples <- nrow(x) / repNum
  # create the (empty) matrices in which the standard deviation and the mean spectrum
  # of each sample will be stored
  rep_means <- matrix(NA, n_samples, ncol(x))
  rep_sds <- matrix(NA, n_samples, ncol(x))
  # Define the sequence of blocks of samples: every columns contains the rownumbers
  # that are replicates of the same soil sample
  blocks <- matrix(seq(from = 1, to = nrow(x)), repNum, n_samples)
  # create a matrix to store the results of the evaluation of the mean
  # standard deviation of each group of 5 samples.
  meanStd <- matrix(NA, n_samples, 1)
  # Create the empty matrix with the new complete names of the spectra
  newID <- matrix(NA, n_samples, 1)
  # Create the empty matrix with the names of the site
  Site <- matrix(NA, n_samples, 1)
  # Create the empty matrix with the number of the replicate
  Rep <- matrix(NA, n_samples, 1)
  # Create the empty matrix with the treatment
  Treatment <- matrix(NA, n_samples, 1)
  # Create the empty matrix with the depth
  Depth <- matrix(NA, n_samples, 1)
  
  
  # loop
  for (i in seq(from = 1, to = n_samples)){
    ith_block <- blocks[ ,i] # Defines the samples to be analyzed based on their row number
    # at each iteration
    rep_means[i, ] <- colMeans(x[ith_block, ])    # The average for all replicates is calculated
    rep_sds[i, ] <- apply(x[ith_block, ],2,sd)    # The standard deviation for all replicates is calculated
    # Use the "if" statement to identify samples with a mean value of standard
    # deviation (of each group of 5 samples) is higher than a
    # predefined threshold.
    std <- (apply(x[ith_block, ],2,sd)) # standard deviation
    mstd <- mean(std) # mean standard deviation of all wavenumbers
    if(mstd > thr){
      meanStd[i, ] <- "higher(check_sample)"}
    else{
      meanStd[i, ] <- "lower"
    }
    # The matrices with the samples names are constructed
    nameRow <-  blocks[1,i]     # The rownumber of the first of the replicates
    Site[i,1] <-  sub("^([^_]*).*$", "\\1", ID[nameRow,1])
    Treatment[i,1] <- sub("^([^_]*_)([^_]*).*$", "\\2",ID[nameRow,1])
    Depth[i,1] <- sub("^([^_]*_)([^_]*_)([^_]*).*$", "\\3",ID[nameRow,1])
    newID[i,1] <- paste(Site[i,1],"_",Treatment[i,1],"_",Depth[i,1],sep = "")
  }
  return (list(SpecMeans = rep_means, SdRep = rep_sds, Check = meanStd, Site = Site, Treatment = Treatment, Depth = Depth, ID = newID))
}

