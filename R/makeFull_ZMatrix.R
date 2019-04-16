###### Function to create full Z-matrix to compute prior ######



# #' Create full Z-matrix for prior computation
# #'
# #' From a list of significant studies, create the full Z-matrix that can be used to compute
# #' the prior.
# #'
# #' @inheritParams bGWAS
# #' @param studies The IDs of the significant studies selected by \code{identify_studiesMR()}
# #'        (numeric vector)
# #'
# #' @return An object containing Log file and pruned Z-Matrix of MR instrument (data table) + create a file if saveFiles=T
# #'




makeFull_ZMatrix <- function(studies=NULL, GWASData,  Z_matrices="~/Z_matrices", save_files=F, verbose=F) {
  platform = .Platform$OS.type
  if(platform=="windows") stop("Windows is not supported yet", call. = FALSE)

  Log = c()
  tmp = paste0("# Loading the ZMatrix... \n")
  Log = update_log(Log, tmp, verbose)


  tmp = paste0("Selecting studies :\n")
  Log = update_log(Log, tmp, verbose)

  if(platform == "unix") {
    ZMatrix=data.table::fread(cmd=paste0("zcat < ",paste0(Z_matrices, "/ZMatrix_Imputed.csv.gz")), select=c(1:5, studies+5 ), showProgress = FALSE, data.table = F)
  }
  tmp = paste0(ncol(ZMatrix)-5, " studies \n")
  Log = update_log(Log, tmp, verbose)

  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific = F), " SNPs \n")
  Log = update_log(Log, tmp, verbose)



  # Add conventional GWAS column, at the end (make sure alleles are aligned)
  GName = attributes(GWASData)$GName
  tmp = paste0("# Adding data from the conventional GWAS : \n \"", GName,
               "\" \n")
  Log = update_log(Log, tmp, verbose)
  
  SNPID = match(colnames(GWASData),c("snpid", "snp", "rnpid", "rs", "rsid"))
  SNPID = which(!is.na(SNPID))[1]
  ALT = match(colnames(GWASData),c("a1", "alts", "alt"))
  ALT = which(!is.na(ALT))[1]
  REF = match(colnames(GWASData),c("a2", "ref", "a0"))
  REF = which(!is.na(REF))[1]
  ZSTAT = match(colnames(GWASData),c("z", "Z", "zscore", "z_obs"))
  ZSTAT = which(!is.na(ZSTAT))[1]
  
  # keep the SNPs in our Z matrix and order them correctly
  GWASData = GWASData[match(ZMatrix$rs, GWASData[,SNPID]),]
  # check alignment
  aligned = which(GWASData[,ALT] == ZMatrix$alt &
                    GWASData[,REF] == ZMatrix$ref)
  swapped = which(GWASData[,REF] == ZMatrix$alt &
                    GWASData[,ALT] == ZMatrix$ref)
  
  GWASData$myZ = NA
  GWASData[aligned, "myZ"] =  GWASData[aligned, ZSTAT]
  GWASData[swapped, "myZ"] = -GWASData[swapped, ZSTAT]
  
  
  
  if(is.null(GName)){
    ZMatrix[, "outcome"] = GWASData$myZ
  } else {
    ZMatrix[, GName] = GWASData$myZ
  }  
  
  tmp = "Done! \n"
  Log = update_log(Log, tmp, verbose)
  
  



  ZMatrix = ZMatrix[complete.cases(ZMatrix),]
  tmp = paste0(format(nrow(ZMatrix), big.mark = ",", scientific=F), " SNPs in common between prior studies and the conventional GWAS \n")
  Log = update_log(Log, tmp, verbose)

  res=list()
  res$log_info = Log
  res$mat = ZMatrix
  return(res)
}

