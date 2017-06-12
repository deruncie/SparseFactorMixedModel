# Write the phenotype data to a file in the format used by GEMMA. Each
# line of the file contains one phenotype observation.
write.gemma.pheno <- function (file, phenotype, pheno) {
  y <- pheno[[phenotype]]
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

# ----------------------------------------------------------------------
# Write the covariate data to a file in the format used by GEMMA. Each
# line corresponds to a sample. We must include an additional
# covariate for the intercept.
write.gemma.covariates <- function (file, covariates, pheno) {
  if (length(covariates) == 0) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  } else {
    round.col <- function (x) {
      if (is.numeric(x))
        round(x,digits = 6)
      else
        x
    }
    write.table(cbind(1,data.frame(lapply(pheno[covariates],round.col))),
                file,sep = " ",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
  }
}

# ----------------------------------------------------------------------
# Write the SNP information to a space-delimited text file in the
# format used by GEMMA. This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.
write.gemma.map <- function (file, map)
  write.table(map[c("snp","pos","chr")],file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# ----------------------------------------------------------------------
# Store the mean genotypes as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three column give the SNP label, and
# the two alleles.
write.gemma.geno <- function (file, geno, map) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(map[c("snp","A1","A2")],geno)
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

# ----------------------------------------------------------------------
# Reads in the association results from GEMMA, and returns a data
# frame containing three columns: chromosome number ("chr"); base-pair
# position ("pos"); and the negative base-10 logarithm of the p-value
# ("log10p").
read.gemma.assoc <- function (file) {
  gwscan <- read.table(file,sep = "\t",header = TRUE,check.names = FALSE,
                       quote = "",stringsAsFactors = FALSE)
  rownames(gwscan) <- gwscan$rs
  gwscan           <- gwscan[c("chr","ps","p_lrt")]
  gwscan           <- transform(gwscan,p_lrt = -log10(p_lrt))
  colnames(gwscan) <- c("chr","pos","log10p")

  return(gwscan)
}

# ----------------------------------------------------------------------
# This function maps QTLs using GEMMA without fully accounting for
# population structure. Here we use a block-diagonal relatedness
# matrix in which the blocks correspond to the genetically identical
# mice from the same strain.
run.gemma.norr <- function (phenotype, covariates, pheno, geno, map,
                            gemmadir, gemma.exe) {

  # I add this to the diagonal of the kinship matrix to make sure that
  # calculations involving this matrix are stable.
  delta <- 0.001

  # Set the local directory to the location of the GEMMA files.
  srcdir <- getwd()
  setwd(gemmadir)

  # Give summary of analysis.
  n <- nrow(pheno)
  cat("Mapping QTLs for",phenotype,"in",n,"mice, ")
  if (!is.null(covariates)) {
    cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
  } else {
    cat("with no covariates included.\n")
  }

  # Write the phenotype and covariate data to separate files.
  cat("Writing phenotype and covariate data to file.\n")
  write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno)
  write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),covariates,pheno)

  # Write out the mean genotypes and map information for all markers.
  cat("Writing SNP and genotype data to file.\n")
  write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno,map)
  write.gemma.map(paste0(gemmadir,"/map.txt"),map)

  # Write out the kinship matrix to file.
  cat("Writing block-identity kinship matrix to file.\n");
  K <- repmat(matrix(as.numeric(pheno$strain),n,1),1,n)
  K <- (K == t(K)) + diag(delta,n)
  write.table(round(K,digits = 5),paste0(gemmadir,"/kinship.txt"),sep = " ",
              quote = FALSE,row.names = FALSE,col.names = FALSE)

  # Now we are finally ready to run GEMMA for all markers.
  cat("Computing p-values for",nrow(map),"candidate markers.\n")
  system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
               "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
               "-lmin 0.01 -lmax 100"),
         ignore.stdout = TRUE)

  # Load the results of the GEMMA association analysis.
  gwscan <- read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))
  class(gwscan) <- c("scanone","data.frame")

  # Restore the working directory.
  setwd(srcdir)

  # Return the genome-wide scan.
  return(gwscan)
}

# ----------------------------------------------------------------------
# This function maps QTLs using GEMMA, writing all the files required
# by GEMMA to the directory specified by "gemmadir". The QTLs are
# mapped separately for each chromosome, in which the kinship matrix
# is computed using all markers except the markers on the given
# chromosome.
run.gemma <- function (phenotype, covariates, pheno, geno, map,
                       gemmadir, gemma.exe, K = NULL) {

  # I add this to the diagonal of the kinship matrix to make sure that
  # calculations involving this matrix are stable.
  delta <- 0.001

  # Set the local directory to the location of the GEMMA files.
  srcdir <- getwd()
  setwd(gemmadir)

  # Take care of optional kinsip matrix input.
  compute.kinship.matrix <- is.null(K)

  # Give summary of analysis.
  n <- nrow(pheno)
  cat("Mapping QTLs for",phenotype,"in",n,"mice, ")
  if (!is.null(covariates)) {
    cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
  } else {
    cat("with no covariates included.\n")
  }

  # Write the phenotype and covariate data to separate files.
  cat("Writing phenotype and covariate data to file.\n")
  write.gemma.pheno(paste0(gemmadir,"/pheno.txt"),phenotype,pheno)
  write.gemma.covariates(paste0(gemmadir,"/covariates.txt"),covariates,pheno)

  # Set up the data structures used to store the results of the QTL
  # mapping.
  chromosomes  <- levels(map$chr)
  scans        <- vector("list",length(chromosomes))
  pve          <- rep(NA,length(chromosomes))
  names(scans) <- chromosomes
  names(pve)   <- chromosomes

  # Repeat for each chromosome.
  for (chr in chromosomes) {

    # Compute the kinship matrix, if necessary.
    cat("Mapping QTLs on chromosome ",chr,".\n",sep="")
    if (compute.kinship.matrix) {
      cat(" * Computing kinship matrix.\n")
      markers <- which(map$chr != chr)
      K <- tcrossprod(center.columns(geno[,markers])) / length(markers)
      K <- K + diag(delta,n)
    }

    # Save the kinship matrix to a file.
    cat(" * Writing kinship matrix to file.\n")
    write.table(round(K,digits = 5),paste0(gemmadir,"/kinship.txt"),sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)

    # Write out the mean genotypes and map information for all markers
    # on the chromosome.
    markers <- which(map$chr == chr)
    cat(" * Writing to file genotypes for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    write.gemma.geno(paste0(gemmadir,"/geno.txt"),geno[,markers],map[markers,])
    cat(" * Writing genetic map for",length(markers),"markers on chromosome",
        chr,"to file.\n")
    write.gemma.map(paste0(gemmadir,"/map.txt"),map[markers,])

    # Now we are finally ready to run GEMMA for all markers on the
    # chromosome using the kinship matrix computed using all the
    # markers *not* on the chromosome.
    cat(" * Computing p-values for ",length(markers),
        " markers on chromosome ",chr,".\n",sep="")
    system(paste(gemma.exe,"-g geno.txt -a map.txt -p pheno.txt",
                 "-c covariates.txt -k kinship.txt -notsnp -lmm 2",
                 "-lmin 0.01 -lmax 100"),
           ignore.stdout = TRUE)

    # Load the results of the GEMMA association analysis.
    scans[[chr]] <-
      read.gemma.assoc(paste0(gemmadir,"/output/result.assoc.txt"))

    # Get the estimate of the proportion of variance explained by the
    # polygenic effects on all chromosomes except the current
    # chromosome.
    out      <- scan(paste0(gemmadir,"/output/result.log.txt"),
                     what = "character",sep = " ",quiet = TRUE)
    pve[chr] <- out[which(out == "pve") + 7]
  }

  # Restore the working directory.
  setwd(srcdir)

  # Merge the mapping results from all chromosomes into a single table.
  gwscan           <- do.call(rbind,scans)
  rownames(gwscan) <- do.call(c,lapply(scans,rownames))
  class(gwscan)    <- c("scanone","data.frame")

  # Return the genome-wide scan and the estimates of the proportion of
  # variance explained for each chromosome.
  return(list(gwscan = gwscan,pve = pve))
}

