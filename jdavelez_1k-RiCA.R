#- Custom functions to filter genotype matrix analized in Arbelaez et al. 2019
#- FUNCTION TO ESTIMATE SNP CALL-RATES
snpCR.calc<- function(gMat = M){
  ## Juan David Arbelaez-Velez
  ## 
  ##
  ## A script to calculate SNP call rates from a genotype Matrix (n \times m)
  ## of unphased genotypes for n lines (in columns) and m biallelic markers (in rows)
  ## coded as {-1, 0, 1} and missing data as NA
  ##
  ## Input:
  ## geno: Sample-by-SNP matrix of genotypes {-1,0,1} and missing data as NA.
  ##
  ## Output: 
  ## vector of SNP call-rate values for each SNP
  ##
  sCR<- function(x = y){
    m = length(!is.na(x))/length(x)
    return(m)
  }
  CR.snp<- apply(gMat, 1, sCR)
  return(CR.snp)
}

#- FUNCTION TO ESTIMATE SNP MINOR ALLELE FREQUENCY
maf.calc<- function(gMat = M){
  ## Juan David Arbelaez-Velez
  ## 
  ##
  ## A script to calculate minor allele frequency from a genotype Matrix (n \times m)
  ## of unphased genotypes for n lines (in columns) and m biallelic markers (in rows)
  ## coded as {-1, 0, 1}
  ##
  ## Input:
  ## geno: Sample-by-SNP matrix of genotypes {-1,0,1}.
  ##
  ## Output: 
  ## vector of minor allele frequency values for each SNP
  ##
  maf<- function(x = y){
    m = as.data.frame(table(x))
    m = min(m[,2])/sum(m[,2])
    return(m)
  }
  maf.snp<- apply(gMat, 1, maf)
  return(maf.snp)
}