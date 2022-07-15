#' Predict expression on genetic reference data
#'
#' Given a list of genes and the expression weights for their eQTLs,
#' calculate the matrix of predicted gene expression values on a
#' genetic reference data set.
#'
#' @param genelist A data frame with information on the genes of interest,
#'                 including the path to the data file containing the expression
#'                 weights; typically output from `read_genelist`
#' @param refdata Genetic reference data, as a list with PLINK-style components;
#'                typically output from `read_refdata`
#'
#' @returns A matrix with subjects in the reference data as rows and genes as
#' columns, containing for each subject / gene the predicted expression level
#' (scaled across ubjects).
#'
#' @seealso [read_genelist()], [read_refdata()]
#' @export
predict_expression <- function(genelist, refdata)
{
  nsnps_ref <- nrow(refdata$bim)
  nsubj_ref <- nrow(refdata$bed)
  ngenes    <- nrow(genelist)

  refsnps_keep <- rep(FALSE, nsnps_ref)
  predmat      <- matrix(0, nrow = nsubj_ref, ncol = ngenes)
  ## FIXME: colnames (rownames), do they make sense & are they safe?

  ## Loop over genes
  for (i in 1:ngenes ) {

    ## Load the gene-specific weight information as list
    ww <- load2list( genelist$FILE[i] )

    ## Zero out undefined weights
    ## FIXME: should we report missing values?
    ww$wgt.matrix[is.na(ww$wgt.matrix)] <- 0

    ## Match up the SNPs and weights, iow only keep SNPs relevant
    ## for current gene prediction
    ## FIXME: use variable names instead of col numbers where possible
    m      <- match( ww$snps[,2] , refdata$bim[,2] )
    m.keep <- !is.na(m)
    ## Reduce the weight data to match
    ww$snps <- ww$snps[m.keep, ]
    ww$wgt.matrix = ww$wgt.matrix[m.keep, ]
    refsnps_keep[ m[m.keep] ] = TRUE
    ## Reduce the reference data to match
    cur.genos <- refdata$bed[, m[m.keep] ]
    cur.bim   <- refdata$bim[m[m.keep], ]
    # Flip WEIGHTS for mismatching alleles
    qc <- allele.qc( ww$snps[,5] , ww$snps[,6] , cur.bim[,5] , cur.bim[,6] )
    ww$wgt.matrix[qc$flip, ] <- -1 * ww$wgt.matrix[qc$flip, ]

    # Predict into reference
    ## FIXME: should we test for non-empty? Technically, we could use the
    ## name directly for column selection otherwise
    mod <- which(colnames(ww$wgt.matrix) == genelist$MODEL[i])
    predmat[, i] <- cur.genos %*% ww$wgt.matrix[ , mod ]
  }
  ## Scale gene expression across su
  predmat <- scale( predmat )

  predmat
}
