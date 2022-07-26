#' Map EQTL gene data to reference data
#'
#' Given a list of genes and the weights for their EQTL-SNPs, this function maps
#' the EQTL-SNPs to SNP-level reference data, and predicts all genetic expression
#' components in the expression data.
#'
#' @param genelist A data frame with information on the genes of interest,
#'                 including the path to the data file containing the expression
#'                 weights; typically output from `read_genelist`
#' @param refdata Genetic reference data, as a list with PLINK-style components;
#'                typically output from `read_refdata`
#'
#' @returns A list with two components:
#'
#'     * `refsnps_keep`: a logical vector that indicates for each SNP in the
#'                       reference data whether it has a weight associated in
#'                       the EQTL data.
#'     * `pred_expression`: a numerical matrix with the predicted gene expression
#'                          component for each gene (as columns) for all subjects
#'                          in the reference data (as rows).
#'
#' @details This is a bit of a mixed function, in that the two result components
#' are only loosely related, and would logically be calculated in separate
#' functions. However, calculation of both requires iterating over the per-gene
#' expression weights, which are stored as files, so that would not be very
#' efficient.
#'
#' @seealso [read_genelist()], [read_refdata()]
#' @export
eqtl_to_refdat <- function(genelist, refdata)
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
    ## FIXME: what if weight SNPs are not part of reference? That would
    ## affect predicted expression? Should that not at least be reported?
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
  ## Scale gene expression across subjects / within genes
  predmat <- scale( predmat )

  list(refsnps_keep = refsnps_keep, pred_expression = predmat)
}

#' Condition per-SNP GWAS statistics on per-gene TWAS statistics
#'
#' Given GWAS summary statistics for a chromosome and a corresponding set of
#' TWAS test statistics, this function will calculate adjusted versions of the
#' summary statistics, conditional on the observed TWAS statistics.
#'
#' @param
#' @param opts List of options used for extracting default values for unspecified
#'             arguments
condTWAS_snp <- function(genelist, opts = opts_rtwas$get())
