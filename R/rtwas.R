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
#' @param summarystats
#' @param genelist
#' @param refdata
#' @param pred_ge
#' @param pred_ge_corr
#' @param cons_loci
#' @param lambda
#' @param opts List of options used for extracting default values for unspecified
#'             arguments
#'
#' @export
condTWAS_snp <- function(summarystats, genelist, refdata, pred_ge, pred_ge_corr,
                         cons_loci, lambda = 0, opts = opts_rtwas$get())
{
  nloci  <- nrow(cons_loci)
  ngenes <- nrow(genelist)
  nsubj  <- nrow(refdata$fam)
  ## Set the selection threshold for TWAS z-statistics
  zthresh <- set_zthresh(ngenes, opts)
  ## Return object
  ret <- vector("list", length = nloci)

  ## Prepare header for output file, if required
  if ( opts$report ) {
    file.report <- paste0(opts$out, ".report")
    cat( "FILE", "CHR", "P0", "P1", "HIT.GENES", "JOINT.GENES", "BEST.TWAS.P",
         "BEST.SNP.P", "COND.SNP.P", "VAR.EXP\n", sep='\t' , file=file.report )
  }

  ## Loop over all loci
  for (i in 1:nloci) {

    ## Identify genes with expression data overlapping with the current locus
    ge_keep <- genelist$P0 < cons_loci$ends[i] & genelist$P1 > cons_loci$starts[i]

    ## Set up vars for gene selection
    marg_z <- genelist$TWAS.Z   ## The original test statistics, preserved
    cond_z <- genelist$TWAS.Z   ## Comb. adjusted test stat / processing flag
    cond_z[ !ge_keep ] <- 0     ## ... 0 for genes not part of selection / conditioning
    joint_keep <- rep(FALSE, ngenes) ## Indicator for final set of conditioning genes

    ## SKIP the locus if none of the original TWAS statistics is above the
    ## specified threshold
    if ( sum(cond_z^2 > zthresh^2) == 0 ) {
      if( opts$verbose > 0 ) {
        cat( "WARNING: no models in CLUMP ", i , " had an absolute marginal association statistic higher than " ,
             zthresh , ". Skipping\n" , sep='' , file=stderr() )
      }
    ## Continue with gene selection
    } else {
      while ( any(cond_z^2 > zthresh^2) ) {

        ## Add most conditionally significant feature
        ndx_max <- which.max(cond_z^2)
        joint_keep[ ndx_max ] <- TRUE
        if( opts$verbose > 1 ) {
          cat( (genelist$FILE)[ ndx_max ] , " added to model with conditional Z-score of ",
               cond_z[ndx_max] , "\n" , sep='' , file=stderr() )
        }

        ## Invert correlation matrix of predicted gene expression values to for
        ## local imputation of TWAS Z_statistics
        ## FIXME: reference to paper / equation?
        cur_dinv <- solve(pred_ge_corr[joint_keep, joint_keep])
        ## Loop over genes still remaining in the selection set
        for ( ii in which(cond_z != 0) ) {
          ## Gene has been added to conditioning set, remove from selection set
          if ( joint_keep[ii] ) {
            cond_z[ii] <- 0
          ## Gene has to too strong correlation with at least one gene in the
          ## conditioning set, remove from selection set
          } else if ( max( pred_ge_corr[ii, joint_keep]^2 ) > opts$max_r2 ) {
            cond_z[ii] <- 0
          ## Working branch: calculate the adjusted (conditional) TWAS statistic
          ## for the current gene relative to the conditioning set, and update
          ## the vector of conditional z-statistics
          ## FIXME: reference to paper / equation?
          } else {

            pred_ge_corr_local <- pred_ge_corr[ii, joint_keep, drop=FALSE]
            cur_b     <- marg_z[ii] - pred_ge_corr_local %*% cur_dinv %*%
                                      marg_z[joint_keep, drop=FALSE]
            cur_b_var <- 1 - pred_ge_corr_local %*% cur_dinv %*% t(pred_ge_corr_local)

            ## Check: if we get negative variance, we drop the current gene from
            ## the selection set (crude but effective way of dealing with presumed
            ## collinearity)
            if ( cur_b_var <= 0 ) {
              cond_z[ii] <- 0
            } else {
              cond_z[ii] <- cur_b / sqrt( cur_b_var )
            }

          }  ## end of if
        }    ## end of for loop over current selection set
      }      ## end of while loop for genes with TWAS over threshold

      ## We invert the correlation matrix of the final conditioning set to
      ## prepare for adjusting the GWAS z-statistics
      ## We add a stabilizing constant on the diagonal, similar to the original
      ## imputation algorithm described by Gusev et al 2016 - this is hand-tuned
      ## by the user (default: zero, no diagonal fortification)
      nkeep <- length(which(joint_keep))
      cur_dinv <- solve( pred_ge_corr[joint_keep, joint_keep] + lambda*diag(nkeep) )

      ## Select SNPs to be retained in current locus
      cur_keep_snp = refdata$bim[, "BP"] > cons_loci$starts[i] &
                     refdata$bim[, "BP"] < cons_loci$ends[i]
      snp_z  <- summarystats$Z[ cur_keep_snp ]
      ## FIXME: why do we recalculate this?
      snp_pv <- 2*(pnorm( abs(snp_z) , lower.tail=FALSE))

      ## Covariance between observed SNPs and predicted gene expression across
      ## reference subjects, resulting in SNP x gene matrix
      ## Note that both contributing matrices are still scaled appropriately
      ##
      ## FIXME: reference?
      snp_ge_ld = t(refdata$bed[, cur_keep_snp]) %*% pred_ge[, joint_keep] / ( nsubj - 1 )

      ## TODO (optional) : impute any missing z-scores

      ## Calculate adjusted GWAS z-stats & standard errors
      ## FIXME: reference?
      snp_cond_b   <- snp_z - snp_ge_ld %*% cur_dinv %*% marg_z[joint_keep]
      snp_cond_var <- 1 - diag( snp_ge_ld %*%  cur_dinv  %*% t(snp_ge_ld) )
      ## Check: invalid variances?
      ndx_neg_var <- which(snp_cond_var <= 0)
      if (any(ndx_neg_var)) {
        if( opts$verbose > 0 ) {
          cat( "WARNING: ", length(ndx_neg_var), " SNPs in CLUMP ", i ,
               " had a negative variance for their adjusted test statistics.",
               "Their statistics will be set to NaN\n", sep='' , file=stderr() )
        }
      }
      snp_cond_se <- sqrt(snp_cond_var)

      ## Zero out SNPs with > max_r2 LD to a gene
      autocor <- apply(snp_ge_ld^2, 1, max, na.rm=TRUE) > opts$max_r2
      if (any(autocor)) {
        if( opts$verbose > 0 ) {
          cat( "WARNING: ", length(which(autocor)), " SNPs in CLUMP ", i ,
               " had squared correlation with predicted expression >", opts$max_r2,
               ". Their statistics will be set to zero\n", sep='' , file=stderr() )
        }
      }
      snp_cond_b[ autocor ]  <- 0
      snp_cond_se[ autocor ] <- 1

      ## Test statistic and p-value for adjusted test statistics
      snp_cond_z <- snp_cond_b / snp_cond_se
      snp_cond_p <- 2*(pnorm( abs(snp_cond_z) , lower.tail=FALSE))

      ## Some summary stats for the current locus
      ##if( opts$verbose > 0 ) {
      ##  cat("Locus " , i , ":\n", sep='', file=stderr())
      ##  cat( "\tbest GWAS Chisq\t" , round(max(snp_z^2, na.rm=TRUE), 2) , '\n' ,
      ##       sep='', file=stderr())
      ##  cat( "\tbest GWAS Chisq conditioned\t" , round(snp_cond_z[which.max(snp_z^2)]^2, 2),
      ##       '\n' , sep='', file=stderr() )
      ##  cat( "\tbest conditioned Chisq\t" , round(max(snp_cond_z^2, na.rm=TRUE), 2) ,
      ##       '\n' , sep='', file=stderr() )
      ##}

      # generate before/after manhattan plot
      # THIS IS WHERE THE OPT PLOT GOES

      ## Create results for return (redundancy for possible reporting)
      ## Start with results for SNPs
      snp_res <- data.frame("SNP" = refdata$bim[cur_keep_snp, 2],
                            "POS" = refdata$bim[cur_keep_snp, 4],
                            "GWAS.Z" = snp_z , "GWAS.P" = snp_pv,
                            "GWAS_cond.Z" = snp_cond_z, "GWAS_cond.P" = snp_cond_p)
      ## FIXME: derive from snp_res
      snp_res2 <- data.frame("SNP" = refdata$bim[cur_keep_snp,2],
                             "POS" = refdata$bim[cur_keep_snp,4],
                             "GWAS.LOGP" = -log10(snp_pv) ,
                             "GWAS_cond.LOGP" = -log10(snp_cond_p) )
      ## Gene results
      gene_res <- genelist[ ge_keep, ]
      gene_res$JOINT <- joint_keep[ ge_keep ]
      ## Compute correlation of each gene to the top SNP
      gene_res$TOP.SNP.COR <- round( t( pred_ge[, ge_keep] ) %*%
                                      (refdata$bed[, cur_keep_snp])[, which.max(snp_z^2), drop=FALSE]/ (nsubj - 1 ), 2)

      ret[[i]] <- list(snp_res = snp_res, gene_res = gene_res)

      ## Conditional reporting
      if ( opts$save_loci ) {
        write.table( format(snp_res, cdigits=2), quote=FALSE, col.names=TRUE,
                     row.names=F , sep='\t' , file=paste0(opts$out,".loc_",i,".cond") )
      }
      if ( opts$report ) {
        ## SNP statistics for the current locus
        write.table( format(snp_res2, digits=2, trim=TRUE) , quote=FALSE,
                     col.names=TRUE, row.names=FALSE, sep=',' ,
                     file=paste0(opts$out,".loc_",i,".cond.csv") )
        ## Gene statistics for the current locus
        write.table( gene_res, quote=FALSE , row.names=FALSE , col.names=TRUE ,
                     sep='\t' , file=paste0(opts$out,".loc_", i, ".genes")  )
        ## Summary statistics for the current locus
        ## FIXME: build df incrementally, drop as file after loop
        best.snp.chisq <- max(snp_z^2, na.rm=TRUE)
        cond.snp.chisq <- snp_cond_z[which.max(snp_z^2)]^2
        cat( paste(opts$out, ".loc_", i, sep=''), opts$chr,
             range( refdata$bim[ cur_keep_snp , 4 ] ),
             length(unique(genelist$ID[ge_keep])) , sum( joint_keep[ ge_keep ] ),
             min(genelist$TWAS.P[ ge_keep ], na.rm=TRUE),
             2*pnorm(sqrt(best.snp.chisq), lower.tail=FALSE),
             2*pnorm(sqrt(cond.snp.chisq), lower.tail=FALSE),
             1 - cond.snp.chisq / best.snp.chisq , '\n' , sep='\t' ,
             file=file.report , append=TRUE)
      }
    } ## End of branch: are there sufficiently large TWAS z-stats to condition on?
  } ## End of locus-loop

  ret

}

#' Run a conditional TWAS analysis for a single chromosome
#'
#' This is a wrapper function that combines all necessary function calls to
#' perform a post-TWAS conditional analysis (at the SNP/locus level).
#'
#' @param chr Number of chromosome to analyse (see Details).
#' @param fn_genelist Name of file listing genes with TWAS statistics for
#'                    conditional analysis
#' @param fn_sumstats Name of file listing the GWAS summary statistics to be
#'                    conditioned
#' @param fn_output  Name of file for output generated by this function; note that
#'                  this is a stub, potentially shared by multiple file names
#'                  differing by an extension of the stub
#' @param lambda Ridge parameter for calculating the conditional GWAS statistics
#'               (experimental); default zero, i.e. same behavior as FUSION
#' @param opts List of options used for extracting default values for unspecified
#'             arguments
#'
#' @details If file names (options `fn_`) are not specified, the function falls
#' back on the corresponding elements of the global option argument, namely
#' `ìnput`, `sumstats` and `out` (see `opts_rtwas`); if these are also not set,
#' the function will fail with an error.
#'
#' Note that strictly speaking, the `chr`-argument is redundant, as this is
#' implicitly determined by specifying chromosome-specific input files; this is
#' more of a fallback for situations where the either the gene list- or summary
#' statistics files are not single-chromosome.
#'
#' @export
#' @returns The same list of results as `condTWAS_snp`(the worker function).
#' @seealso [condTWAS_snp()] [opts_rtwas]
postTWAS_conditional_chr <- function(chr, fn_genelist, fn_sumstats, fn_output,
                                     lambda = 0, opts = opts_rtwas$get())
{
  ## Check and fill in as required
  stopifnot( !missing(chr) )
  chr <- round(chr)
  stopifnot( (1<=chr) & (chr <=22) )

  ## Check that a valid gene list file exists
  if (missing(fn_genelist)) fn_genelist <- opts[["input"]]
  stopifnot( file.exists( fn_genelist ))
  genelist <- read_genelist(fn_genelist, chr = chr, opts = opts)

  ## Check that valid reference data has been specified (via global options)
  stopifnot( !is.na( opts[["ref_ld_chr"]]) )
  generef <- read_refdata(chr = chr, opts = opts)

  ## Check that a valid summary statistics file exists
  if (missing(fn_sumstats)) {
    stopifnot( !is.na( opts[["sumstats"]]) )
  } else {
    stopifnot( is.na( opts[["sumstats"]] ) )
    opts[["sumstats"]] <- fn_sumstats
  }
  sumstat <- read_sumstats(generef, opts = opts)

  ## Set the output stub
  if (!missing(fn_output)) {
    stopifnot( is.na(opts[["out"]]) )
    opts[["out"]] <- fn_output
  }

  ## Construct the matrix of predicted gene expressions, and the vector that indicates
  ## which reference SNPs are part of the EQTL
  ## Known as "ge_g.matrix" and "geno.keep" in FUSION
  expr2ref <- eqtl_to_refdat(genelist, generef)

  ## Correlation matrix of predicted gene expressions
  ## Known as "ge_g.ld" in FUSION
  pred_ge_corr <- calcCorr_predGE( expr2ref$pred_expression , opts = opts)

  ## Calculate the contiguous loci
  ## Known as "cons.loc.starts" and "cons.loc.end" in FUSION
  cons.loc <- find_loci(generef, expr2ref$refsnps_keep, opts = opts)

  ## Do the conditional analysis
  cond_snp <- condTWAS_snp(sumstat, genelist, generef, expr2ref$pred_expression,
                           pred_ge_corr, cons.loc, lambda = lambda, opts = opts)
  cond_snp$Call <- match.call()

  cond_snp
}


#' Run a genome-wide conditional TWAS analysis f
#'
#' This is a wrapper function that combines all necessary function calls to
#' perform a post-TWAS conditional analysis (at the SNP/locus level) across
#' a set of chromosomes (by default, all chromosomes 1-22, i.e. genome-wide).
#'
#' @param fn_genelist Name of file listing genes with TWAS statistics for
#'                    conditional analysis; see Details
#' @param fn_sumstats Name of file listing the GWAS summary statistics to be
#'                    conditioned; see Details
#' @param fn_output  Name of file for output generated by this function; note that
#'                  this is a stub shared by multiple output file names
#'                  differing by an extension of the stub
#' @param chr Numbers of chromosomes to analyse, by default `1:22`
#' @param lambda Ridge parameter for calculating the conditional GWAS statistics
#'               (experimental); default zero, i.e. same behavior as FUSION
#' @param opts List of options used for extracting default values for unspecified
#'             arguments
#'
#' @details File names specifying where to read the genes with a TWAS statistic
#'   of interest (`fn_genelist`), the corresponding GWAS summary statistics
#'   (`fn_sumstats`) and a stub for output file names (`fn_output`) need to
#'   be provided via the named arguments, whereas the location of the reference
#'   data need to be specified via `opts_rtwas`.
#'
#'   There is some flexbility in how TWAS- and GWAS statistics can be
#'   specified, however:
#'
#'   1. As character vectors: this puts the responsibility for naming files
#'      completely on the user - they need to provide matching vectors of
#'      fully specified `genelist`- and `sumstats` file names that correspond
#'      to the list of chromosomes to be analysed.
#'
#'   2. As character expressions of length 1: these are interpreted as file name
#'      _templates_: the list of file names will be constructed by replacing
#'      the special string `!!!!` with the numbers of the chromosomes specified
#'      in argument `chr`.
#'
#'   Note that the same special string can (and generally _should_) be used for
#'   the file stub for output file names. All copies of `!!!!` will be replaced,
#'   so this can be used to e.g. have both chromosome-specific file names and
#'   chromosome-specific subfolders.
#'
#' @export
#' @returns A list of per-chromosome post-conditional TWAS results, as returned
#'   by `postTWAS_conditional_chr` (which are essentially the same results as
#'   returned by `condTWAS_snp`).
#' @seealso [postTWAS_conditional_chr()] [condTWAS_snp()] [opts_rtwas]
postTWAS_conditional_gw <- function(fn_genelist, fn_sumstats, fn_output,
                                    chr=1:22, lambda = 0, opts = opts_rtwas$get())
{
  chr <- unique(chr)

  fn_genelist <- expand_special(fn_genelist, chr)
  fn_sumstats <- expand_special(fn_sumstats, chr)
  fn_output   <- expand_special(fn_output, chr)

  ## FIXME: check these all agree as they should

  nchr <- length(chr)
  ret <- vector("list", length = nchr)
  for (i in 1:nchr) {

    ret[[i]] <- postTWAS_conditional_chr(chr = chr[i],
                                         fn_genelist = fn_genelist[i],
                                         fn_sumstats = fn_sumstats[i],
                                         fn_output   = fn_output[i],
                                         lambda = lambda, opts = opts)
  }
  names(ret) <- chr
  ret$Call <- match.call()

  ret
}
