#' Options for the current RTWAS analysis
#'
#' The FUSION code base on which RTWAS is built uses an option list to control
#' the execution of the component scripts: for specifying input data files,
#' output directories, reference data, type of analysis, window sizes for genetic
#' loci, plot legends etc. RTWAS initializes and maintains a compatible
#' list of options, both defaults and required inputs, but uses only a subset
#' of options itself (e.g. names of input files), while type of analysis or plot
#' are determined by calling the corresponding function rather than setting an
#' option.
#'
#' The management of the option list is built on the `new_defaults`-mechanism
#' implemented in package `knitr`: the list of default options is an R6 object
#' with main methods `get` and `set` to either read or modify the default values.
#' The basic list of options `opts_rtwas` is a global object that lives in the
#' package namespace, though local instances (for multiple concurrent analyses)
#' can be created via calls to `knitr:::new_defaults` and the `$restore`-method.
#'
#' Normally, we would start an analysis by setting the data- and output options
#' in a call to `opts_rtwas$set`, possibly followed by a second call to modify
#' default options for the desired function calls.
#'
#' FIXME: a list of options
#'
#' @export
#' @examples
#' opts_rtwas$get("chr")
#' opts_rtwas$set(chr = 17)
#' opts_rtwas$get("chr")
#' ## Complete list
#' opts_rtwas$get()
opts_rtwas <- knitr:::new_defaults(list(
  ## Data
  sumstats   = NA, ## Summary statistics
  input      = NA, ## List of pre-calculated expression weights
  out        = NA, ## Stub name for output files
  ref_ld_chr = NA, ## Stub name for LD reference files
  chr        = NA, ## Chromosome for current analysis
  ## Inclusions/exclusions
  minp_input      = 1.0,
  max_r2          = 0.90,
  min_r2          = 0.05,
  locus_win       = 100000, ##	[default=100e3]
  max_cz_increase = 1.96,
  zthresh         = FALSE,
  ## Plotting
  plot            = FALSE,
  plot_legend     = NA,
  plot_corr       = FALSE,
  plot_individual = FALSE,
  plot_eqtl       = FALSE,
  plot_scatter    = FALSE,
  ## Alternative tests
  omnibus      = FALSE,
  omnibus_corr = NA,     ## options: top1, blup, bslmm, enet, lasso or best
  eqtl_model   = "top1",
  ldsc         = FALSE,
  ## Output
  report    = FALSE,
  save_loci = FALSE,
  verbose   = 1       ##	0 = nothing, 1 = minimal, 2 = all
))
