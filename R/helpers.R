#' Identify incorrectly coded and reversed alleles
#'
#' When comparing effect / reference alleles for vectors of SNPs from two sources,
#' identify a. SNPs that are either not ACGT or strand-ambiguous, b. where effect
#' / reference allele have been reversed between sources.
#'
#' @param a1,a2 vectors of effect / reference allele from first source
#' @param ref1,ref2 vectors of effect / reference alleles from second source
#'                  (implied: the reference, but not necessary)
#'
#' @returns A list with two logical vectors, each of the same length as the
#'          put vectors: `keep` indicating SNPs that have valid coding and are
#'          not strand-ambiguous, and `flip` indicating SNPs with reversed
#'          effect / reference definition between sources.
#'
#' @export
allele.qc = function(a1,a2,ref1,ref2)
{
  ## FIXME: check vector lengths, NAs?
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ## FIXME: de-duplicate code
  ref  = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip

  snp = list()

  ## Do not keep SNPs that are distinct biological variants, but may be just
  ## complementary (strand-ambiguous)
  snp[["keep"]] = !( (a1=="A" & a2=="T") | (a1=="T" & a2=="A") |
                     (a1=="C" & a2=="G") | (a1=="G" & a2=="C") )
  ## Do not keep SNPs that are not clean single-nucleotide variants
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = FALSE
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = FALSE

  ## Flipped variants when there is cross-allele agreement between sources
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

  snp
}
