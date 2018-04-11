#' Get the SNP probs nearest to the given chr and position.
#'
#' chr: character indicating the chromosome to use. This must match
#'      one of the names of genoprobs and map.
#' pos: numeric value in Mb indicating the position on the chromosome.
#' genoprobs: list containing 3-dimensional arrays of allele probs, one 
#'            for each chromosome. In qtl2 format.
#'            Samples in rows, 8 founders in columns, markers in slices.
#' map: list containing numeric vectors of marker positions in Mb.
#' query_fxn: function for querying the CC_SNP data base. Obtain this
#            by calling qtl2:::create_variant_query_func().
#' @return an n x 1 numeric matrix.
#' @examples
#' query_fxn = qtl2::create_variant_query_func(cc_dbfile, filter = "type=='snp'")
#' get_snp_probs = function(chr = "1", pos = 5.0, genoprobs = genoprobs, map = map, query_fxn = query_fxn)
#' @export
get.snp.probs = function(chr, pos, genoprobs, map, query_fxn) {
  
  # Get genoprobs closest to the requested position.
  gp = get.allele.probs(chr, pos, genoprobs, map)
  
  # Get SNPs in a window around the requested position.
  variants = query_fxn(chr, pos - 0.25, pos + 0.25)
  
  # Get the SNP closest to the requested position.
  wh = which.min(abs(variants$pos - pos))
  variants = variants[wh,]
  
  strains = c("A_J",       "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ",
              "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ",     "WSB_EiJ")
  alleles = as.numeric(variants[,strains])
  sp = t(tcrossprod(alleles, gp)) - 1.0
  
  return(sp)
  
} # get.snp.probs()

