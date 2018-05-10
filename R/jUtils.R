#' Simulate random points from a genome
#'
#' Weighted sampling by chromosome length
#' @param n Number of points to sample
#' @param genome Seqinfo object with genome parameters
#' @return GRanges with \code{n} simulated points (width=1)
#' @note This is much faster than \code{gr.rand}, but can only simulate points. (33k points / sec)
gr.randpoints <- function(n, genome) {

  if (!is(genome, "Seqinfo")) 
    genome = seqinfo(genome)
  sl = seqlengths(genome)
  available = si2gr(genome)

  ## sample a single linear "chromosome"
  cs <- cumsum(as.numeric(sl))
  pos <- sample(sum(as.numeric(sl)), n)

  # get the chromosome
  s <- sapply(pos, function(i) which(ceiling(i / cs) == 1)[1])

  return(GRanges(seqlevels(genome)[s], IRanges(pos - c(0,cs)[s], width=1), seqinfo=genome))

}
