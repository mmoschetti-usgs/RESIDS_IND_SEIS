RNsubset <- function(x, col.R, col.EQID, maxR, minN){
  # Description:
  # Dates: 12/28/13 - Written by Eric Thompson
  #        10/01/14 - DMB: Corrected mistake in the next to last line
  #        10/16/14 - DMB: Uncorrected the "mistake"
  #        08/12/15 - DMB: Generalize by specifying columns with R, EQID 
  #                   (as used before, I had to modify the colnames when
  #                   using this with different datasets, because the
  #                   colnames weren't always the same).
  # ------------
  # - Function for filtering dataframe to only keep EQs
  #   with at least minN records within maxR
  # - Records are retained with R > maxR for qualifying
  #   events. 
  # - Assumes x is a dataframe, where two of the columns
  #   are EQID and R: 
  #   (NOTE: these might need to be renamed to
  #   agree with the names in the input dataframe "x")
  #
  xR <- x[x[,col.R] <= maxR, ]
  NrR <- table(xR[,col.EQID])
  MID <- match(xR[,col.EQID], names(NrR))
  xR$Nrecs <- NrR[MID]
  xRN <- xR[xR$Nrecs >= minN, ]
  EQIDkeep <- unique(xRN[,col.EQID]) 
  xRNext <- x[!is.na(match(x[,col.EQID], EQIDkeep)), ]
#  xRNext <- xR[!is.na(match(xR[,col.EQID], EQIDkeep)), ]  # The mistake I made; this does NOT retain records with R > maxR
  return(xRNext)
}
