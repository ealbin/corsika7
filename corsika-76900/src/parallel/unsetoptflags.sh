#
# unsetoptflags.sh:
# =================
#          erase high compiler optimization flags for parallel corsika.
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# usage: . unsetoptflags.sh
# ------------------------------------------------------------------------
#
  flagge=OMPI_CXXFLAGS
  export -n $flagge
  flagge=CXXFLAGS
  export -n $flagge
#
  flagge=OMPI_CFLAGS
  export -n $flagge
  flagge=CFLAGS
  export -n $flagge
#
  flagge=OMPI_FFLAGS
  export -n $flagge
  flagge=FFLAGS
  export -n $flagge
#
  flagge=OMPI_FCFLAGS
  export -n $flagge
  flagge=FCFLAGS
  export -n $flagge
#
