#
# set compiler FLAGS to debug mode -g by poghosyan@kit.edu
# ------------------------------------------------------------------------
# usage: source resetopt2gflags.sh
# The reason for using command "source" is to have the current shell execute the commands
# Otherwise an export through a shell script will not work because a shell script runs in a child shell process, and only as long as the child shell is running children processes of chld shell would inherit the export
# 
# ------------------------------------------------------------------------
#
	export CXXFLAGS=-g
	export OMPI_CXXFLAGS=-g
	export CFLAGS=-g
	export OMPI_CFLAGS=-g
	export FFLAGS=-g
	export OMPI_FFLAGS=-g
	export FCFLAGS=-g
	export OMPI_FCFLAGS=-g
