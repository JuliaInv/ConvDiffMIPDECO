#!/bin/bash

##############################################################################
# runFile.sh: solve AMPL FDM instance
##############################################################################
ampl SrcInvertFDM.mod N08.ampl SrcInvertFDM200N10.dat solver.ampl PlotUW.ampl &> amplFDM08.txt 
#ampl SrcInvertFDM.mod N16.ampl SrcInvertFDM200N10.dat solver.ampl PlotUW.ampl &> amplFDM16.txt 
#ampl SrcInvertFDM.mod N32.ampl SrcInvertFDM200N10.dat solver.ampl PlotUW.ampl &> amplFDM32.txt 

