#!/bin/csh
#
#

# QUICK REFERENCE
# -in      Read in the fid file
# -fn SP   Window, (watch first and last points)
# -fn ZF   Zero Fill (Double Size, Round Up to Power of 2)
# -fn FT   Fourier Transform
# -fn PS   Phase Correction
# -fn TP   Transpose
# -fn LP   Linear Prediction
# -fn EXT  Extract Left Half
# -fn MED  Median Baseline Filter
# -fn POLY Polynomial Baseline Filter 
#
#
nmrPipe  -in ascr6_1D.fid				\
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto -verb				\
| nmrPipe  -fn PS -p0 -204.0  -p1 14.0 -di 	\
| nmrPipe -out ascr6_1D.ft -ov -verb			\



