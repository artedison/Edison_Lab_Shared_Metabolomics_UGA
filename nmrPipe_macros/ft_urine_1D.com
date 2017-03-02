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


nmrPipe  -in ../fid/$1.fid				\
#| nmrPipe  -fn SOL -po 1					\
| nmrPipe  -fn EM -lb 0.5         \
| nmrPipe  -fn ZF -auto                                 \
| nmrPipe  -fn FT -auto -verb                           \
| nmrPipe  -fn PS -p0 -70.5  -p1 43.0 -di       \
| nmrPipe  -fn POLY -ord 5      \
| nmrPipe -out ../ft/$1.ft -ov -verb			\

echo done


