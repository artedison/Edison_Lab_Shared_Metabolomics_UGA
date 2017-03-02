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
nmrPipe  -in ascr6.fid				\
#| nmrPipe  -fn POLY -time				\
| nmrPipe  -fn SOL 					\
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto -verb				\
| nmrPipe  -fn PS -p0 -26.0  -p1 19.0 -di 	\
| nmrPipe  -fn TP                                       \
| nmrPipe  -fn LP -fb                                      \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 1.0	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -alt -verb				\
| nmrPipe  -fn PS -p0 -76.0 -p1 166.0 -di	\
#| nmrPipe -fn POLY	-auto				\
| nmrPipe -fn MED		\
| nmrPipe  -fn TP                                       \
| nmrPipe -fn MED					\
#| nmrPipe -fn POLY	-auto				\
| nmrPipe -out ascr6.ft -ov -verb			\



