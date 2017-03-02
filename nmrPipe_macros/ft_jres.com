#!/bin/csh

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
| nmrPipe -fn SP -off 0.0 -end 1 -pow 2.0 -c 1.0 \
| nmrPipe -fn ZF -auto \
| nmrPipe -fn FT -verb \
#| nmrPipe -fn EXT -x1 0% -xn 100% -sw \
| nmrPipe -fn TP \
| nmrPipe -fn LP -fb \
| nmrPipe -fn SP -off 0.0 -end 1 -pow 2.0 -c 1.0 \
| nmrPipe -fn ZF -auto \
| nmrPipe -fn FT -auto -verb \
| nmrPipe -fn MC \
| nmrPipe -fn TP \
 -out spec.ft2 -ov

# Tilt and symmetrize:

set ySize = `getParm -in spec.ft2 -parm NDSIZE -dim CUR_YDIM -fmt %.0f`
set xSW   = `getParm -in spec.ft2 -parm NDSW   -dim CUR_XDIM -fmt %f`
set ySW   = `getParm -in spec.ft2 -parm NDSW   -dim CUR_YDIM -fmt %f`

set mid   = `MATH "integer($ySize/2) + 1"`
set slope = `MATH "$ySW/$xSW"`

echo "J-Tilt. Midpoint: $mid Slope: $slope"

nmrPipe -in spec.ft2 \
| nmrPipe -fn MAC -macro $NMRTXT/jTilt.M -var slope $slope mid $mid \
  -out tilt.ft2 -ov -verb

#proj2D.tcl -in tilt.ft2

###############################################

nmrPipe -in ./tilt.ft2 \
| nmrPipe -fn TP     \
| nmrPipe -fn MAC -macro $NMRTXT/sym.M -var mid $mid -verb \
| nmrPipe -fn TP     \
| nmrPipe -out ../ft/$1.ft -ov -verb			\

echo done


