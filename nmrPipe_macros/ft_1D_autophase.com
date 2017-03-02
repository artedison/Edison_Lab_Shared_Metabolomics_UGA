#!/bin/csh
#
set name = "$1.fid"
set outname = "$1.ft"

nmrPipe -in ./fid/$name \
#| nmrPipe  -fn SOL -po 1                                        \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 1.0      \
#| nmrPipe -fn EM -lb 1.0 -c 1.0 \
| nmrPipe -fn ZF -zf 1 -auto \
| nmrPipe -fn FT \
| nmrPipe -fn MAC -macro $NMRTXT/aph1D.M -all \
-reg -str psX1 9.0ppm psX2 6.0ppm -var xP1 0 \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn MAC -macro $NMRTXT/aph1D.M -all \
-reg -str psX1 9.0ppm psX2 6.0ppm -var xP1 0 \
| nmrPipe -fn MAC -macro $NMRTXT/aphP1_1D.M -all -str \
p0X1 9.0ppm p0X2 6.0ppm p1X1 3.0ppm p1X2 -0.1ppm \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn MAC -macro $NMRTXT/aph1D.M -all \
-reg -str psX1 9.0ppm psX2 6.0ppm -var xP1 0 \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn MAC -macro $NMRTXT/aphP1_1D.M -all -str \
p0X1 9.0ppm p0X2 6.0ppm p1X1 3.0ppm p1X2 -0.1ppm \
| nmrPipe -fn POLY -auto -nl 10.5ppm 10.0ppm 6.2ppm 5.8ppm 1.2ppm -0.5ppm -nw 2 \
| nmrPipe -fn PS -p0 0.0 -p1 0.0 -di \
-out ./ft/$outname -ov


