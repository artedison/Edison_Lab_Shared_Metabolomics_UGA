#!/bin/csh
#
#
nmrPipe  -in ./4_tocsy.fid				\
| nmrPipe  -fn SOL -po 1    \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 	\
| nmrPipe  -fn ZF -size 4096 \
| nmrPipe  -fn FT -auto -verb				\
| nmrPipe  -fn PS -p0 -43.2  -p1 -25.0 -di 			\
| nmrPipe -fn POLY -auto \
| covNMR -N1 512 -N2 4096 -QUAD1 1 -CALC 1 			\
| nmrPipe -fn POLY -auto \
| nmrPipe -out ./cov_4_toc.ft -ov -verb                    \

