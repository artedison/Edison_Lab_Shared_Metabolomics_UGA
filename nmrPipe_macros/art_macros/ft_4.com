#!/bin/csh
#
#
nmrPipe  -in ./4_tocsy.fid				\
| nmrPipe  -fn SOL -po 1    \
#| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c .5  	\
| nmrPipe  -fn EM -lb 5 -c 0.5   \
| nmrPipe  -fn ZF -auto 					\
| nmrPipe  -fn FT -auto  				\
| nmrPipe  -fn PS -p0 -43.2  -p1 -25.0  -di 			\
| nmrPipe  -fn TP   \
| nmrPipe  -fn LP -fb   \
#| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c .5	\
| nmrPipe  -fn EM -lb 5 -c 0.5   \
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT  				\
| nmrPipe  -fn PS -p0 1.6 -p1 -3  -di	\
#| nmrPipe -fn POLY	-auto \
| nmrPipe -fn MED					\
| nmrPipe  -fn TP					 \
| nmrPipe -fn MED 					\
#| nmrPipe -fn POLY	-auto				\
| nmrPipe -out ./4_tocsy_ase.ft -ov -verb			\

pipe2xyz -nv -in ./4_tocsy.ft -out ./tocsy_4.nv
