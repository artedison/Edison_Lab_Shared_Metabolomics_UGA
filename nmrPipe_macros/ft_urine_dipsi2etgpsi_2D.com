#!/bin/csh
#
#

nmrPipe  -in ../fid/$1.fid				\
| nmrPipe  -fn SOL -po 1    \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 1  	\
| nmrPipe  -fn ZF -auto 					\
| nmrPipe  -fn FT -auto  				\
| nmrPipe  -fn PS -p0 -29.2  -p1 -29.0  -di 			\
| nmrPipe  -fn TP   \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 1	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT  				\
| nmrPipe  -fn PS -p0 1.6 -p1 -3  -di	\
| nmrPipe -fn MED					\
| nmrPipe  -fn TP					 \
| nmrPipe -fn MED 					\
| nmrPipe -out ../ft/$1.ft -ov -verb			\

echo done
