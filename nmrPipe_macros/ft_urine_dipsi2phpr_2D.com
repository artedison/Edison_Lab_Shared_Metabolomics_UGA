#!/bin/csh
#
#

nmrPipe  -in ../fid/$1.fid				\
| nmrPipe  -fn SOL -po 1    \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto -verb				\
| nmrPipe  -fn PS -p0 -64.0  -p1 34.0 -di 			\
| nmrPipe  -fn TP   \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 1.0	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -alt					\
| nmrPipe  -fn PS -p0 -90 -p1 180.0 -di	\
| nmrPipe -fn MED		\
| nmrPipe  -fn TP                                       \
| nmrPipe -fn MED					\
| nmrPipe -out ../ft/$1.ft -ov -verb			\

echo done
