#!/bin/csh
#
#
  nmrPipe  -in H_inad.fid \
| nmrPipe  -fn SP -off 0.0 -end 0.98 -pow 2 -c .5 	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto	-verb				\
| nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di			\
| nmrPipe  -fn TP 					\
#| nmrPipe  -fn LP -fb					\
| nmrPipe  -fn SP -off 0.0 -end 0.98 -pow 2 		\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto				\
| nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di			\
#| nmrPipe  -fn MED 					\
| nmrPipe  -fn POLY -auto			\
| nmrPipe  -fn TP                                    	\
#| nmrPipe  -fn MED					\
| nmrPipe  -fn POLY -auto			\
| nmrPipe -out H_inad.ft -ov -verb 			\

