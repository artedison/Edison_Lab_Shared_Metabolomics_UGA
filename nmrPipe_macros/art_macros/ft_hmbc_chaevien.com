#!/bin/csh
#
#
  nmrPipe  -in hmbc_ascr6.fid \
| nmrPipe  -fn SP -off 0.0 -end 0.98 -pow 2 -c .5 	\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto	-verb				\
| nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di			\
| nmrPipe  -fn TP 					\
| nmrPipe  -fn LP -fb					\
| nmrPipe  -fn SP -off 0.0 -end 0.98 -pow 2 		\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto					\
| nmrPipe  -fn MC					\
| nmrPipe  -fn MED 					\
| nmrPipe  -fn TP                                    	\
| nmrPipe  -fn MED					\
| nmrPipe -out ascr6_hmbc.ft -ov -verb 			\

