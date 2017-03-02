#!/bin/csh
#
#
  nmrPipe  -in fid/11_2.fid \
| nmrPipe  -fn SOL -po 1                  \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 2.0 	\
| nmrPipe  -fn ZF -auto 					\
| nmrPipe  -fn FT -auto \
| nmrPipe  -fn PS -p0 -23.0 -p1 0.0 -di			\
| nmrPipe  -fn TP 					\
#| nmrPipe  -fn LP -fb					\
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 		\
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto \
| nmrPipe  -fn PS -p0 97.0 -p1 0.0 -di			\
| nmrPipe  -fn MED 					\
| nmrPipe  -fn TP                                    	\
| nmrPipe  -fn MED					\
| nmrPipe -out ./ft/11_2_HSQC.ft -ov -verb 			\

