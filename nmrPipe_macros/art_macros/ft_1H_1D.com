#!/bin/csh
#
#
  nmrPipe  -in H_1H_1D.fid \
| nmrPipe  -fn EM -lb 2 \
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto	-verb				\
| nmrPipe  -fn PS -p0 86.0 -p1 -11.0 -di			\
| nmrPipe  -fn POLY -auto					\
| nmrPipe -out H_1H_1D.ft -ov -verb 			\

