#!/bin/csh
#
#
  nmrPipe  -in H_13C_1D.fid \
| nmrPipe  -fn EM -lb 4 \
| nmrPipe  -fn ZF -auto					\
| nmrPipe  -fn FT -auto	-verb				\
| nmrPipe  -fn PS -p0 -47.0 -p1 -68.0 -di			\
| nmrPipe  -fn POLY -auto					\
| nmrPipe -out H_13C_1D.ft -ov -verb 			\

