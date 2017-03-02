#!/bin/csh

bruk2pipe -in $1/fid \
  -bad 0.0 -aswap -AMX -decim 2773.33333333333 -dspfvs 20 -grpdly 67.9858856201172  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW         7211.538  \
  -xOBS         600.233  \
  -xCAR           4.800  \
  -xLAB              1H  \
  -ndim               1  \
  -out convfiddata/hs$1.fid -verb -ov \
#
nmrPipe  -in convfiddata/hs$1.fid                                \
| nmrPipe  -fn EM -lb 1.0 -c 0.5        \
| nmrPipe  -fn ZF -auto                                 \
| nmrPipe  -fn FT -auto -verb                           \
| nmrPipe  -fn PS -p0 130.0  -p1 0.0 -di        \
| nmrPipe  -fn POLY -auto               \
| nmrPipe -out procdat/hs$1.ft -ov -verb                 \


sleep 5
