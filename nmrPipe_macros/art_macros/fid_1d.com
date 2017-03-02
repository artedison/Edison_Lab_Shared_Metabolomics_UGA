#!/bin/csh

bruk2pipe -in 1/fid \
  -bad 0.0 -aswap -DMX -decim 3328 -dspfvs 20 -grpdly 67.9842681884766  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW         6009.615  \
  -xOBS         600.233  \
  -xCAR           4.89  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./ascr6_1D.fid -verb -ov

sleep 5
