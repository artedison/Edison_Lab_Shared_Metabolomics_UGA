#!/bin/csh

set name = "$1_1D.fid"

bruk2pipe -in ../$1/fid \
  -bad 0.0 -aswap -DMX -decim 3029.33333333333 -dspfvs 20 -grpdly 67.9857482910156  \
  -xN             30208  \
  -xT             14986  \
  -xMODE            DQD  \
  -xSW         6602.113  \
  -xOBS         600.233  \
  -xCAR           4.772  \
  -xLAB              1H  \
  -ndim               1  \
  -out ../fid/$name -verb -ov

sleep 1
