#!/bin/csh

bruk2pipe -in 2/ser \
  -bad 0.0 -aswap -DMX -decim 3029.33333333333 -dspfvs 20 -grpdly 67.9857482910156  \
  -xN              2048  -yN               512  \
  -xT              1024  -yT               512  \
  -xMODE            DQD  -yMODE           Real  \
  -xSW         6602.113  -ySW         6602.839  \
  -xOBS         600.233  -yOBS         600.233  \
  -xCAR           4.752  -yCAR           4.752  \
  -xLAB             1Hx  -yLAB              1H  \
  -ndim               2  -aq2D          Magnitude	\
  -out ./ascr6_cosy.fid -verb -ov

sleep 5
