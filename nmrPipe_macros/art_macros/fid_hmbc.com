#!/bin/csh

bruk2pipe -in 8/ser \
  -bad 0.0 -aswap -DMX -decim 3029.33333333333 -dspfvs 20 -grpdly 67.9857482910156  \
  -xN              2048  -yN               256  \
  -xT              1024  -yT               128  \
  -xMODE            DQD  -yMODE        Echo-Antiecho\
  -xSW         6602.113  -ySW        36231.884  \
  -xOBS         600.233  -yOBS         150.945  \
  -xCAR           4.752  -yCAR         117.494  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./hmbc_ascr6.fid -verb -ov

sleep 5
