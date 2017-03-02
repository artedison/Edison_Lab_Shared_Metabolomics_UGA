#!/bin/csh

bruk2pipe -in ./4/ser \
  -bad 0.0 -aswap -DMX -decim 2218.66666666667 -dspfvs 20 -grpdly 67.9862518310547  \
  -xN              4096  -yN              1024  \
  -xT              2048  -yT               512  \
  -xMODE            DQD  -yMODE        Echo-AntiEcho  \
  -xSW         9014.423  -ySW         9014.423  \
  -xOBS         600.233  -yOBS         600.233  \
  -xCAR           4.772  -yCAR           4.772  \
  -xLAB             1Hx  -yLAB              1H  \
  -ndim               2  -aq2D          States  \
  -out ./4_tocsy.fid -verb -ov

sleep 5
