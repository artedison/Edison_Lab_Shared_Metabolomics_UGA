#!/bin/csh

set name = "$1_jres.fid"

bruk2pipe -in ../$1/ser \
  -bad 0.0 -aswap -DMX -decim 2005.33333333333 -dspfvs 20 -grpdly 67.9864501953125  \
  -xN              8192  -yN                64  \
  -xT              4096  -yT                32  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW         9973.404  -ySW           78.125  \
  -xOBS         600.233  -yOBS         600.233  \
  -xCAR           4.754  -yCAR           4.754  \
  -xLAB             1Hx  -yLAB              1H  \
  -ndim               2  -aq2D          States  \
  -out ../fid/$name -verb -ov

sleep 1
