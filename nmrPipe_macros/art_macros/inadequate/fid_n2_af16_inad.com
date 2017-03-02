#!/bin/csh

var2pipe -in N2_AF16_1_inad.fid/fid \
 -noaswap  \
  -xN              8192  -yN              1024  \
  -xT              4096  -yT               512  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW        30487.805  -ySW        30487.805  \
  -xOBS         150.807  -yOBS         150.807  \
  -xCAR           100.000  -yCAR           100.00 \
  -xLAB            C13x  -yLAB            C13y  \
  -ndim               2  -aq2D          States  \
  -out ./n2_af16_inad.fid -verb -ov

sleep 5
