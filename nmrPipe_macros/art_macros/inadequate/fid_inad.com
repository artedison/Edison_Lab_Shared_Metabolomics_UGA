#!/bin/csh

var2pipe -in His_20mM_INAD_128.fid/fid \
 -noaswap  \
  -xN              8192  -yN              1024  \
  -xT              4096  -yT               512  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW        30487.805  -ySW        30487.805  \
  -xOBS         150.807  -yOBS         150.807  \
  -xCAR         109.000  -yCAR         109.000  \
  -xLAB            C13  -yLAB            C13dq  \
  -ndim               2  -aq2D          States  \
  -out ./H_inad.fid -verb -ov

sleep 5
