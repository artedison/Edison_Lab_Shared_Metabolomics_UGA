#!/bin/csh

var2pipe -in His_20mM_C.fid/fid \
 -noaswap  \
  -xN             96154  \
  -xT             48077  \
  -xMODE        Complex  \
  -xSW        32051.282  \
  -xOBS         150.805  \
  -xCAR         100.000  \
  -xLAB             C13  \
  -ndim               1  \
  -out ./H_13C_1D.fid -verb -ov

sleep 5
