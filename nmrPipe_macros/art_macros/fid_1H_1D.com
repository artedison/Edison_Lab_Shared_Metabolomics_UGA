#!/bin/csh

var2pipe -in His_20mM_H.fid/fid \
 -noaswap  \
  -xN             32768  \
  -xT             16384  \
  -xMODE        Complex  \
  -xSW         7183.908  \
  -xOBS         599.690  \
  -xCAR           4.773  \
  -xLAB              H1  \
  -ndim               1  \
  -out ./H_1H_1D.fid -verb -ov

sleep 5
