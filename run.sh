# !/bin/bash

PROCESSING_TIMES="data/ta100x20/processing_times_1.txt"
MACHINES="data/ta100x20/machines_1.txt"
MIN_TEMP=0.1
RHO=1
NSWEEP=100
TIME=600
C=0.9
ETA=0.69

./jssp_solver \
  --processing-times "$PROCESSING_TIMES" \
  --machines "$MACHINES" \
  --min-temp "$MIN_TEMP" \
  --rho "$RHO" \
  --nsweep "$NSWEEP" \
  --time "$TIME" \
  --C "$C" \
  --eta "$ETA"