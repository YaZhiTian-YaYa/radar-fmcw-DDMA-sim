#!/bin/bash
R="F:/Work3/CH_NEW_95_射频_DDMA"
cd /c/Users/NJ/AppData/Local/Temp/probe95
declare -a KN=(enablePhaseNoise enableHopErr pnTrackCoef enablePhaseErr enableTxCoupling enablePaIM3 enableRxCoupling enableRipple enableNonlin enableFeedDisp enableLeakNear enableLeakDC enablePhaseQuant)
for K in "${KN[@]}"; do
  case $K in
    pnTrackCoef) PAT='^pnTrackCoef    = 0;'; REP='pnTrackCoef    = 7.6e-9;'; RPAT='^pnTrackCoef    = 7.6e-9;'; RREP='pnTrackCoef    = 0;' ;;
    *) PAT="^$K *= *false; %IDEAL"; REP="$K = true;"; RPAT="^$K *= *true;" ; RREP="$K = false; %IDEAL" ;;
  esac
  sed -i -E "s/$PAT/$REP/" "$R/rf_sim_core.m"
  if grep -qE "$RPAT" "$R/rf_sim_core.m"; then
    W=$(matlab -batch "run('diagWall2.m')" 2>&1 | iconv -f GBK -t UTF-8 2>/dev/null | grep "排线±4" | grep -o "+[0-9.]*$")
    echo "单独开 $K → 墙高 $W"
    sed -i -E "s/$RPAT/$RREP/" "$R/rf_sim_core.m"
  else
    echo "单独开 $K → sed 未命中"
  fi
done
echo BISECT_ALL_DONE
