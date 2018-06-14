#!/bin/bash
cat ../650GTEX/population1 | xargs -P10 -I% -n1 ./runfiltervcf.sh %

