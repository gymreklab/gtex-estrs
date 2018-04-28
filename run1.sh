#!/bin/bash
cat ../650GTEX/population2 | xargs -P10 -I% -n1 ./runfiltervcf.sh %

