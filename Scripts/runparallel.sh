cat chromosomes | xargs -P3 -I% -n1 ./get_strs.sh %
