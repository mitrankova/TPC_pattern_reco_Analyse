#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}


source /opt/sphenix/core/bin/sphenix_setup.sh new

#printenv

nEvents="$1"
runnumber=79513
segment="$2"
nSkip="$3"

OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513

echo ========================================================
echo Running Fun4All_raw_hit_ZeroField_frog.C  p+p2025Zero field
echo "  nEvents             = ${nEvents}"
echo "  nSkip               = ${nSkip}"
echo "  Run Number          = ${runnumber}"
echo "  Segment             = ${segment}"
echo "  OutDir              = ${OutDir}"
echo ========================================================


 root.exe -l -b << EOF
        .x Fun4All_raw_hit_ZeroField_frog.C(${nEvents},${runnumber}, ${segment},"${OutDir}",${nSkip})
EOF

echo all done