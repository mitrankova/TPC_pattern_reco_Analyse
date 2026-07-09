#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}


source /opt/sphenix/core/bin/sphenix_setup.sh new

#printenv

nEvents="$1"
#runnumber=79513
#collision="run3pp"
#cdbtag="ana532_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513
runnumber=75405
collision="run3auau"
cdbtag="ana514_nocdbtag_v001"
OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405
#runnumber=76905
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_19t0
segment="$2"
nSkip="$3"





echo ========================================================
echo Running Fun4All_raw_hit_ZeroField_frog.C
echo "  nEvents             = ${nEvents}"
echo "  nSkip               = ${nSkip}"
echo "  Run Number          = ${runnumber}"
echo "  Segment             = ${segment}"
echo "  OutDir              = ${OutDir}"
echo "  collision           = ${collision}"
echo "  cdbtag              = ${cdbtag}"
echo ========================================================


 root.exe -l -b << EOF
        .x Fun4All_raw_hit_ZeroField_frog.C(${nEvents},${runnumber}, ${segment},"${OutDir}",${nSkip},"${collision}","${cdbtag}")
EOF

echo all done