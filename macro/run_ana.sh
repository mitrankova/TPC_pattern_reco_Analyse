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
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_keff_4
#InDir=/sphenix/user/mitrankov/patterns/run/output_DST/keff_4
#runnumber=75405
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405
#InDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405
runnumber=76905
collision="run3auau"
cdbtag="ana514_nocdbtag_v001"
OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_19t0
InDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905_19t0/output_DST
segment="$2"
nSkip="$3"


echo ========================================================
echo Running Fun4All_TPC_Tracking.C  
echo "  nEvents             = ${nEvents}"
echo "  nSkip               = ${nSkip}"
echo "  Run Number          = ${runnumber}"
echo "  Segment             = ${segment}"
echo "  OutDir              = ${OutDir}"
echo "  InDir               = ${InDir}"
echo "  collision           = ${collision}"
echo "  cdbtag              = ${cdbtag}"
echo ========================================================


 root.exe -l -b << EOF
        .x Fun4All_TPC_Tracking.C(${nEvents},${runnumber}, ${segment},"${OutDir}","${InDir}",${nSkip},"${collision}","${cdbtag}")
EOF

echo all done