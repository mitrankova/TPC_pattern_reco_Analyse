#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}


source /opt/sphenix/core/bin/sphenix_setup.sh new

#printenv

nEvents="$1"
runnumber=79513
collision="run3pp"
cdbtag="ana532_nocdbtag_v001"
OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/79513_Skip
name="HITS_pp"
#runnumber=75405
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/75405
#runnumber=76905
#collision="run3auau"
#cdbtag="ana514_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/76905
#name="HITS_6x6_1mrad_AuAu"
#runnumber=82626
#collision="run3oo"
#cdbtag="ana537_nocdbtag_v001"
#OutDir=/sphenix/tg/tg01/hf/mitrankova/PatternReco/82626
#name="HITS_OO"
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
echo "  name                = ${name}"
echo ========================================================


 root.exe -l -b << EOF
        .x Fun4All_raw_hit_ZeroField_frog.C(${nEvents},${runnumber}, ${segment},"${OutDir}",${nSkip},"${collision}","${cdbtag}","${name}")
EOF

echo all done