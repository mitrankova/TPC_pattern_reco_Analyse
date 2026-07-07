# TPC_pattern_reco_Analyse
All the things for reconstruction

Steps
1. Connecting hits within the modules - primary pattern recogintion.
 Hardware coordinates (pad, layer,timebin)
 IdealPadMap is used to get phi, radius from the hardware coordinates by lookup table from cdb. Now it is idealized, but in future we can upload aligned version (TpcPadMapBuilder must not be used)
2. FullTrackConnector
Connects hardware pieces of inmodule tracks to make a full track in these pseudohardware coordinates - phi, radius, timebin
3. TpcPolyCluster and TpcPolyTrack - both legit, though I'm focusing on TpcPolyCluster as it also provides clusters
The tracks on Polylines - coorect position in x,y,z space!! PHGarfield applied


4. A bunch of vertexers - not used for pattern recogintion itself, but probably some containers assume this information is there...