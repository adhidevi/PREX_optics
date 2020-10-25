#!/bin/bash
rm -f ./temp/all_EdgeP_L_run*.pdf
rm -f ./TextFiles/ped_EdgeP_L.csv
while read dateL targetL runL
do
   echo $dateL $targetL $runL
   root -l -b -q "plotmom.C($runL,\"$dateL\",\"$targetL\")"
done<TextFiles/good_runL1.list
rm -f ./temp/all_EdgeP_R_run*.pdf
rm -f ./TextFiles/ped_EdgeP_R.csv
while read dateR targetR runR
do
   echo $dateR $targetR $runR
   root -l -b -q "plotmom.C($runR,\"$dateR\",\"$targetR\")"
done<TextFiles/good_runR1.list
pdfunite ./temp/all_EdgeP_L_*.pdf ./temp/all_EdgeP_R_*.pdf ./plots/ped_EdgeP_all_runs.pdf
rm -f ./temp/all_EdgeP_*
