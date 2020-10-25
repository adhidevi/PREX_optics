#!/bin/bash
rm -f ./TextFiles/septum_scan_Ycoll_scan_noefactsim_run2052.csv
rm -f ./TextFiles/vertex_septum_scan_Ycoll_scan_noefactsim_run2052.csv
peakpos=0.9514
simMomCut=0.0022
while read run sept YColl_offset ADC th0 dp offset date target Ebeam
do
   echo $run $sept $YColl_offset $ADC $th0 $dp $offset $date $target $Ebeam
   root -l -b -q "AsymL.C($run,$sept,$ADC,$th0,$peakpos,$simMomCut,$Ebeam,$YColl_offset)"
done<TextFiles/run_listL.list
