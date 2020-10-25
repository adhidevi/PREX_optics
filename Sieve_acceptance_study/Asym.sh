#!/bin/bash
rm -f ./TextFiles/septum_scan_Ycoll_scan_oldDatabase_run2052.csv
#peakpos=0.9512625
peakpos=0.9492875
simMomCut=0.0022
while read run sept YColl_offset ADC th0 dp offset date target Ebeam
do
   echo $run $sept $YColl_offset $ADC $th0 $dp $offset $date $target $Ebeam
#   root -l -b -q "AsymL_outerHole_test.C($run,$sept,$ADC,$th0,$dp,$offset)"
#   root -l -b -q "kAsymL.C($run,$sept,$ADC,$th0,$peakpos,$simMomCut,$Ebeam,$YColl_offset)"
   root -l -b -q "AsymL.C($run,$sept,$ADC,$th0,$peakpos,$simMomCut,$Ebeam,$YColl_offset)"
done<TextFiles/run_listL.list
