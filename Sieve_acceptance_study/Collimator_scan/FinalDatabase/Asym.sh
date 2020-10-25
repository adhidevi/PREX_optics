#!/bin/bash
#rm -f ./TextFiles/finalDB_953MeVsimulation_app_vs_data_sept_scan_Ycoll_scan_run2052.csv
#rm -f ./TextFiles/finalDB_953MeVsimulation_ver_vs_app_sept_scan_Ycoll_scan_run2052.csv
rm -f ./TextFiles/finalDB_953MeVsimulation_app_vs_data_sept_scan_Ycoll_scan_run21185.csv
rm -f ./TextFiles/finalDB_953MeVsimulation_ver_vs_app_sept_scan_Ycoll_scan_run21185.csv
peakpos=0.9495
simMomCut=0.0022
while read run sept YColl_offset ADC th0 dp offset date target Ebeam
do
   echo $run $sept $YColl_offset $ADC $th0 $dp $offset $date $target $Ebeam
#   root -l -b -q "AsymL.C($run,$sept,$ADC,$th0,$peakpos,$simMomCut,$Ebeam,$YColl_offset)"
   root -l -b -q "AsymR.C($run,$sept,$ADC,$th0,$peakpos,$simMomCut,$Ebeam,$YColl_offset)"
#done<TextFiles/run_listL.list
done<TextFiles/run_listR.list
