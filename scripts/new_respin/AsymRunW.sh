#!/bin/bash
rm -f TextFiles/Qsq_edgeCut_W_allRunsL.csv
rm -f TextFiles/Qsq_edgeCut_W_allRunsR.csv
while read date target runL runR Ebeam
do
   echo $date $target $runL $runR $Ebeam
   rm -f TextFiles/Qsq_pCut_W_run$runL.csv
   rm -f TextFiles/Qsq_pCut_W_run$runR.csv
   root -l -b -q "Qsq_only_pCut_weighted.C($runL,$Ebeam)"
   root -l -b -q "Qsq_only_pCut_weighted.C($runR,$Ebeam)"
done<TextFiles/good_run.list
pdfunite ./plots/Qsq_pCut_W_run*.pdf ./plots/all_Qsq_W_pCut.pdf
rm -f ./plots/Qsq_pCut_W_*
