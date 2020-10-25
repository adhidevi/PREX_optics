#!/bin/bash
#first run the CheckACCshape.C script for a given shape for all th_shift then run this script
rm -f acceptance_gaus_shiftTH.csv
#rm -f acceptance_box_shiftTH.csv
#rm -f acceptance_shift_shiftTH.csv
IFS=$'\t'
while read shape shift
do
	echo $shape $shift
        root -l -b -q "CheckVertexAsym.C(\"$shape\",$shift)"
done < ./TextFiles/gaus_th_shift.list
#done < ./TextFiles/box_th_shift.list
#done < ./TextFiles/original_th_shift.list
pdfunite ./temp1/*gaus_0.0* ./temp1/*gaus_0.1* ./temp1/*gaus_-0.1* ./temp1/*gaus_0.2* ./temp1/*gaus_-0.2* ./temp1/*gaus_0.3* ./temp1/*gaus_-0.3* ./temp1/*gaus_0.4* ./temp1/*gaus_-0.4* ./temp1/*gaus_0.5* ./temp1/*gaus_-0.5* ./plots/acceptance_gaus_shiftTH_953MeV_4by6.pdf
#pdfunite ./temp1/*box_0.0* ./temp1/*box_0.1* ./temp1/*box_-0.1* ./temp1/*box_0.2* ./temp1/*box_-0.2* ./temp1/*box_0.3* ./temp1/*box_-0.3* ./temp1/*box_0.4* ./temp1/*box_-0.4* ./temp1/*box_0.5* ./temp1/*box_-0.5* ./plots/acceptance_box_shiftTH_953MeV_4by6.pdf
#pdfunite ./temp1/*shift_0.0* ./temp1/*shift_0.1* ./temp1/*shift_-0.1* ./temp1/*shift_0.2* ./temp1/*shift_-0.2* ./temp1/*shift_0.3* ./temp1/*shift_-0.3* ./temp1/*shift_0.4* ./temp1/*shift_-0.4* ./temp1/*shift_0.5* ./temp1/*shift_-0.5* ./plots/acceptance_original_shiftTH_953MeV_4by6.pdf
cp ./temp1/* ./temp/
rm -rf ./temp1/*
