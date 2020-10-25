#!/bin/bash
rm -rf ./temp/*
IFS=$'\t'
#
while read date target runL runR
do
	echo $date $target $runL $runR
#	root -l -b -q plotThetaPhi.C\($runL,\"$target\"\)
#	root -l -b -q plotThetaPhi.C\($runR,\"$target\"\)
#	root -l -b -q theta.C\($runL,\"$target\"\)
#	root -l -b -q theta.C\($runR,\"$target\"\)
#	root -l -b -q Qsq_ana_cal.C\($runL,\"$target\"\)
#	root -l -b -q Qsq_ana_cal.C\($runR,\"$target\"\)
	root -l -b -q "plotMomentum.C($runL, $runR,\"$target\")"
done < run_list3.list
#pdfunite ./temp/*.pdf ./plots/transport_angle_theta_phi_LR.pdf
#pdfunite ./temp/*.pdf ./plots/scattering_angle_LR.pdf
#pdfunite ./temp/*.pdf ./plots/Qsquare_compare_LR.pdf
pdfunite ./temp/*.pdf ./plots/momentum_distribution_LR.pdf
#No need to run the following script
#while read runL runR target
#do
#	echo $runL $runR $target
#	root -l -b -q "plotMomentum.C($runL, $runR,\"$target\")"
#done < run_list2.list
#pdfunite ./temp/*.pdf ./plots/momentum_distribution_LR.pdf
#
rm -rf ./temp/*
