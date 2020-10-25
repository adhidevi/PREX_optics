#!/bin/bash
run=21438
while read cut
do
	echo $cut
#	root -l -b -q "AsymForAccXcut.C($run,$cut)"
	root -l -b -q "AsymForAccPcut.C($run,$cut)"
#	root -l -b -q "AsymForAccYcutLo.C($run,$cut)"
#	root -l -b -q "AsymForAccYcutHi.C($run,$cut)"
done<TextFiles/pCutR.list
#   done<TextFiles/pCutL.list
