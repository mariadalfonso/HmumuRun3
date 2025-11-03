#!/bin/bash

#cardDIR=$1
#wsDIR=$2

wsDIR="WS"
cardDIR="workspaces"
resultDir="workspaces"
bin=""

rm -rf $dataCardDIR

rm -rf $cardDIR
mkdir $cardDIR

##########

echo $cardDIR
echo $wsDIR

#for mycat in "ggHcat" "VBFcat" "Zinvcat" "VHcat" "VLcat" "TTHcat" "TTLcat"; 
for mycat in "ggHcat";
	     
do
    echo $mycat

    python bwsHrare.py --whichCat=$mycat --inputFileSig=$wsDIR/Signal_$mycat\_2024_workspace.root --inputFileBKG=$wsDIR/Bkg_$mycat\_2024_workspace.root --output=$cardDIR/workspace_$mycat\_workspace.root --datCardName=$cardDIR/datacard_$mycat\_workspace.txt

done
