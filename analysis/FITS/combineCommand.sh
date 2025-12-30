#!/bin/bash

#./combineCommand.sh --combine-only
#./combineCommand.sh --build-only

#cardDIR=$1
#wsDIR=$2

wsDIR="WS"
cardDIR="workspaces"
resultDir="workspaces"
resultFile="results_DEC16.txt"
bin=""

##########
##########
##########

MODE="all"   # default

while [[ $# -gt 0 ]]; do
    case "$1" in
        --combine-only)
            MODE="combine"
            shift
            ;;
        --build-only)
            MODE="build"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

##########
##########
##########


if [[ "$MODE" == "all" || "$MODE" == "build" ]]; then

rm -rf $cardDIR
mkdir $cardDIR


echo $cardDIR
echo $wsDIR

# Years and categories
#years=("Run3")
#cats=("VLcat" "VHcat" "Zinvcat" "TTLcat" "TTHcat")
#bins=("")

years=("12022" "22022" "12023" "22023" "2024" "2025")
cats=("ggHcat" "VBFcat")
bins=("bdt0" "bdt1" "bdt2")

# Generate datacards for each category & year
for year in "${years[@]}"; do
    for mycat in "${cats[@]}"; do
	for mybin in "${bins[@]}"; do
            echo "Processing: $mycat  $year $mybin"

            python bwsHrare.py \
		   --whichCat="$mycat" \
		   --whichBIN="$mybin" \
		   --inputFileSig="$wsDIR/Signal_${mycat}_${mybin}_${year}_workspace.root" \
		   --inputFileBKG="$wsDIR/Bkg_${mycat}_${mybin}_${year}_workspace.root" \
		   --whichYear="${year}" \
		   --output="$cardDIR/workspace_${mycat}_${mybin}_${year}_workspace.root" \
		   --datCardName="$cardDIR/datacard_${mycat}_${mybin}_${year}.txt"
	done
    done
done

for mycat in "${cats[@]}"; do
    args=()
    found_any=false

    for mybin in "${bins[@]}"; do
        for year in "${years[@]}"; do

            f="$cardDIR/datacard_${mycat}_${mybin}_${year}.txt"

            if [[ -f "$f" ]]; then
                # label uses both bin and year
                label="${mybin}_${year}"
                args+=("${label}=${f}")
                found_any=true
            else
                echo "  WARNING: missing datacard, skipping: $f"
            fi

        done
    done

    if [[ "$found_any" == true ]]; then
        out="$cardDIR/datacard_${mycat}_Run3_combined.txt"

        echo "  Running: combineCards.py ${args[*]} > $out"
        combineCards.py "${args[@]}" > "$out"
        echo "  Created: $out"
    else
        echo "  No datacards found for category $mycat — no combined file created."
        continue
    fi

    # Temporary file name
    tmp="${out}.TMP"

    mv "$out" "$tmp"

    # Fix path duplication
    sed -e 's#workspaces/workspaces/#workspaces/#g' "$tmp" > "$out"

    rm -f "$tmp"

done


fi

if [[ "$MODE" == "all" || "$MODE" == "combine" ]]; then

echo ' **** GFcat ****' > $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_ggHcat_Run3_combined.txt -n _DEC16_ggH --run expected >> $resultFile
mv higgsCombine_DEC16_ggH.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VBFcat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VBFcat_Run3_combined.txt -n _DEC16_VBF --run expected >> $resultFile
mv higgsCombine_DEC16_VBF.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VL cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VLcat_Run3_combined.txt -n _DEC16_VL --run expected >> $resultFile
mv higgsCombine_DEC16_VL.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VH cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VHcat_Run3_combined.txt -n _DEC16_VH --run expected >> $resultFile
mv higgsCombine_DEC16_VH.AsymptoticLimits.mH125.root $cardDIR

echo ' **** Zinv cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Zinvcat_Run3_combined.txt -n _DEC16_Zinv --run expected >> $resultFile
mv higgsCombine_DEC16_Zinv.AsymptoticLimits.mH125.root $cardDIR

echo ' **** TTL cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_TTLcat_Run3_combined.txt -n _DEC16_TTL --run expected >> $resultFile
mv higgsCombine_DEC16_TTL.AsymptoticLimits.mH125.root $cardDIR

echo ' **** TTH cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_TTHcat_Run3_combined.txt -n _DEC16_TTH --run expected >> $resultFile
mv higgsCombine_DEC16_TTH.AsymptoticLimits.mH125.root $cardDIR

fi
