#!/bin/bash

#TODO) extend the RMin -10,10

#./combineCommand.sh --combine-only
#./combineCommand.sh --build-only

#./combineCommand.sh --combine-only --limits-only
#./combineCommand.sh --combine-only --significance-only

#-M AsymptoticLimits → limits
#-M Significance → sigma

#cardDIR=$1
#wsDIR=$2

day="apr27"

wsDIR="WS_APR27"
cardDIR="workspaces"
resultDir="workspaces"
resultFileLim="results_"$day"_ul.txt"
resultFileSig="results_"$day"_sig.txt"
bin=""


##########
##########
##########

if [[ "$1" == "--help" ]]; then
    echo "Usage:"
    echo "  --build-only           Build datacards only"
    echo "  --combine-only         Run combine only"
    echo "  --limits-only          Run only limits"
    echo "  --significance-only    Run only significance"
    exit 0
fi

MODE="all"   # default
RUN_MODE="both"   # default

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
        --limits-only)
            RUN_MODE="limits"
            shift
            ;;
        --significance-only)
            RUN_MODE="significance"
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

#rm -rf $cardDIR
#mkdir $cardDIR

echo $cardDIR
echo $wsDIR

# Years and categories
years=("Run3")

# Define multiple configurations
cases=(
#    "VHcat,TTHcat,Zinv|bdt0,bdt1,bdt2"   # full
#    "VLcat,TTLcat,VBFcat|bdt0,bdt1"
    "ggHcat,VBFcat,VHcat,TTHcat,Zinvcat,VLcat,TTLcat|incl"
)

for case in "${cases[@]}"; do
    IFS='|' read -r cats_str bins_str <<< "$case"
    IFS=',' read -r -a cats <<< "$cats_str"
    IFS=',' read -r -a bins <<< "$bins_str"

    echo "=== Running case: cats=${cats[*]} bins=${bins[*]} ==="

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

done # close the big for case in "${cases[@]}";

fi


if [[ "$MODE" == "all" || "$MODE" == "combine" ]]; then

echo "" > $resultFileLim
echo "" > $resultFileSig

COMMON_OPTS="--rMin 0 \
--cminDefaultMinimizerStrategy 1 \
--cminDefaultMinimizerTolerance 0.01 \
--X-rtd MINIMIZER_freezeDisassociatedParams \
--X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
#COMMON_OPTS_overkill="--rMin 0 --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.01 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd FITTER_NEWER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd MINIMIZER_MaxCalls=9999999 --X-rtd FAST_VERTICAL_MORPH"

cats=("ggH" "VL" "VH" "Zinv" "TTL" "TTH" "VBF")

for cat in "${cats[@]}"; do

    echo " **** ${cat} cat ****" >> $resultFileLim
    echo " **** ${cat} cat ****" >> $resultFileSig
    card="$cardDIR/datacard_${cat}cat_Run3_combined.txt"
    t2w="$cardDIR/datacard_${cat}cat_Run3_combined.root"

    # -------------------------
    # T2W
    # -------------------------

    text2workspace.py -m 125 $card -o $t2w # this create t2w

    # -------------------------
    # LIMITS
    # -------------------------
    if [[ "$RUN_MODE" == "both" || "$RUN_MODE" == "limits" ]]; then

        echo " **** ${cat} cat LIMIT ****" | tee -a $resultFileLim

        combine -M AsymptoticLimits -m 125 -t -1 \
		--run expected -v 1 \
		$COMMON_OPTS \
		-d $t2w \
		-n _$day\_${cat} >> $resultFileLim

        mv higgsCombine_$day\_${cat}.AsymptoticLimits.mH125.root $cardDIR

    fi

    # -------------------------
    # SIGNIFICANCE
    # -------------------------
    if [[ "$RUN_MODE" == "both" || "$RUN_MODE" == "significance" ]]; then

        echo " **** ${cat} cat SIGNIFICANCE ****" | tee -a $resultFileSig

        combine -M Significance -m 125 -t -1 \
		--expectSignal=1 \
		$COMMON_OPTS \
		-d $t2w \
		-n _$day\_${cat} >> $resultFileSig

    fi

done

# to debug crash  -v 3
# to check model dependency
#--setParameters pdfindex=0
#--setParameters pdfindex=1


fi

exit


echo ' **** GFcat ****' > $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_ggHcat_Run3_combined.txt -n _$day\_ggH --run expected >> $resultFile
mv higgsCombine_$day\_ggH.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VBFcat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VBFcat_Run3_combined.txt -n _$day\_VBF --run expected >> $resultFile
mv higgsCombine_$day\_VBF.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VL cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VLcat_Run3_combined.txt -n _$day\_VL --run expected >> $resultFile
mv higgsCombine_$day\_VL.AsymptoticLimits.mH125.root $cardDIR

echo ' **** VH cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_VHcat_Run3_combined.txt -n _$day\_VH --run expected >> $resultFile
mv higgsCombine_$day\_VH.AsymptoticLimits.mH125.root $cardDIR

echo ' **** Zinv cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Zinvcat_Run3_combined.txt -n _$day\_Zinv --run expected >> $resultFile
mv higgsCombine_$day\_Zinv.AsymptoticLimits.mH125.root $cardDIR

echo ' **** TTL cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_TTLcat_Run3_combined.txt -n _$day\_TTL --run expected >> $resultFile
mv higgsCombine_$day\_TTL.AsymptoticLimits.mH125.root $cardDIR

echo ' **** TTH cat ****' >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_TTHcat_Run3_combined.txt -n _$day\_TTH --run expected >> $resultFile
mv higgsCombine_$day\_TTH.AsymptoticLimits.mH125.root $cardDIR

fi
