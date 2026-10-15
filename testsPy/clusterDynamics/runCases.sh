#!/bin/bash
# Runs the cases that check the kinetics of the immobile clusters (CVODE).
# Each case is a copy of tutorials/spatialCDcvode with a few input variables changed. They need a build with SUNDIALS.
#
# Usage, from this folder, after building MoDELib with the default species counts (4 mobile, 8 immobile):
#   ./runCases.sh                      runs all the cases
#   ./runCases.sh cvode                runs the cases named
#   python3 checkCases.py              compares the results with the expected values
#
# BUILD   build folder (default ../../build)

set -e
ROOT=$(cd ../.. && pwd)
BUILD=${BUILD:-$ROOT/build}
MG=$BUILD/tools/MicrostructureGenerator/microstructureGenerator
DD=$BUILD/tools/DDomp/DDomp

setvar () { # file key value. Stops if the key is not in the file
    grep -q "^$2 *=" "$1" || { echo "variable $2 not found in $1"; exit 1; }
    sed -i "s|^$2 *=[^;]*;|$2=$3;|" "$1"
}

CASES=${@:-cvode cvodeRestart cvodeStress}
for CASE in $CASES
do
    rm -rf $CASE
    mkdir $CASE
    cp $ROOT/tutorials/spatialCDcvode/generateInputFiles.py $CASE/
    cd $CASE
    sed -i 's|"../../python/"|"'$ROOT'/python/"|; s|\.\./\.\./Library/|'$ROOT'/Library/|g' generateInputFiles.py
    python3 generateInputFiles.py > generateInputFiles.log 2>&1
    FIRSTCONFIGURATION=
    DC=inputFiles/DefectiveCrystal.txt
    case $CASE in
        cvode)    # the D1/M1 model: four steps of 0.1 dpa
            setvar $DC Nsteps 4
            ;;
        cvodeRestart) # a new study that starts from the third configuration of the case cvode
            setvar $DC Nsteps 2
            FIRSTCONFIGURATION=../cvode/evl/evl_2.txt
            ;;
        cvodeStress)  # a tensile stress of 2 MPa along a1 orients the nucleation
            setvar $DC Nsteps 2
            setvar $DC physics "ClusterDynamics ElasticDeformation"
            cp $ROOT/Library/ElasticDeformation/ElasticDeformation.txt inputFiles/ElasticDeformation.txt
            setvar inputFiles/ElasticDeformation.txt ExternalStress0 " 6.0606e-5 0.0 0.0 0.0 0.0 0.0"
            ;;
        *)
            echo "unknown case $CASE"; exit 1
            ;;
    esac
    $MG . > microstructureGenerator.log 2>&1
    if [ -n "$FIRSTCONFIGURATION" ]; then cp $FIRSTCONFIGURATION evl/evl_0.txt; fi
    $DD . > DDomp.log 2>&1
    echo "case $CASE: $(ls evl | grep -c evl_) configurations written"
    cd ..
done
