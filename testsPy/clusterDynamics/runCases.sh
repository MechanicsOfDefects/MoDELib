#!/bin/bash
# Runs the cases that check the nucleation of immobile clusters and their conversion to discrete loops.
# Each case is a copy of tutorials/spatialCDtest with a few input variables changed.
#
# Usage, from this folder, after building MoDELib with the default species counts (4 mobile, 8 immobile):
#   ./runCases.sh                      runs all the cases
#   ./runCases.sh nuc disc             runs the cases named
#   python3 checkCases.py              compares the results with the expected values
#
# BUILD   build folder (default ../../build)
# Each case takes a few minutes; "reg" takes about ten.

set -e
ROOT=$(cd ../.. && pwd)
BUILD=${BUILD:-$ROOT/build}
MG=$BUILD/tools/MicrostructureGenerator/microstructureGenerator
DD=$BUILD/tools/DDomp/DDomp

setvar () { # file key value. Stops if the key is not in the file
    grep -q "^$2 *=" "$1" || { echo "variable $2 not found in $1"; exit 1; }
    sed -i "s|^$2 *=[^;]*;|$2=$3;|" "$1"
}

discreteLoops () { # adds dislocation dynamics, with climb only, and converts the fields at the first step
    setvar inputFiles/DefectiveCrystal.txt physics "ClusterDynamics DislocationDynamics"
    setvar inputFiles/Zr4_Fitted.txt clusterDiscretizationTime "0.0"
    cp $ROOT/Library/DislocationDynamics/DD.txt inputFiles/DD.txt
    setvar inputFiles/DD.txt glideSolverType none
    setvar inputFiles/DD.txt climbSolverType Galerkin
    setvar inputFiles/DD.txt Lmin 25
    setvar inputFiles/DD.txt Lmax 150
    setvar inputFiles/DD.txt outputQuadraturePoints 0
}

CASES=${@:-reg nuc nuc0 clus clusRef disc discAll discMin}
for CASE in $CASES
do
    rm -rf $CASE
    mkdir $CASE
    cp $ROOT/tutorials/spatialCDtest/generateInputFiles.py $CASE/
    cd $CASE
    sed -i 's|"../../python/"|"'$ROOT'/python/"|; s|\.\./\.\./Library/|'$ROOT'/Library/|g' generateInputFiles.py
    python3 generateInputFiles.py > generateInputFiles.log 2>&1
    MAT=inputFiles/Zr4_Fitted.txt
    DC=inputFiles/DefectiveCrystal.txt
    case $CASE in
        reg)      # the tutorial as it is: nucleation off
            setvar $DC Nsteps 6
            ;;
        nuc)      # clusters born in cascades and by clustering
            setvar $DC Nsteps 2
            setvar $MAT loopCascadeFractions "0.005 0.000704 0.000704 0.000704"
            setvar $MAT mobileSpeciesCascadeFractions "0.995 0.597888 0.3212 0.0788"
            setvar $MAT loopClusteringNucleation 1
            ;;
        nuc0)     # the mobile production of nuc, without nucleation
            setvar $DC Nsteps 2
            setvar $MAT mobileSpeciesCascadeFractions "0.995 0.597888 0.3212 0.0788"
            ;;
        clus)     # clustering alone, with few clusters so that the source is resolved by the six digits of the output
            setvar $DC Nsteps 2
            setvar $MAT loopClusteringNucleation 1
            setvar $MAT initLoopSinks_SI "4e19 4e19 4e19 4e19 2e-7 12e-7 12e-7 12e-7"
            ;;
        clusRef)  # clus without nucleation
            setvar $DC Nsteps 2
            setvar $MAT initLoopSinks_SI "4e19 4e19 4e19 4e19 2e-7 12e-7 12e-7 12e-7"
            ;;
        disc)     # conversion with one loop per 100 clusters, then three climb steps
            setvar $DC Nsteps 4
            setvar $MAT clusterDiscretizationFactor 100
            discreteLoops
            ;;
        discAll)  # conversion with one loop per cluster
            setvar $DC Nsteps 1
            setvar $MAT clusterDiscretizationFactor 1
            discreteLoops
            ;;
        discMin)  # the <c> clusters are smaller than minimumLoopSize and stay in the fields
            setvar $DC Nsteps 3
            setvar $MAT clusterDiscretizationFactor 100
            setvar $MAT minimumLoopSize "30e-9"
            discreteLoops
            ;;
        *)
            echo "unknown case $CASE"; exit 1
            ;;
    esac
    $MG . > microstructureGenerator.log 2>&1
    $DD . > DDomp.log 2>&1
    echo "case $CASE: $(ls evl | grep -c evl_) configurations written"
    cd ..
done
