#!/bin/bash

DIRECTORY="Results/Run$1/"
DIREC="Results/Run$1"

if [ "$#" != "1" ]; then
    echo ""
    echo "Need run number"
    echo "Exiting..."
    echo ""
    exit
fi

if [ -d "$DIRECTORY" ]; then
    echo ""
    echo "Run number already exists. Please choose another."
    echo "Exiting..."
    exit 
fi

if [ ! -d "Results/" ]; then
    mkdir Results/
fi

if [ ! -d "Results/Plots/" ]; then
    mkdir Results/Plots
fi

if [ ! -f "Program/Simulation" ]; then
    cd Program/
    make
    rm *.o *.gch
    cd ..
fi

mkdir $DIRECTORY
./Program/Simulation $1
mv Run$1_Parameters.dat $DIRECTORY
mv Simulation.dat $DIRECTORY
mv Solenoid.dat $DIRECTORY
mv Helmholtz1.dat $DIRECTORY
mv Helmholtz2.dat $DIRECTORY
sh ./Program/B_Magnitude.gscript $DIRECTORY $1
sh ./Program/VecFields.gscript $DIRECTORY $1
cp $DIREC/Run$1_Solenoid-magnitude.ps Results/Plots/
cp $DIREC/Run$1_Solenoid-Vec.ps Results/Plots/
cp $DIREC/Run$1_Helmholtz1-Vec.ps Results/Plots/
cp $DIREC/Run$1_Helmholtz2-Vec.ps Results/Plots/