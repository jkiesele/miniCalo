#!/bin/bash

number=$1

sourcedir="/afs/cern.ch/user/j/jkiesele/work/HGCal/DisplacedCalo/build/" #end with "/"
eosdir="/eos/home-j/jkiesele/DeepNtuples/HGCalToys/" #end with "/"

THISDIR=`pwd`
cd $sourcedir
source env.sh
cd $THISDIR

cp $sourcedir/exampleB4a .
cp $sourcedir/batchrun.mac .

./exampleB4a -f $number -m batchrun.mac 2>&1 > "${number}.txt"

eoscp "out${number}.root" "${eosdir}minbias_${number}.root"
eoscp "${number}.txt" "${eosdir}minbias_${number}.txt"
rm -f "out${number}.root" "${number}.txt"
