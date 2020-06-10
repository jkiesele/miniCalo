#!/bin/bash

number=$1

sourcedir="/afs/cern.ch/work/j/jkiesele/muonCalo/miniCalo/batch"
eosdir="/eos/cms/store/cmst3/group/dehep/miniCalo/muoncalo/"

THISDIR=`pwd`
cd $sourcedir
source env.sh
cd $THISDIR

cp $sourcedir/exampleB4a .
cp $sourcedir/batchrun.mac .

./exampleB4a -f $number -m batchrun.mac 2>&1 > "${number}.txt"

eoscp `ls *.root` "${eosdir}${number}.root"

rm -f *.root "${number}.txt" 
