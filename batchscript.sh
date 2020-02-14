#!/bin/bash

number=$1

sourcedir="<afs.... minicalo build>/" #end with "/"
cp $sourcedir/exampleB4a .
cp $sourcedir/batchrun.mac .

./exampleB4a -f $number -m batchrun.mac

eoscp "out${number}.root" "/eos/home-j/jalimena/DisplacedCalo/minbias_${number}.root"
rm -f "out${number}.root"