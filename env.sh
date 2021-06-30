
#export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/bin/:$PATH

export PATH=$PATH:/cvmfs/sft.cern.ch/lcg/releases/lcgenv/latest

cd /cvmfs/sft.cern.ch/lcg/views/LCG_96c_LS/x86_64-centos7-gcc9-opt/
source setup.sh
cd -
#exit

#export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases/
#export PATH=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/latest/:$PATH
#export PATH=/usr/bin:$PATH
#cd /cvmfs/sft.cern.ch/lcg/views/

#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt CMake`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt cmaketools`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt Geant4`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt Qt`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt Qt5`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt fastjet`"
#eval "`lcgenv_py2 -p LCG_96c_LS x86_64-centos7-gcc9-opt pythia8 244`"
#cd -

export PYTHIA8_HOME=$PYTHIA8
export PYTHIA8__HOME=$PYTHIA8_HOME  
export LD_LIBRARY_PATH=$PYTHIA8__HOME/lib:$LD_LIBRARY_PATH

export FASTJET_HOME=/cvmfs/sft.cern.ch/lcg/views/LCG_96c_LS/x86_64-centos7-gcc9-opt/
export FASTJET__HOME=$FASTJET_HOME

export G4RHADRONS=/afs/cern.ch/work/j/jkiesele/public/software/G4Rhadrons_with_Regge_model


export G4DATAINSTALL=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Geant4/10.04/x86_64-slc6-gcc7-opt/share/Geant4-10.4.0
export G4LEVELGAMMADATA=$G4DATAINSTALL/data/PhotonEvaporation5.2
export G4NEUTRONXSDATA=$G4DATAINSTALL/data/G4NEUTRONXS1.4
export G4LEDATA=$G4DATAINSTALL/data/G4EMLOW7.3
export G4SAIDXSDATA=$G4DATAINSTALL/data/G4SAIDDATA1.1
export G4RADIOACTIVEDATA=$G4DATAINSTALL/data/RadioactiveDecay5.2
export G4NEUTRONHPDATA=$G4DATAINSTALL/data/G4NDL4.5
export G4ENSDFSTATEDATA=$G4DATAINSTALL/data/G4ENSDFSTATE2.2
