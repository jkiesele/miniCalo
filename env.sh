export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-centos7-gcc8-opt/bin/:$PATH

source /cvmfs/sft.cern.ch/lcg/views/LCG_96c_LS/x86_64-centos7-gcc9-opt/setup.sh

export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases/
export PATH=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/latest/:$PATH
#export PATH=/usr/bin:$PATH
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt CMake`"
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt cmaketools`"
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt Geant4`"
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt Qt5`"
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt fastjet`"
eval "`lcgenv -p LCG_96c_LS x86_64-centos7-gcc9-opt pythia8 244`"


export LD_LIBRARY_PATH=$PYTHIA8__HOME/lib:$LD_LIBRARY_PATH

export G4DATAINSTALL=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Geant4/10.04/x86_64-slc6-gcc7-opt/share/Geant4-10.4.0
export G4LEVELGAMMADATA=$G4DATAINSTALL/data/PhotonEvaporation5.2
export G4NEUTRONXSDATA=$G4DATAINSTALL/data/G4NEUTRONXS1.4
export G4LEDATA=$G4DATAINSTALL/data/G4EMLOW7.3
export G4SAIDXSDATA=$G4DATAINSTALL/data/G4SAIDDATA1.1
export G4RADIOACTIVEDATA=$G4DATAINSTALL/data/RadioactiveDecay5.2
export G4NEUTRONHPDATA=$G4DATAINSTALL/data/G4NDL4.5
export G4ENSDFSTATEDATA=$G4DATAINSTALL/data/G4ENSDFSTATE2.2