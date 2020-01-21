export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Python/2.7.13/x86_64-centos7-gcc7-opt/bin/:$PATH
export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases/
export PATH=/cvmfs/sft.cern.ch/lcg/releases/lcgenv/latest/:$PATH
export PATH=/usr/bin:$PATH
eval "`lcgenv -p LCG_93 x86_64-centos7-gcc7-opt Geant4`"
eval "`lcgenv -p LCG_93 x86_64-centos7-gcc7-opt Qt5`"

export G4INSTALL=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Geant4/10.04/x86_64-centos7-gcc7-opt/share/Geant4-10.4.0
export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation5.2
export G4NEUTRONXSDATA=$G4INSTALL/data/G4NEUTRONXS1.4
export G4LEDATA=$G4INSTALL/data/G4EMLOW7.3
export G4SAIDXSDATA=$G4INSTALL/data/G4SAIDDATA1.1
export G4RADIOACTIVEDATA=$G4INSTALL/data/RadioactiveDecay5.2
export G4NEUTRONHPDATA=$G4INSTALL/data/G4NDL4.5
export G4ENSDFSTATEDATA=$G4INSTALL/data/G4ENSDFSTATE2.2
