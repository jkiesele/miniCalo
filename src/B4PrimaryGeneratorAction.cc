//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4PrimaryGeneratorAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4INCLRandom.hh"
#include <G4INCLGeant4Random.hh>
#include <G4INCLRandomSeedVector.hh>
#include<ctime>
#include<sys/types.h>
#include <fstream>

template<class T>
static G4String createString(const T& i){
    std::stringstream ss;
    ss << i;
    std::string number=ss.str();
    return number;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : B4PartGeneratorBase(),
   particleid_(particles_size),
   fParticleGun(nullptr)
{

    //read particles from file
    std::ifstream file("particles.txt");
    std::string a,b,c;//dummy
    int pid;
    if(!file)
        throw std::runtime_error("file particles.txt could not be opened");

    while (file >> pid >> a >> b >> c){
        G4cout << pid <<" " << c << G4endl;
        availParticles_.push_back({c,pid});
    }



    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    // default particle kinematic
    //
    auto particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("Neutralino");
  //  G4ParticleTable::GetParticleTable()->DumpTable();
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(100.*GeV);


    G4INCL::Random::setGenerator( new G4INCL::Geant4RandomGenerator());

    for(int i=0;i<seedsoffset_;i++)
        G4double rand =  G4INCL::Random::shoot();


    from_beamspot_=false;

#ifdef FROMBEAMSPOT
    from_beamspot_=true;
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}
std::vector<G4String> B4PrimaryGeneratorAction::generateAvailableParticles()const{
	std::vector<G4String> out;
	for(const auto& p:availParticles_){
	    auto s = p.first;
	    if(s.at(0) == '-')
	        s[0]='m';
		out.push_back(p.first);
	}
	return out;
}


G4String B4PrimaryGeneratorAction::getParticleName(enum particles p)const{

    return "isInvalid"+createString((int)p);
}


void B4PrimaryGeneratorAction::setParticle(const G4String & part){
    auto particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle(part);
    fParticleGun->SetParticleDefinition(particleDefinition);
    primaryPart_={part,particleDefinition->GetPDGEncoding()};
}

G4String B4PrimaryGeneratorAction::setParticleID(enum particles p){
    return getParticleName(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
 // G4cout << "shooting.. ";
  // Set gun position


  setParticle("~g_rho0");



  G4double energy_max=10;//GeV
  G4double energy_min=2.;

  energy_=10001;
  while(energy_>energy_max){//somehow sometimes the random gen shoots >1??
      G4double rand =  G4INCL::Random::shoot();
      energy_=(energy_max-energy_min)*rand+energy_min;
  }


  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  fParticleGun->SetParticleEnergy(energy_ * GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
  fParticleGun->GeneratePrimaryVertex(anEvent);


 // G4cout << "energy: " << energy_ <<  G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

