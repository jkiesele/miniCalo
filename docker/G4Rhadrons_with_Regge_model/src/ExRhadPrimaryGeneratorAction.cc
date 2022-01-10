//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExRhadPrimaryGeneratorAction.cc,v 1.8 2006/07/11 08:26:06 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadPrimaryGeneratorAction.hh"

#include "ExRhadDetectorConstruction.hh"
#include "ExRhadPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CustomPDGParser.h"
#include "Randomize.hh"

using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPrimaryGeneratorAction::ExRhadPrimaryGeneratorAction(
                                             ExRhadDetectorConstruction* ExRhadDC)
  :ExRhadDetector(ExRhadDC),rndmFlag("off"),customGunPos(false),flaten(false)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //create a messenger for this class
  gunMessenger = new ExRhadPrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="pi+");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(10.*GeV);
  G4double position = -1.0*(ExRhadDetector->GetWorldSizeX());
  G4cout<<"World Size is: "<<ExRhadDetector->GetWorldSizeX() / cm<<" cm"<<G4endl;
  particleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPrimaryGeneratorAction::~ExRhadPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double x0 = -1.0*(ExRhadDetector->GetWorldSizeX());
  if (customGunPos) x0=xGun;
  G4double y0 = 0.*cm, z0 = 0.*cm;
  if (rndmFlag == "on")
     {y0 = (ExRhadDetector->GetWorldSizeYZ())*(G4UniformRand()-0.5);
      z0 = (ExRhadDetector->GetWorldSizeYZ())*(G4UniformRand()-0.5);
     } 
  if(flaten) particleGun->SetParticleEnergy(G4UniformRand()*emax);

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
  G4cout<<"Generating energy: "<<particleGun->GetParticleEnergy()/GeV<<" GeV"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPrimaryGeneratorAction::SetFlatEn(G4double val)
{
  if (val == 0){
    flaten=false;
  } else {
    flaten=true;
    emax = val;
  }
}

void ExRhadPrimaryGeneratorAction::SetBeta(G4double val)
{
  
  double Enew = particleGun->GetParticleDefinition()->GetPDGMass()*1/sqrt(1-val*val);
  Enew -= particleGun->GetParticleDefinition()->GetPDGMass();
  G4cout<<"Fisse"<<G4endl;
  particleGun->SetParticleEnergy(Enew);
  return;
}




