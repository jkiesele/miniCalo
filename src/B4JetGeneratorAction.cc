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
// $Id: B4JetGeneratorAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file B4JetGeneratorAction.cc
/// \brief Implementation of the B4JetGeneratorAction class

#include "B4JetGeneratorAction.hh"

#ifndef NOPYTHIA

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
#include <G4INCLRanecu.hh>
#include<ctime>
#include<sys/types.h>
#include <cstdlib>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool WriteTruthTree = false;

B4JetGeneratorAction::B4JetGeneratorAction(particles p) :
// B4PrimaryGeneratorAction(),
  pythia_("/Applications/pythia8/share/Pythia8/xmldoc"),
  particle_(p),
nPart_(1),
energy_(0)
{

    G4INCL::Random::setGenerator(new G4INCL::Ranecu());

    for(int i=0;i<seedsoffset_;i++)
        G4double rand =  G4INCL::Random::shoot();

    pythia_.readString("Beams:eCM = 8000.");
    pythia_.readString("HardQCD:all = on");
    pythia_.readString("PhaseSpace:pTHatMin = 20.");

  pythia_.init();
  pythia_.rndm.init(seedsoffset_);


  nPU_=1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4JetGeneratorAction::~B4JetGeneratorAction()
{
  pythia_.stat();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B4JetGeneratorAction::GeneratePrimaries(G4Event* anEvent){
    G4ThreeVector vertex_position(xorig_,yorig_,0);
    G4double vertex_time(0.);

    G4PrimaryVertex* vertex = new G4PrimaryVertex(vertex_position, vertex_time);


    GenerateSingleVertex(vertex);

    G4cout << "npart: " << vertex->GetNumberOfParticle() << G4endl;
    anEvent->AddPrimaryVertex(vertex);

}

void B4JetGeneratorAction::GenerateSingleVertex(G4PrimaryVertex* vertex)
{


  if (firstEvent_) {
    firstEvent_ = false;
  }

  double const etaTargetMin = -4.0;
  double const etaTargetMax = 4.0;
  double const eMin = 2.;
  double const eMax = 20.;

  // jets from this position at eta ~ 3.6 will hit the center of the detector
  xorig_=0;
  yorig_=0;
  G4double zorig=430;


  // // make a dummy primary proton
  // G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(2212);

  G4double zsign = 1.;
  std::vector<int> primaries;
  primaries.reserve(1024 * 16);

  G4ThreeVector vertex_position(xorig_,yorig_,zorig);
  G4double vertex_time(0.);
 
  while (true) {
      // loop until an event passing a gen filter is generated

      unsigned iTry(0);
      while (true) {
          // loop until pythia generates a consistent event (should succeed)
          if (!pythia_.next()) {
              if (++iTry < 10)
                  continue;
              pythia_.stat();
          }
          break;
      }

      if (iTry == 10) {
          // failed event generation - make a dummy vertex (or just crash?)
          G4cerr << "EVENT GENERATION FAILED!!!!" << G4endl;

          //G4PrimaryVertex* vertex = new G4PrimaryVertex(vertex_position, vertex_time);
          // let there be light
          G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(22);
          vertex->SetPrimary(primaryParticle);


          return;
      }


      primaries.clear();
      double totalen=0;

      // cluster ak4 jets from final state particles
      for (int i{0}; i < pythia_.event.size(); ++i) {
          auto& part(pythia_.event[i]);

          if (part.isFinal()) {
              // if(part.eta()<4. && part.eta() > 1.0){
              //

              primaries.push_back(i);
              totalen+=part.e();
              //  }
          }
      }
      energy_ = 0;
      G4cout << "total calo energy " << totalen << std::endl;
      break;
  }

  // create G4PrimaryVertex object

  for (int ipart : primaries) {
    auto& pj(pythia_.event[ipart]);

    int pdgId = pj.id();
    auto* partDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdgId);
    if (partDefinition == nullptr)
      continue; //throw std::runtime_error(std::string("Unknown particle ") + std::to_string(pdgId));

    G4PrimaryParticle* particle = new G4PrimaryParticle(pdgId);
    particle->SetMass(pj.m()*GeV);
    particle->SetMomentum(pj.px()*GeV, pj.py()*GeV, pj.pz()*GeV*zsign);
    particle->SetCharge(partDefinition->GetPDGCharge());

    vertex->SetPrimary(particle);

    //      primaryParticle->SetDaughter(particle);

    if (WriteTruthTree) {
      partPid_[nPart_] = pdgId;
      partE_[nPart_] = pj.e();
      partEta_[nPart_] = pj.eta();
      partPhi_[nPart_] = pj.phi();
      ++nPart_;
    }
  }


}

#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
