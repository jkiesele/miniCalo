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

#include "TFile.h"
#include "TTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool WriteTruthTree = false;

B4JetGeneratorAction::B4JetGeneratorAction(particles p) :
// B4PrimaryGeneratorAction(),
  pythia_("/cvmfs/sft.cern.ch/lcg/releases/LCG_latest/MCGenerators/pythia8/244/x86_64-centos7-gcc9-opt/share/Pythia"),
  particle_(p),
  jetDef_(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4)),
  fjinputs_()
{
  G4INCL::Random::setGenerator(new G4INCL::Ranecu());

  pythia_.readString("Beams:eCM = 14000.");
  pythia_.readString("Init:showChangedParticleData = off");
  pythia_.readString("Init:showChangedSettings = on");
  pythia_.readString("Next:numberShowLHA = 0");
  pythia_.readString("Next:numberShowInfo = 0");
  pythia_.readString("Next:numberShowProcess = 0");
  pythia_.readString("Next:numberShowEvent = 0");

  switch (particle_) {
  case minbias:
    pythia_.readString("SoftQCD:nonDiffractive = on");
    pythia_.readString("SoftQCD:singleDiffractive = on");
    pythia_.readString("SoftQCD:doubleDiffractive = on");
    break;
  case displacedjet:
    pythia_.readString("NewGaugeBoson:ffbar2gmZZprime = on");
    pythia_.readString("Zprime:gmZmode = 3"); // only pure Z'
    pythia_.readString("Zprime:ve = 0.");
    pythia_.readString("Zprime:ae = 0.");
    pythia_.readString("Zprime:vnue = 0.");
    pythia_.readString("Zprime:anue = 0.");
    pythia_.readString("32:tau0 = 100.");
    pythia_.readString("32:oneChannel = 1 1. 100 1 -1");
    break;
  default:
    break;
  }
  pythia_.init();

  fjinputs_.reserve(4096);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4JetGeneratorAction::~B4JetGeneratorAction()
{
  pythia_.stat();

  delete jetDef_;

  if (truthTree_ != nullptr) {
    auto* truthFile(truthTree_->GetCurrentFile());
    truthFile->cd();
    truthTree_->Write();
    delete truthFile;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4JetGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (WriteTruthTree && truthTree_ == nullptr) {
    auto* truthFile(TFile::Open("truth.root", "recreate"));
    truthFile->cd();
    truthTree_ = new TTree("truth", "truth tree");

    truthTree_->Branch("jetType", &jetType_, "jetType/I");
    truthTree_->Branch("jetE", &jetE_, "jetE/F");
    truthTree_->Branch("jetEta", &jetEta_, "jetEta/F");
    truthTree_->Branch("nPart", &nPart_, "nPart/i");
    truthTree_->Branch("partPid", partPid_, "partPid[nPart]/I");
    truthTree_->Branch("partE", partE_, "partE[nPart]/F");
    truthTree_->Branch("partEta", partEta_, "partEta[nPart]/F");
    truthTree_->Branch("partPhi", partPhi_, "partPhi[nPart]/F");
  }

  if (firstEvent_) {
    pythia_.rndm.init(G4INCL::Random::getSeeds()[0]);
    firstEvent_ = false;
  }

  double const etaTargetMin = 1.5;
  double const etaTargetMax = 3.0;
  double const eMin = 10.;
  double const eMax = 4000.;

  // jets from this position at eta ~ 3.6 will hit the center of the detector
  xorig_=0;
  yorig_=0;

  G4ThreeVector vertex_position(xorig_,yorig_,0);
  G4double vertex_time(0.);

  // // make a dummy primary proton
  // G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(2212);

  G4double zsign = 1.;
  std::vector<int> primaries;
  primaries.reserve(1024 * 16);
 
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

      G4PrimaryVertex* vertex = new G4PrimaryVertex(vertex_position, vertex_time);
      // let there be light
      G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(22);
      vertex->SetPrimary(primaryParticle);

      anEvent->AddPrimaryVertex(vertex);

      return;
    }

    if (particle_ == minbias) {
      primaries.clear();
      fjinputs_.clear();

      // cluster ak4 jets from final state particles
      for (int i{0}; i < pythia_.event.size(); ++i) {
        auto& part(pythia_.event[i]);

        if (part.isFinal()) {
          if(part.eta()>3.2 || part.eta() < 1.3)continue;
          fjinputs_.emplace_back(part.px(), part.py(), part.pz(), part.e());
          fjinputs_.back().set_user_index(i);

          primaries.push_back(i);
        }
      }
      energy_ = 0;

    }
    else if (particle_ == displacedjet) {
      // ak4 doesn't make sense for displaced jets - simply check if one of the immediate daughters (down quarks) of the Z' points to the detector

      int ipart{0};
      for (; ipart < pythia_.event.size(); ++ipart) {
        auto& part(pythia_.event[ipart]);
        auto& mother(pythia_.event[part.mother1()]); // mother1 is minimum 0 (not -1) so there won't be segfaults
        if (mother.id() != 32)
          continue;

        // I'm lazy so I'm not going to implement an actual filter - just checking that it's "roughly in the detector direction"
        if (part.e() > eMin && part.e() < eMax && std::abs(part.eta()) < etaTargetMax + 0.2 && std::abs(part.eta()) > etaTargetMin - 0.2) {
          if (part.pz() > 0.) {
            vertex_position = G4ThreeVector(mother.xProd(), mother.yProd(), mother.zProd());
          }
          else {
            zsign = -1.;
            vertex_position = G4ThreeVector(mother.xProd(), mother.yProd(), -mother.zProd());
          }
          break;
        }
      }

      if (ipart == pythia_.event.size())
        continue;

      primaries.clear();

      ipart = 0;
      for (; ipart < pythia_.event.size(); ++ipart) {
        auto& part(pythia_.event[ipart]);

        if (!part.isFinal())
          continue;

        auto* mother(&pythia_.event[part.mother1()]);
        while (mother->id() != 32 && mother->mother1() != 0)
          mother = &pythia_.event[mother->mother1()];

        if (mother->id() == 32)
          primaries.push_back(ipart);
      }
    }

    break;
  }

  // auto* partDefinition(G4ParticleTable::GetParticleTable()->FindParticle(2212));
  // primaryParticle->SetMass(jet->m()*GeV);
  // primaryParticle->SetMomentum(jet->px()*GeV, jet->py()*GeV, jet->pz()*GeV);
  // primaryParticle->SetCharge(partDefinition->GetPDGCharge());

  // if (WriteTruthTree) {
  //   jetE_ = jet->e();
  //   jetEta_ = jet->eta();
  //   nPart_ = 0;
  // }

  // create G4PrimaryVertex object
  G4PrimaryVertex* vertex = new G4PrimaryVertex(vertex_position, vertex_time);

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

  if (WriteTruthTree)
    truthTree_->Fill();

  anEvent->AddPrimaryVertex(vertex);
}

#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
