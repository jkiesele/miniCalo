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

B4JetGeneratorAction::B4JetGeneratorAction()
 :// B4PrimaryGeneratorAction(),
   pyqq_("/cvmfs/sft.cern.ch/lcg/releases/LCG_latest/MCGenerators/pythia8/244/x86_64-centos7-gcc9-opt/share/Pythia"),
   pygg_("/cvmfs/sft.cern.ch/lcg/releases/LCG_latest/MCGenerators/pythia8/244/x86_64-centos7-gcc9-opt/share/Pythia"),
   jetDef_(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4)),
   fjinputs_()
{
  G4INCL::Random::setGenerator(new G4INCL::Ranecu());

  for (auto* py : {&pygg_, &pyqq_}) {
    py->readString("Beams:eCM = 14000.");
    py->readString("Init:showChangedParticleData = off");
    py->readString("Init:showChangedSettings = on");
    py->readString("Next:numberShowLHA = 0");
    py->readString("Next:numberShowInfo = 0");
    py->readString("Next:numberShowProcess = 0");
    py->readString("Next:numberShowEvent = 0");
    py->readString("PhaseSpace:pTHatMin = 10.");
  }
  pygg_.readString("HardQCD:gg2gg = on");
  pyqq_.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  pygg_.init();
  pyqq_.init();

  fjinputs_.reserve(4096);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4JetGeneratorAction::~B4JetGeneratorAction()
{
  pygg_.stat();
  pyqq_.stat();

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
    pyqq_.rndm.init(G4INCL::Random::getSeeds()[0]);
    pygg_.rndm.init(G4INCL::Random::getSeeds()[1]);
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

  // make a dummy primary proton
  G4PrimaryParticle* primaryParticle = new G4PrimaryParticle(2212);

  Pythia8::Pythia* pythia(&pygg_);//G4INCL::Random::shoot() > 0.5 ? &pyqq_ : &pygg_);
  if (pythia == &pyqq_) {
    eventType_ = kQQ;
    jetType_ = 0;
  }
  else {
    eventType_ = kGG;
    jetType_ = 1;
  }

  unsigned iTry(0);
  while (true) {
    if (!pythia->next()) {
      if (++iTry < 10)
        continue;
      pythia->stat();
      break;
    }

    fjinputs_.clear();

    // Find number of all final charged particles and fill histogram.
    for (int i(0); i < pythia->event.size(); ++i) {
      auto& part(pythia->event[i]);
      if (part.isFinal()) {
        fjinputs_.emplace_back(part.px(), part.py(), part.pz(), part.e());
        fjinputs_.back().set_user_index(i);
      }
    }

    fastjet::ClusterSequence seq(fjinputs_, *jetDef_);
    auto jets(fastjet::sorted_by_pt(seq.inclusive_jets(15.)));

    fastjet::PseudoJet* jet(nullptr);

    for (auto& j : jets) {
      if (j.e() > eMin && j.e() < eMax && fabs(j.eta())<etaTargetMax && fabs(j.eta())>etaTargetMin) {
        jet = &j;
        break;
      }
    }

    if (jet == nullptr)
      continue;

    energy_ = jet->e()*GeV;

    auto constituents(seq.constituents(*jet));

    // bring all particles around +z, +y axis
    if (jet->pz() < 0.) {
      for (auto& pj : constituents)
        pj.reset_momentum(pj.px(), pj.py(), -pj.pz(), pj.E());

      jet->reset_momentum(jet->px(), jet->py(), -jet->pz(), jet->E());
    }


    auto* partDefinition(G4ParticleTable::GetParticleTable()->FindParticle(2212));
    primaryParticle->SetMass(jet->m()*GeV);
    primaryParticle->SetMomentum(jet->px()*GeV, jet->py()*GeV, jet->pz()*GeV);
    primaryParticle->SetCharge(partDefinition->GetPDGCharge());

    if (WriteTruthTree) {
      jetE_ = jet->e();
      jetEta_ = jet->eta();
      nPart_ = 0;
    }

    //
    // Write not just the jet but ALL particles
    //
    //  for (auto& pj : constituents) {
    //   int pdgId(pythia->event[pj.user_index()].id());
    for (int i(0); i < pythia->event.size(); ++i) {
        auto& pj(pythia->event[i]);
        int pdgId= pythia->event[i].id();
        partDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdgId);
        if (partDefinition == nullptr)
            continue; //throw std::runtime_error(std::string("Unknown particle ") + std::to_string(pdgId));

        G4PrimaryParticle* particle = new G4PrimaryParticle(pdgId);
        particle->SetMass(pj.m()*GeV);
        particle->SetMomentum(pj.px()*GeV, pj.py()*GeV, pj.pz()*GeV);
        particle->SetCharge(partDefinition->GetPDGCharge());

        primaryParticle->SetDaughter(particle);

        if (WriteTruthTree) {
            partPid_[nPart_] = pj.id();
            partE_[nPart_] = pj.e();
            partEta_[nPart_] = pj.eta();
            partPhi_[nPart_] = pj.phi();
            ++nPart_;
        }
    }

    if (WriteTruthTree)
      truthTree_->Fill();

    break;
  }

  // create G4PrimaryVertex object
  G4PrimaryVertex* vertex = new G4PrimaryVertex(vertex_position, vertex_time);
  vertex->SetPrimary(primaryParticle);

  anEvent->AddPrimaryVertex(vertex);
}

/*static*/
G4String
B4JetGeneratorAction::eventTypeName(unsigned t)
{
  switch (t) {
  case kQQ:
    return "QQ";
  case kGG:
    return "GG";
  default:
    throw std::runtime_error("Invalid event type");
  };
}

#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
