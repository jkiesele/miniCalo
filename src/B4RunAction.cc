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
// $Id: B4RunAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B4PrimaryGeneratorAction.hh"

#include "B4aEventAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction(B4PrimaryGeneratorAction *gen, B4aEventAction* ev, G4String fname)
 : G4UserRunAction()
{ 
	fname_=fname;
	eventact_=ev;
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  

  // Creating ntuple
  //
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  generator_=gen;
  G4cout << "creating particle entries" << G4endl;
  auto parts=generator_->generateAvailableParticles();
  for(const auto& p:parts){
	  analysisManager->CreateNtupleIColumn(p);
  }
  analysisManager->CreateNtupleDColumn("true_energy");
  analysisManager->CreateNtupleDColumn("true_x");
  analysisManager->CreateNtupleDColumn("true_y");
  analysisManager->CreateNtupleDColumn("track_momentum");
  analysisManager->CreateNtupleDColumn("total_dep_energy");

//if(false){
  analysisManager->CreateNtupleDColumn("rechit_energy",eventact_->rechit_energy_);
 // analysisManager->CreateNtupleDColumn("rechit_absorber_energy",eventact_->rechit_absorber_energy_);
  analysisManager->CreateNtupleDColumn("rechit_x",eventact_->rechit_x_);
  analysisManager->CreateNtupleDColumn("rechit_y",eventact_->rechit_y_);
  analysisManager->CreateNtupleDColumn("rechit_z",eventact_->rechit_z_);
  analysisManager->CreateNtupleDColumn("rechit_layer",eventact_->rechit_layer_);
  analysisManager->CreateNtupleDColumn("rechit_varea",eventact_->rechit_varea_);
  analysisManager->CreateNtupleDColumn("rechit_vz",eventact_->rechit_vz_);
  analysisManager->CreateNtupleDColumn("rechit_vxy",eventact_->rechit_vxy_);
  analysisManager->CreateNtupleIColumn("rechit_detid",eventact_->rechit_detid_);
//}
  analysisManager->FinishNtuple();

  G4cout << "run action initialised" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = fname_;
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if (false && analysisManager->GetH1(1) ) {
    
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
