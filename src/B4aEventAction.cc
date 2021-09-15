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
// $Id: B4aEventAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4aEventAction.cc
/// \brief Implementation of the B4aEventAction class

#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "B4Analysis.hh"
#include "G4RootAnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "G4INCLRandom.hh"
#include "B4PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   generator_(0),
   detector_(0)
{

    LArEnergy=0;
    rodEnergy=0;
    nevents_=0;
    nsteps_=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}


bool B4aEventAction::checkConstruct(){
    if(!detector_)
        return false;

    return true;
}

void B4aEventAction::accumulateVolumeInfo(G4VPhysicalVolume * volume,const G4Step* step){

	const auto& activesensors=detector_->getActiveSensors();

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4int copyNo = theTouchable->GetVolume()->GetCopyNo();

	size_t idx=activesensors->size();
	for(size_t i=0;i<activesensors->size();i++){
		if(volume == activesensors->at(i).getVol()){
		    //if(copyNo == activesensors->at(i).getCopyNo()){
		        idx=i;
		        break;

		   // }
		}
	}

	if(idx>=activesensors->size())return;//not active volume

	auto & sensor = activesensors->at(idx);

	if(sensor.getLayer() > 0){
	    rodEnergy+= step->GetTotalEnergyDeposit()/1000.; //in GeV
	}
	else{
	    LArEnergy+= step->GetTotalEnergyDeposit()/1000.;//in GeV
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event


  LArEnergy=0;
  rodEnergy=0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

    //check if event is fine

    if(event->IsAborted()){//don't fill information
        G4cout << "event was aborted, not writing output" << G4endl;
        return;
    }

    //fill ntuple

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(0,generator_->getEnergy());
    analysisManager->FillNtupleDColumn(1,LArEnergy);
    analysisManager->FillNtupleDColumn(2,rodEnergy);
    analysisManager->FillNtupleDColumn(3,generator_->x_component_);
    analysisManager->FillNtupleDColumn(4,generator_->y_component_);
    analysisManager->FillNtupleDColumn(5,generator_->z_component_);

    analysisManager->AddNtupleRow();


    nevents_++;

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
