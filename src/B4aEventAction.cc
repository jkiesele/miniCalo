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

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   generator_(0)
{
	//create vector ntuple here
//	auto analysisManager = G4AnalysisManager::Instance();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}


void B4aEventAction::accumulateVolumeInfo(G4VPhysicalVolume * volume,const G4Step* step){

	const auto& activesensors=detector_->getActiveSensors();

	bool issensor=true;

	size_t idx=activesensors->size();
	for(size_t i=0;i<activesensors->size();i++){
		if(volume == activesensors->at(i).getVol()){
			idx=i;
			break;
		}
		if(volume == activesensors->at(i).getAbsorberVol()){
			idx=i;
			issensor=false;
			break;
		}
	}

	if(idx==activesensors->size())return;//not active volume
	size_t currentindex=0;

	size_t hitidx=allvolumes_.size();
	hitidx = std::find(allvolumes_.begin(),allvolumes_.end(),activesensors->at(idx).getVol())-allvolumes_.begin();

	if(hitidx != allvolumes_.size()){
		currentindex=hitidx;
	}
	else{
		currentindex=allvolumes_.size();
		allvolumes_.push_back(activesensors->at(idx).getVol());
		rechit_energy_.push_back(0);
		rechit_absorber_energy_.push_back(0);
		rechit_x_.push_back(0);
		rechit_y_.push_back(0);
		rechit_z_.push_back(0);
		rechit_layer_.push_back(0);
		rechit_varea_.push_back(0);
		rechit_vr_.push_back(0);
		rechit_vl_.push_back(0);
	}

	auto energy=step->GetTotalEnergyDeposit();
	if(issensor){
		energy *= activesensors->at(idx).getEnergyscalefactor();
		rechit_energy_.at(currentindex)+=energy;
	}
	else{
		rechit_absorber_energy_.at(currentindex)+=energy;
	}
	rechit_x_.at     (currentindex)=activesensors->at(idx).getPosx();
	rechit_y_.at     (currentindex)=activesensors->at(idx).getPosy();
	rechit_z_.at     (currentindex)=activesensors->at(idx).getPosz();
	rechit_layer_.at (currentindex)=activesensors->at(idx).getLayer();
	rechit_varea_.at (currentindex)=activesensors->at(idx).getArea();
	rechit_vr_.at    (currentindex)=activesensors->at(idx).getDimxy();
	rechit_vl_.at    (currentindex)=activesensors->at(idx).getDimxy();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
  // initialisation per event
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
  clear();

  //set generator stuff
//random particle
  //random energy
  /*
  do this in the generator
  */

  //
  //


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //



  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  
  // fill ntuple
  int i=0;
  int ispart[B4PrimaryGeneratorAction::particles_size];
  for( ;i<B4PrimaryGeneratorAction::particles_size;i++){
	  ispart[i]=B4PrimaryGeneratorAction::globalgen->isParticle(i);
	  analysisManager->FillNtupleIColumn(i,ispart[i]);
  }
  analysisManager->FillNtupleDColumn(i,B4PrimaryGeneratorAction::globalgen->getEnergy());

  //filling deposits and volume info for all volumes automatically..
  for(auto& e:rechit_energy_){
	  if(e<0.01)e=0; //threshold
  }

  analysisManager->AddNtupleRow();  

  clear();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
