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
#include "G4INCLRandom.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   generator_(0),
   nsteps_(0)
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

	nsteps_++;

	G4StepPoint* preStepPoint = step->GetPreStepPoint();
	G4TouchableHistory* theTouchable =
	(G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4int copyNo = theTouchable->GetVolume()->GetCopyNo();
//	G4int motherCopyNo
//	= theTouchable->GetVolume(1)->GetCopyNo();

	size_t idx=activesensors->size();
	for(size_t i=0;i<activesensors->size();i++){
		if(volume == activesensors->at(i).getVol()){
		    if(copyNo == activesensors->at(i).getCopyNo()){

		        idx=i;
		        break;

		    }
		}
	}

	if(idx>=activesensors->size())return;//not active volume


	auto energy=step->GetTotalEnergyDeposit();
	if(idx<rechit_energy_.size()){
	    rechit_energy_.at(idx)+=energy/1000.; //GeV
	}


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
  totalen_=0;
  nsteps_=0;
  //set generator stuff
//random particle
  //random energy
  /*
  do this in the generator
  */

  //
  //

  const auto& activesensors=detector_->getActiveSensors();

  rechit_absorber_energy_ = std::vector<double>(activesensors->size(),0);//reset
  rechit_energy_ = std::vector<double>(activesensors->size(),0);//reset
  rechit_x_.resize(activesensors->size(),0);
  rechit_y_.resize(activesensors->size(),0);
  rechit_z_.resize(activesensors->size(),0);
  rechit_layer_.resize(activesensors->size(),0);
  rechit_varea_.resize(activesensors->size(),0);
  rechit_vz_.resize(activesensors->size(),0);
  rechit_vxy_.resize(activesensors->size(),0);
  rechit_detid_.resize(activesensors->size(),0);
  for(size_t i=0;i<activesensors->size();i++){
      rechit_x_.at     (i)=activesensors->at(i).getPosx();//eta
      rechit_y_.at     (i)=activesensors->at(i).getPosy();//phi
      rechit_z_.at     (i)=activesensors->at(i).getPosz();//z
      rechit_layer_.at (i)=activesensors->at(i).getLayer();
      rechit_varea_.at (i)=activesensors->at(i).getArea();//0
      rechit_vz_.at    (i)=activesensors->at(i).getDimz();//lengthz
      rechit_vxy_.at   (i)=activesensors->at(i).getDimxy();//0
      rechit_detid_.at (i)=activesensors->at(i).getGlobalDetID();
  }


}

double getTrackMomentum(double pt, bool isgamma){
    double smearing=(pt/100.)*(pt/100.)*0.04 +0.01;
    if(!isgamma)
        return pt + pt*smearing*G4INCL::Random::gauss();
    else
        return 0;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::EndOfEventAction(const G4Event* event)
{
  // Accumulate statistics
  //

    G4cout << "nsteps_ "<<nsteps_ <<G4endl;

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
  analysisManager->FillNtupleDColumn(i+1,B4PrimaryGeneratorAction::globalgen->getX());
  analysisManager->FillNtupleDColumn(i+2,B4PrimaryGeneratorAction::globalgen->getY());
  analysisManager->FillNtupleDColumn(i+3,0); //maybe this could be displacement?
  analysisManager->FillNtupleDColumn(i+4,totalen_);


  //filling deposits and volume info for all volumes automatically..


  analysisManager->AddNtupleRow();  

  clear();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
