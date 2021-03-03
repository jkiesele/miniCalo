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
#include "B4PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::B4aEventAction()
 : G4UserEventAction(),
   fEnergyAbs(0.),
   fEnergyGap(0.),
   fTrackLAbs(0.),
   fTrackLGap(0.),
   generator_(0),
   detector_(0),
   nsteps_(0),
   navail_parts(0)
{
    totalen_=0;

	//create vector ntuple here
//	auto analysisManager = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}


bool B4aEventAction::checkConstruct(){
    if(!detector_ || bsmparticles_.size()<1)
        return false;


  //  const auto& activesensors=detector_->getActiveSensors();
 //   G4cout << "B4aEventAction::checkConstruct: " << activesensors->size() << " sensors active" << G4endl;
    if(hit_stopped_.size()){
        pdgids_.clear();
        for(const auto& ps: generator_->availParticles()){
            pdgids_.push_back(ps.second);

        }
        for(auto & p:hit_stopped_)
            p.second.clear();
        for(auto & p:hit_stopped_)
            p.second.resize(NACTIVELAYERS);

        return true;
    }

    for(const auto& ps: generator_->availParticles()){
        int s = ps.second;
        bool neg=s<0;
        if(neg)
            s*=-1;
        auto ss= std::to_string(s);
        if(neg)
            ss="m"+ss;
        hit_stopped_.push_back({ss, std::vector<int>(NACTIVELAYERS)});
    }
    for(const auto& ps: generator_->availParticles()){
        pdgids_.push_back(ps.second);
        std::cout << "particle " << ps.second << " registered"<<std::endl;
    }

    return true;
}

void B4aEventAction::accumulateVolumeInfo(G4VPhysicalVolume * volume,const G4Step* step){

	const auto& activesensors=detector_->getActiveSensors();

	bool issensor=true;

	nsteps_++;

	auto track = step->GetTrack();
	G4int trackid= track->GetTrackID();
	G4double momentum = track->GetMomentum().mag();
	size_t secondaries = step->GetNumberOfSecondariesInCurrentStep();
	if(false && !secondaries && !momentum){//primary
	    std::cout << trackid<< " " <<  track->GetDefinition()->GetPDGEncoding() << " status " << track->GetTrackStatus()<< " step "<< track->GetCurrentStepNumber() << std::endl;
	    std::cout << "secondaries: " << step->GetNumberOfSecondariesInCurrentStep()  << std::endl;
	    std::cout << "momentum "<< track->GetMomentum().mag() << std::endl;
	}


	G4int pdgid = track->GetDefinition()->GetPDGEncoding();

	//check if it's in the interesting particle ids
	//fStopAndKill

	auto stat = track->GetTrackStatus() ;

	bool done = stat == fStopButAlive || (!secondaries && !momentum);//stopped and no secondaries (stop and kill would be after decay)



	if(!done)
	    return; //only stopped interesting


	auto partit = std::find(bsmparticles_.begin(),bsmparticles_.end(),pdgid);
	if(partit == bsmparticles_.end()){
	    //if it has already decayed to SM particles, we don't care anymore, kill them
	    track->SetTrackStatus(fStopAndKill);// doesn't seem to bring many performance improvements

	    return; //not interesting
	}

   // std::cout << "BSM particle stopped " << pdgid << std::endl;

	size_t partidx = partit - bsmparticles_.begin();

    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(preStepPoint->GetTouchable());
    G4int copyNo = theTouchable->GetVolume()->GetCopyNo();

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

	int laynum = activesensors->at(idx).getLayer();

	if(laynum>=0){
	    //std::cout << "BSM particle stopped! " << pdgid ;
	    //std::cout << " --> in in detector layer "<< laynum<< "/"<<
	    //        hit_stopped_.at(partidx).second.size() <<std::endl;
        //
	    //auto energy=step->GetTotalEnergyDeposit();
	    hit_stopped_.at(partidx).second.at(laynum)++;
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

    //check if event is fine

    if(event->IsAborted()){//don't fill information
        G4cout << "event was aborted, not writing output" << G4endl;
        clear();
        return;
    }

    const auto& activesensors=detector_->getActiveSensors();
    if(!hit_layer_.size()){
        hit_layer_.resize(activesensors->size());
        for(size_t i=0;i<hit_layer_.size();i++){
            hit_layer_.at(i) = activesensors->at(i).getLayer();
        }
    }
    auto analysisManager = G4AnalysisManager::Instance();

  //  for(const auto& b: analysisManager->GetNtuple()->branches())
  //      std::cout << b->name() <<" " <<  std::endl;

    analysisManager->FillNtupleDColumn(0,generator_->getEnergy());




   // G4cout << "nsteps_ "<<nsteps_ <<G4endl;

    analysisManager->AddNtupleRow();

    clear();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
