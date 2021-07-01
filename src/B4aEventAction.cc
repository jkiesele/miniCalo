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
 : G4UserEventAction()

{

    nevents_=0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aEventAction::~B4aEventAction()
{}


bool B4aEventAction::checkConstruct(){


    return true;
}


void B4aEventAction::accumulateParticleInfo(const G4Track * track,
        const G4ThreeVector& pos){

    auto pol = track->GetPolarization();
    auto mom = track->GetMomentum();


  //  std::cout << "step at " << pos << " momentum " << mom << std::endl;

    if(ps_.size()){
        if(mom == ps_.at(ps_.size()-1)) //no change, just boundary crossing
            return;
    }
    ps_.push_back(mom);
    pols_.push_back(pol);
    poss_.push_back(pos);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static void fillflatvector(std::vector<std::vector<double> * > tofill,
        const std::vector<G4ThreeVector>& data){
    if(tofill.size()!=3)
        throw std::out_of_range("fillflatvector");
    for(const auto& v:data){
        tofill.at(0)->push_back(v.x());
        tofill.at(1)->push_back(v.y());
        tofill.at(2)->push_back(v.z());
    }

}

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

    //fill into flat vectors
   // t_px, t_py, t_pz, t_polx, t_poly, t_polz, t_posx, t_posy, t_posz

   fillflatvector({&t_px,   &t_py,   &t_pz}, ps_);
   fillflatvector({&t_polx, &t_poly, &t_polz}, pols_);
   fillflatvector({&t_posx, &t_posy, &t_posz}, poss_);


    auto analysisManager = G4AnalysisManager::Instance();


    analysisManager->AddNtupleRow();

    nevents_++;

    clear();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
