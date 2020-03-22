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
// $Id: B4aSteppingAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "globals.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
		const B4DetectorConstruction* detectorConstruction,
		B4aEventAction* eventAction)
: G4UserSteppingAction(),
  fDetConstruction(detectorConstruction),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
	// Collect energy and track length step by step

	// get volume of the current step
	auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

	//step->GetPreStepPoint()->GetTouchableHandle()->Get

	// energy deposit

    auto  transportMgr = G4TransportationManager::GetTransportationManager() ;

    auto fFieldPropagator = transportMgr->GetPropagatorInField() ;

    G4FieldManager* fieldMgr = fFieldPropagator->FindAndSetFieldManager( volume);
    G4Track* track = step->GetTrack();
    fieldMgr->ConfigureForTrack( track );


    bool looper = fFieldPropagator->IsParticleLooping() ;
    if(looper){//if anything is looping it's in vacuum
        //G4cout << "looper " <<G4endl;
        track->SetTrackStatus(fStopAndKill);
       // track->

    }

	fEventAction->accumulateVolumeInfo(volume, step);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
