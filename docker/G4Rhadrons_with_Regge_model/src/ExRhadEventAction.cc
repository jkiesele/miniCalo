//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExRhadEventAction.cc,v 1.4 2006/07/11 08:26:06 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadEventAction.hh"
#include "ExRhadEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//***LOOKHERE***
#include "MyAnalysis.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "MySensitiveTracker.hh"
#include "MyTrackerHit.hh"
//**************

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadEventAction::ExRhadEventAction()
:drawFlag("all"),printModulo(1),eventMessenger(0)
{
  eventMessenger = new ExRhadEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadEventAction::~ExRhadEventAction()
{
  delete eventMessenger;
  MyAnalysis* analysis = MyAnalysis::getInstance();
  analysis->~MyAnalysis();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadEventAction::BeginOfEventAction(const G4Event* evt)
{  
 evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) { 
   G4cout << "\n---> Begin of event: " << evtNb << G4endl;
   CLHEP::HepRandom::showEngineStatus();
 }
 
 // initialisation per event
 EnergyAbs = EnergyGap = 0.;
 TrackLAbs = TrackLGap = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadEventAction::EndOfEventAction(const G4Event* evt)
{

  MyAnalysis* analysis = MyAnalysis::getInstance();

  // Get the IDs of the collections of tracker hits. 
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int trackerCollID = SDman->GetCollectionID(colNam="trackerCollection");

  // Extract the hits.
  G4HCofThisEvent* collectionEventHits = evt->GetHCofThisEvent();
  MyTrackerHitsCollection* collectionTrackerHits = 0;

  if ( collectionEventHits ) {
    collectionTrackerHits = dynamic_cast< MyTrackerHitsCollection* >
      ( collectionEventHits->GetHC(trackerCollID) );
  }
  if ( collectionTrackerHits ) {
    analysis->fill(collectionTrackerHits);
  }

  // extract the trajectories and draw them

  // You can get a default drawing without this code by using, e.g.,
  // /vis/scene/add/trajectories 1000
  // The code here adds sophistication under control of drawFlag.

  // See comments in G4VTrajectory::DrawTrajectory for the
  // interpretation of the argument, 1000.
  
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager)
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") pVisManager->Draw(*trj);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  pVisManager->Draw(*trj);
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                  pVisManager->Draw(*trj);
        }
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
