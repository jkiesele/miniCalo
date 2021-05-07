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
// $Id: B4aEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "defines.h"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "B4PartGeneratorBase.hh"
#include "G4StepPoint.hh"
#include "sensorContainer.h"
#include "B4DetectorConstruction.hh"
#include "G4Step.hh"
#include "B4RunAction.hh"
#include "B4PrimaryGeneratorAction.hh"
#include <algorithm>

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()
class G4VPhysicalVolume;
class B4aEventAction : public G4UserEventAction
{
	friend B4RunAction;
  public:
    B4aEventAction();
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddEnergy(G4double de, G4double dl);
    

    void accumulateVolumeInfo(G4VPhysicalVolume *,const G4Step* );

    void clear(){
    	allvolumes_.clear();
    	//for(auto& v: hit_stopped_)
    	 //   for(auto& vv:v.second)
    	  //      vv=0;
    	//hit_layer_.clear();
        nsteps_=0;
        totalen_=0;
    }


    void setGenerator(B4PartGeneratorBase  * generator){
    	generator_=generator;
    	auto checkparts = generator_->availParticles();
    	for(const auto p: checkparts){
    	    bsmparticles_.push_back(p.second);
    	    std::cout << "added from generator " << p.second << std::endl;
    	}
    	navail_parts=bsmparticles_.size();
    	checkConstruct();
    }
    void setDetector(B4DetectorConstruction * detector){
    	detector_=detector;
    	checkConstruct();
    }

    bool checkConstruct();

    size_t nevents_;

  private:
    G4double  fEnergyAbs;
    G4double totalen_;

    //particle in layer


    std::vector<int> pdgids_;
    std::vector<std::pair<G4String, std::vector< int> > > hit_stopped_;
    std::vector<int>  hit_layer_;


    std::vector<const G4VPhysicalVolume * > allvolumes_;

    G4double  fEnergyGap;
    G4double  fTrackLAbs; 
    G4double  fTrackLGap;


    B4PartGeneratorBase * generator_;
    B4DetectorConstruction * detector_;

    size_t nsteps_;
    size_t navail_parts;

    std::vector<int> bsmparticles_;

};

// inline functions



inline void B4aEventAction::AddEnergy(G4double de, G4double dl) {
  fEnergyGap += de; 
  fTrackLGap += dl;
}
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
