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

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "G4StepPoint.hh"
#include "sensorContainer.h"
#include "B4DetectorConstruction.hh"
#include "G4Step.hh"
#include "B4RunAction.hh"
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
    	rechit_energy_.clear();
    	allvolumes_.clear();
    	rechit_absorber_energy_.clear();
        rechit_x_.clear();
        rechit_y_.clear();
        rechit_z_.clear();
        rechit_vz_.clear();
        rechit_varea_.clear();
        rechit_vxy_.clear();
        rechit_layer_.clear();
        rechit_detid_.clear();
        nsteps_=0;
        totalen_=0;
    }


    void setGenerator(B4PrimaryGeneratorAction * generator){
    	generator_=generator;
    }
    void setDetector(B4DetectorConstruction * detector){
    	detector_=detector;
    }

  private:
    G4double  fEnergyAbs;
    G4double totalen_;
    std::vector<G4double>  rechit_energy_,rechit_absorber_energy_;
    std::vector<G4double>  rechit_x_;
    std::vector<G4double>  rechit_y_;
    std::vector<G4double>  rechit_z_;
    std::vector<G4double>  rechit_layer_;
    std::vector<G4double>  rechit_vz_;
    std::vector<G4double>  rechit_varea_;
    std::vector<G4double>  rechit_vxy_;
    std::vector<int>       rechit_detid_;
    std::vector<const G4VPhysicalVolume * > allvolumes_;

    G4double  fEnergyGap;
    G4double  fTrackLAbs; 
    G4double  fTrackLGap;


    B4PrimaryGeneratorAction * generator_;
    B4DetectorConstruction * detector_;

    size_t nsteps_;

};

// inline functions



inline void B4aEventAction::AddEnergy(G4double de, G4double dl) {
  fEnergyGap += de; 
  fTrackLGap += dl;
}
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
