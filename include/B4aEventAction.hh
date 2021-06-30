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
#include "G4ParticleGun.hh"

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
    
    inline void readInitialParticleInfo(){
        auto mom = generator_->getGun()->GetParticleMomentumDirection();
        mom *= generator_->getGun()->GetParticleMomentum();
        setInitialMomentum(mom);

        setInitialPol(generator_->getGun()->GetParticlePolarization());

    }

    void accumulateParticleInfo(const G4Track *);

    inline void setInitialMomentum(const G4ThreeVector& vec){
        in_px = vec.x();
        in_py = vec.y();
        in_pz = vec.z();
    }
    inline void setInitialPol(const G4ThreeVector& vec){
        in_polx = vec.x();
        in_poly = vec.y();
        in_polz = vec.z();
    }
    inline void setOutMomentum(const G4ThreeVector& vec){
        out_px = vec.x();
        out_py = vec.y();
        out_pz = vec.z();
    }
    inline void setOutPol(const G4ThreeVector& vec){
        out_polx = vec.x();
        out_poly = vec.y();
        out_polz = vec.z();
    }



    void clear(){

    }


    void setGenerator(B4PartGeneratorBase  * generator){
        generator_=generator;
        generator_->evtact = this;
    }
    void setDetector(B4DetectorConstruction * detector){
        detector_=detector;
        checkConstruct();
    }

    bool checkConstruct();

    size_t nevents_;

  private:

    G4double out_px,out_py,out_pz;
    G4double out_polx,out_poly,out_polz;

    G4double in_polx,in_poly,in_polz;
    G4double in_px,in_py,in_pz;

    B4PartGeneratorBase * generator_;
    B4DetectorConstruction * detector_;


};

// inline functions


                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
