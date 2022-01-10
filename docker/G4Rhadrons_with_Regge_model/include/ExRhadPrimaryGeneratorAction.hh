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
// $Id: ExRhadPrimaryGeneratorAction.hh,v 1.2 2005/11/02 09:01:37 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExRhadPrimaryGeneratorAction_h
#define ExRhadPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class ExRhadDetectorConstruction;
class ExRhadPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExRhadPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ExRhadPrimaryGeneratorAction(ExRhadDetectorConstruction*);    
  ~ExRhadPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetGunPos(G4double val) { customGunPos = true; xGun = val; }
  void SetBeta(G4double val);
  void SetFlatEn(G4double val);

private:
  G4ParticleGun*                particleGun;	  //pointer a to G4  class
  ExRhadDetectorConstruction*    ExRhadDetector;  //pointer to the geometry
  
  ExRhadPrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  G4String                      rndmFlag;	  //flag for a rndm impact point
  G4bool                        customGunPos;
  G4bool                        flaten;
  G4double                      emax;
  G4double                      xGun;

  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


