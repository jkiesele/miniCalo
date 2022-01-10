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
// $Id: ExRhadSteppingAction.hh,v 1.2 2006/07/11 08:26:05 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExRhadSteppingAction_h
#define ExRhadSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4String.hh"

class ExRhadDetectorConstruction;
class ExRhadEventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExRhadSteppingAction : public G4UserSteppingAction
{
  public:
    ExRhadSteppingAction(ExRhadDetectorConstruction*, ExRhadEventAction*);
   ~ExRhadSteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
  ExRhadDetectorConstruction* detector;
  ExRhadEventAction*          eventaction;  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
