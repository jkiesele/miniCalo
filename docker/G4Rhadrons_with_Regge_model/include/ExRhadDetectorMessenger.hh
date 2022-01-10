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
// $Id: ExRhadDetectorMessenger.hh,v 1.2 2005/07/23 11:59:49 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExRhadDetectorMessenger_h
#define ExRhadDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ExRhadDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExRhadDetectorMessenger: public G4UImessenger
{
public:
  ExRhadDetectorMessenger(ExRhadDetectorConstruction* );
  ~ExRhadDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  ExRhadDetectorConstruction* ExRhadDetector;
    
  G4UIdirectory*             RhadDir;
  G4UIdirectory*             detDir;
  
  G4UIcmdWithADoubleAndUnit*   WorldSizeCmd;

  G4UIcmdWithADoubleAndUnit*   Box1SizeCmd;
  G4UIcmdWithADoubleAndUnit*   Box2SizeCmd;
  G4UIcmdWithADoubleAndUnit*   Box3SizeCmd;
  G4UIcmdWithADoubleAndUnit*   Box4SizeCmd;
  G4UIcmdWithADoubleAndUnit*   DetectorSpacingCmd;

  G4UIcmdWithADoubleAndUnit*   Box1PosCmd;
  G4UIcmdWithADoubleAndUnit*   Box2PosCmd;
  G4UIcmdWithADoubleAndUnit*   Box3PosCmd;
  G4UIcmdWithADoubleAndUnit*   Box4PosCmd;

  G4UIcmdWithAString*        Box1MaterCmd;
  G4UIcmdWithAString*        Box2MaterCmd;
  G4UIcmdWithAString*        Box3MaterCmd;
  G4UIcmdWithAString*        Box4MaterCmd;



  G4UIcmdWithAString*        AbsMaterCmd;
  G4UIcmdWithAString*        GapMaterCmd;
  G4UIcmdWithADoubleAndUnit* AbsThickCmd;
  G4UIcmdWithADoubleAndUnit* GapThickCmd;
  G4UIcmdWithADoubleAndUnit* SizeYZCmd;
  G4UIcmdWithAnInteger*      NbLayersCmd;    
  G4UIcmdWithADoubleAndUnit* MagFieldCmd;
  G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

