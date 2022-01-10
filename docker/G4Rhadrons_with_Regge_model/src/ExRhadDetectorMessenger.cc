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
// $Id: ExRhadDetectorMessenger.cc,v 1.2 2005/07/23 11:59:49 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadDetectorMessenger.hh"

#include "ExRhadDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadDetectorMessenger::ExRhadDetectorMessenger(
					       ExRhadDetectorConstruction* ExRhadDet)
  :ExRhadDetector(ExRhadDet)
{ 
  RhadDir = new G4UIdirectory("/Rhad/");
  RhadDir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/Rhad/det/");
  detDir->SetGuidance("detector control");

  WorldSizeCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setWorldSize",this);
  WorldSizeCmd->SetGuidance("Set size of World volume. Necessary if you enlargen volumes.");
  WorldSizeCmd->SetParameterName("Size",false);
  WorldSizeCmd->SetRange("Size>=0.");
  WorldSizeCmd->SetUnitCategory("Length");
  WorldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  Box1SizeCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox1Size",this);
  Box1SizeCmd->SetGuidance("Set size of Box 1");
  Box1SizeCmd->SetParameterName("Size",false);
  Box1SizeCmd->SetRange("Size>=0.");
  Box1SizeCmd->SetUnitCategory("Length");
  Box1SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box2SizeCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox2Size",this);
  Box2SizeCmd->SetGuidance("Set size of Box 2");
  Box2SizeCmd->SetParameterName("Size",false);
  Box2SizeCmd->SetRange("Size>=0.");
  Box2SizeCmd->SetUnitCategory("Length");
  Box2SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box3SizeCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox3Size",this);
  Box3SizeCmd->SetGuidance("Set size of Box 3");
  Box3SizeCmd->SetParameterName("Size",false);
  Box3SizeCmd->SetRange("Size>=0.");
  Box3SizeCmd->SetUnitCategory("Length");
  Box3SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  Box4SizeCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox4Size",this);
  Box4SizeCmd->SetGuidance("Set size of Box 4");
  Box4SizeCmd->SetParameterName("Size",false);
  Box4SizeCmd->SetRange("Size>=0.");
  Box4SizeCmd->SetUnitCategory("Length");
  Box4SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box1PosCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox1Pos",this);
  Box1PosCmd->SetGuidance("Set position of Box 1");
  Box1PosCmd->SetParameterName("Pos",false);
  Box1PosCmd->SetUnitCategory("Length");
  Box1PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box2PosCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox2Pos",this);
  Box2PosCmd->SetGuidance("Set position of Box 2");
  Box2PosCmd->SetParameterName("Pos",false);
  Box2PosCmd->SetUnitCategory("Length");
  Box2PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box3PosCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox3Pos",this);
  Box3PosCmd->SetGuidance("Set position of Box 3");
  Box3PosCmd->SetParameterName("Pos",false);
  Box3PosCmd->SetUnitCategory("Length");
  Box3PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  Box4PosCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setBox4Pos",this);
  Box4PosCmd->SetGuidance("Set position of Box 4");
  Box4PosCmd->SetParameterName("Pos",false);
  Box4PosCmd->SetUnitCategory("Length");
  Box4PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DetectorSpacingCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setDetectorSpacing",this);
  DetectorSpacingCmd->SetGuidance("Set detector spacing and resize world accordingly");
  DetectorSpacingCmd->SetParameterName("Size",false);
  DetectorSpacingCmd->SetRange("Size>=0.");
  DetectorSpacingCmd->SetUnitCategory("Length");
  DetectorSpacingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
         
  Box1MaterCmd = new G4UIcmdWithAString("/Rhad/det/setBox1Mat",this);
  Box1MaterCmd->SetGuidance("Select Material of Box1.");
  Box1MaterCmd->SetGuidance("Default: Scintillator");
  Box1MaterCmd->SetParameterName("choice",false);
  Box1MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box2MaterCmd = new G4UIcmdWithAString("/Rhad/det/setBox2Mat",this);
  Box2MaterCmd->SetGuidance("Select Material of Box2.");
  Box2MaterCmd->SetGuidance("Default: Iron");
  Box2MaterCmd->SetParameterName("choice",false);
  Box2MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box3MaterCmd = new G4UIcmdWithAString("/Rhad/det/setBox3Mat",this);
  Box3MaterCmd->SetGuidance("Select Material of Box3.");
  Box3MaterCmd->SetGuidance("Default: Air");
  Box3MaterCmd->SetParameterName("choice",false);
  Box3MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Box4MaterCmd = new G4UIcmdWithAString("/Rhad/det/setBox4Mat",this);
  Box4MaterCmd->SetGuidance("Select Material of Box4.");
  Box4MaterCmd->SetGuidance("Default: Iron");
  Box4MaterCmd->SetParameterName("choice",false);
  Box4MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  //Remember to keep this
  UpdateCmd = new G4UIcmdWithoutParameter("/Rhad/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/det/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadDetectorMessenger::~ExRhadDetectorMessenger()
{
  delete WorldSizeCmd;
  delete Box1SizeCmd;
  delete Box2SizeCmd;
  delete Box3SizeCmd;
  delete Box4SizeCmd;
  delete Box1PosCmd;
  delete Box2PosCmd;
  delete Box3PosCmd;
  delete Box4PosCmd;
  delete Box1MaterCmd;
  delete Box2MaterCmd;
  delete Box3MaterCmd;
  delete Box4MaterCmd;
  delete UpdateCmd;
  delete MagFieldCmd;
  delete detDir;
  delete RhadDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == WorldSizeCmd ) ExRhadDetector->SetWorldSize(WorldSizeCmd->GetNewDoubleValue(newValue));

  if( command == Box1SizeCmd ) ExRhadDetector->SetBox1Size(Box1SizeCmd->GetNewDoubleValue(newValue));

  if( command == Box2SizeCmd ) ExRhadDetector->SetBox2Size(Box2SizeCmd->GetNewDoubleValue(newValue));

  if( command == Box3SizeCmd ) ExRhadDetector->SetBox3Size(Box3SizeCmd->GetNewDoubleValue(newValue));

  if( command == Box4SizeCmd ) ExRhadDetector->SetBox4Size(Box4SizeCmd->GetNewDoubleValue(newValue));

  if( command == Box1PosCmd ) ExRhadDetector->SetBox1Pos(Box1PosCmd->GetNewDoubleValue(newValue));

  if( command == Box2PosCmd ) ExRhadDetector->SetBox2Pos(Box2PosCmd->GetNewDoubleValue(newValue));

  if( command == Box3PosCmd ) ExRhadDetector->SetBox3Pos(Box3PosCmd->GetNewDoubleValue(newValue));

  if( command == Box4PosCmd ) ExRhadDetector->SetBox4Pos(Box4PosCmd->GetNewDoubleValue(newValue));

  if( command == DetectorSpacingCmd ) ExRhadDetector->SetDetectorSpacing(DetectorSpacingCmd->GetNewDoubleValue(newValue));

  if( command == Box1MaterCmd ) ExRhadDetector->SetBox1Material(newValue);

  if( command == Box2MaterCmd ) ExRhadDetector->SetBox2Material(newValue);

  if( command == Box3MaterCmd ) ExRhadDetector->SetBox3Material(newValue);

  if( command == Box4MaterCmd ) ExRhadDetector->SetBox4Material(newValue);

  if( command == UpdateCmd ) ExRhadDetector->UpdateGeometry();

  if( command == MagFieldCmd ) ExRhadDetector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
