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
// $Id: ExRhadPrimaryGeneratorMessenger.cc,v 1.2 2005/11/02 09:01:37 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadPrimaryGeneratorMessenger.hh"

#include "ExRhadPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPrimaryGeneratorMessenger::ExRhadPrimaryGeneratorMessenger(
                                          ExRhadPrimaryGeneratorAction* ExRhadGun)
:ExRhadAction(ExRhadGun)
{
  gunDir = new G4UIdirectory("/Rhad/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  RndmCmd = new G4UIcmdWithAString("/Rhad/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/gun/position",this);
  PosCmd->SetGuidance("Set x-position of gun (helps to align with specific volume)");
  PosCmd->SetParameterName("Pos",false);
  PosCmd->SetUnitCategory("Length");
  PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BetaCmd = new G4UIcmdWithADouble("/Rhad/gun/beta",this);
  BetaCmd->SetGuidance("Set velocity");
  BetaCmd->SetParameterName("Beta",false);
  BetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FlatEnCmd = new G4UIcmdWithADoubleAndUnit("/Rhad/gun/flatenergy",this);
  FlatEnCmd->SetGuidance("Shoots particles randomly distributed in energy up to the given number");
  FlatEnCmd->SetParameterName("Energy",false);
  FlatEnCmd->SetUnitCategory("Energy");
  FlatEnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPrimaryGeneratorMessenger::~ExRhadPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete PosCmd;
  delete BetaCmd;
  delete FlatEnCmd;
  delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == RndmCmd ) ExRhadAction->SetRndmFlag(newValue);
  if( command == PosCmd ) ExRhadAction->SetGunPos(PosCmd->GetNewDoubleValue(newValue));
  if( command == BetaCmd ) ExRhadAction->SetBeta(BetaCmd->GetNewDoubleValue(newValue));
  if( command == FlatEnCmd ) ExRhadAction->SetFlatEn(FlatEnCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

