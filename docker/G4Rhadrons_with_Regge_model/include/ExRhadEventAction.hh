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
// $Id: ExRhadEventAction.hh,v 1.3 2006/07/11 08:26:05 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExRhadEventAction_h
#define ExRhadEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class ExRhadEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExRhadEventAction : public G4UserEventAction
{
 public:
   ExRhadEventAction();
  ~ExRhadEventAction();

 public:
   void  BeginOfEventAction(const G4Event*);
   void    EndOfEventAction(const G4Event*);
    
   void AddAbs(G4double de, G4double dl) {EnergyAbs += de; TrackLAbs += dl;};
   void AddGap(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
                     
   void SetDrawFlag   (G4String val)  {drawFlag = val;};
   void SetPrintModulo(G4int    val)  {printModulo = val;};
    
 private:
   G4double  EnergyAbs, EnergyGap;
   G4double  TrackLAbs, TrackLGap;
                     
   G4String  drawFlag;
   G4int     printModulo;

   G4int evtNb;
                             
   ExRhadEventActionMessenger*  eventMessenger;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
