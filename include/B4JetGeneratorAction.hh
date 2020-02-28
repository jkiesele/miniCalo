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
// $Id: B4JetGeneratorAction.hh 94808 2015-12-10 08:22:26Z gcosmo $
//
/// \file B4JetGeneratorAction.hh
/// \brief Definition of the B4JetGeneratorAction class

#ifndef B4JetGeneratorAction_h
#define B4JetGeneratorAction_h 1

#include "defines.h"

#ifndef NOPYTHIA

#include "B4PrimaryGeneratorAction.hh"
#include "G4PrimaryVertex.hh"
#include "B4PartGeneratorBase.hh"
#include "globals.hh"
#include <vector>

#include "Pythia8/Pythia.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

unsigned const NMAX(1024);

class G4ParticleGun;
class G4Event;
class TTree;

/// The primary generator action class with pythia + fastjet.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).



class B4JetGeneratorAction : public B4PartGeneratorBase
{
public:
  B4JetGeneratorAction(particles p);
  virtual ~B4JetGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  void GenerateSingleVertex(G4PrimaryVertex* vertex);

  G4double getEnergy()const{return energy_;}

  G4double getX()const{return xorig_;}
  G4double getY()const{return yorig_;}
  G4double getR()const{return std::sqrt(yorig_*yorig_+xorig_*xorig_);}

  bool isJetGenerator(){return true;}

  std::vector<G4String> generateAvailableParticles()const{return {"isMinbias","isDisplacedJet"};}

  particles getParticle() const{
    return particle_;
  }

  int isParticle(int i)const{
    if (i == int(particle_))
      return 1;
    else
      return 0;
  }

private:
  G4double energy_;
  G4double xorig_,yorig_;

  Pythia8::Pythia pythia_;
  particles particle_;
  fastjet::JetDefinition* jetDef_;
  std::vector<fastjet::PseudoJet> fjinputs_;

  bool firstEvent_{true};

  TTree* truthTree_{nullptr};
  //mkbranch
  int jetType_;
  float jetE_;
  float jetEta_;
  unsigned nPart_;
  int partPid_[NMAX];
  float partE_[NMAX];
  float partEta_[NMAX];
  float partPhi_[NMAX];
  //mkbranch

  int nPU_;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
#endif
