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
// $Id: B4PrimaryGeneratorAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4INCLRandom.hh"
#include <G4INCLGeant4Random.hh>
#include <G4INCLRandomSeedVector.hh>
#include<ctime>
#include<sys/types.h>

template<class T>
static G4String createString(const T& i){
    std::stringstream ss;
    ss << i;
    std::string number=ss.str();
    return number;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : B4PartGeneratorBase(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
 // G4ParticleTable::GetParticleTable()->DumpTable();
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(100.*GeV);


  G4INCL::Random::setGenerator( new G4INCL::Geant4RandomGenerator());

  for(int i=0;i<seedsoffset_;i++)
      G4double rand =  G4INCL::Random::shoot();

  xorig_=0;
  yorig_=0;
  setParticleID(gamma);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}
std::vector<G4String> B4PrimaryGeneratorAction::generateAvailableParticles()const{
	std::vector<G4String> out;
	for(int i=0;i<particles_size;i++)
		out.push_back(getParticleName((B4PrimaryGeneratorAction::particles)i));

	return out;
}


G4String B4PrimaryGeneratorAction::getParticleName(enum particles p)const{
    if(p==elec){
        return "isElectron";
    }
    if(p==positron){
        return "isPositron";
    }
    if(p==muon){
        return "isMuon";
    }
    if(p==pioncharged){
        return "isPionCharged";
    }
    if(p==pionneutral){
        return "isPionNeutral";
    }
    if(p==klong){
        return "isK0Long";
    }
    if(p==kshort){
        return "isK0Short";
    }
    if(p==gamma){
        return "isGamma";
    }

    return "isInvalid"+createString((int)p);
}

G4String B4PrimaryGeneratorAction::setParticleID(enum particles p){
	particleid_=p;
	G4ParticleDefinition * particleDefinition =0;

	if(p==elec){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("e-");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==positron){
	    particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("e+");
	    fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==muon){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==pioncharged){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==pionneutral){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==klong){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("kaon0L");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==kshort){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("kaon0S");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}
	if(p==gamma){
		particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
		fParticleGun->SetParticleDefinition(particleDefinition);
	}

	return getParticleName(p);
	return "isInvalid"+createString((int)p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void generatePosition(G4double dRmin,G4double dRmax, G4double&x , G4double&y){

    while(true){
        x = dRmax - 2*dRmax*G4INCL::Random::shoot();
        y = dRmax - 2*dRmax*G4INCL::Random::shoot();
        G4ThreeVector dir(x,y,0);
        G4double dr = dir.r();
        if(dr > dRmin && dr <= dRmax)
            break;
    }
}
double gen_etaToR(const G4double& eta, const G4double& z){
    return z * exp(-eta);
}

//just make sure it hits the first calo layer
G4ThreeVector generateDirection(G4double etamin,G4double etamax, const G4ThreeVector& position){

    G4double calo_z = 30*cm+300*cm;//should hit all layers

    G4ThreeVector initialdir = position.unit();

    G4double dz = calo_z-position.z();

    while(true){
        G4double x_dir= 1.-2.*G4INCL::Random::shoot();
        G4double y_dir= 1.-2.*G4INCL::Random::shoot();
        G4double z_dir= G4INCL::Random::shoot()+0.01;

        G4ThreeVector dir(x_dir,y_dir,z_dir);
        dir=dir.unit();
        dir *= dz/dir.z();

        G4ThreeVector newpos = position + dir;

        G4double eta = newpos.eta();
        if(eta < etamax && eta > etamin)
            return dir.unit();
    }
}

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  G4cout << "shooting.. " <<G4endl;
  // Set gun position


  double energy_max=200;
  double energy_min=10;
  //energy_=15;

  //iterate
//  if(particleid_==gamma)
//      G4cout << setParticleID(elec) <<G4endl;
//  else if(particleid_==elec)
//      G4cout << setParticleID(gamma) <<G4endl;
 // else if(particleid_==gamma)
  //    G4cout << setParticleID(pioncharged) <<G4endl;

  setParticleID(gamma);
  //positron


// good metrics are: angle w.r.t. projective and energy



  //particleid_=pioncharged;

  energy_=10001;
  while(energy_>energy_max){//somehow sometimes the random gen shoots >1??
      G4double rand =  G4INCL::Random::shoot();
      energy_=(energy_max-energy_min)*rand+energy_min;
  }


  G4cout << "get position " << G4endl;
  generatePosition(20.*cm,60*cm,xorig_,yorig_);//15.*cm,66*cm,xorig_,yorig_);

  G4ThreeVector position(xorig_,yorig_,299*cm);//xorig_,yorig_,0);
  G4double x_dir,y_dir,z_dir;


  G4ThreeVector direction = generateDirection(1.55,2.95,position);

  dirx_ = direction.x();
  diry_ = direction.y();
  dirz_ = direction.z();

  G4cout << "projective direction " << position.unit() << G4endl;


  diff_proj_phi_=direction.deltaPhi(position);
  diff_proj_theta_=direction.theta(position);

  angle_ = direction.angle(position.unit());

  G4cout << "Dphi,DTheta " <<  diff_proj_phi_ << ", " << diff_proj_theta_ << " angle " << angle_ << G4endl;

  G4cout << "new direction " << direction << G4endl;

  //get eta

//0.2*G4INCL::Random::shoot(),0.2*G4INCL::Random::shoot(),0);


  fParticleGun->SetParticleMomentumDirection(direction);
  fParticleGun->SetParticleEnergy(energy_ * GeV);
  fParticleGun->SetParticlePosition(position);
  fParticleGun->GeneratePrimaryVertex(anEvent);


  G4cout << "Energy: " << energy_ << ", position: "<< position <<" direction "<<direction<< " field "<<1*tesla <<  G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

