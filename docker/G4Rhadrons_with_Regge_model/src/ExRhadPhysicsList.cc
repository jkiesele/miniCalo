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
// $Id: ExRhadPhysicsList.cc,v 1.12 2006/07/11 08:26:06 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"




using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPhysicsList::ExRhadPhysicsList():  G4VPhysicsConstructor()
{
// defaultCutValue = 1.0*mm;
 SetVerboseLevel(1);

 G4cout << "initialised ExRhadPhysicsList" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadPhysicsList::~ExRhadPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructNuclei();
  ConstructExotics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

void ExRhadPhysicsList::ConstructNuclei()
{
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
}
#include "CustomParticleFactory.h"

void ExRhadPhysicsList::ConstructExotics()
{
  CustomParticleFactory::loadCustomParticles();     
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructProcess()
{
 // AddTransportation();
  //  ConstructEM();
  addCustomPhysics();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::ConstructEM()
{
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxAlongStep,3);      
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);

    } else if (particleName == "e+") {
    //positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxAlongStep,3);      
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
    //muon
      G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
      G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction     = new G4MuPairProduction();
      G4VProcess* anIonisation        = new G4MuIonisation();
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxAlongStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxAlongStep,4);      
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxPostStep,4);

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino" &&
		particle->GetParticleType() !="rhadron" &&
		particle->GetParticleType() !="Custom" &&
		particle->GetParticleName() !="pi-" &&
		particle->GetParticleName() !="pi+" 
		// These pion lines used with ConstructCustom to cross check with pion behaviour
		)) {
      G4cout<<particle->GetParticleName()<<G4endl;
     // all others charged particles except geantino
     G4VProcess* aMultipleScattering = new G4hMultipleScattering();
     G4VProcess* anIonisation        = new G4hIonisation();
     //
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     //
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
     //
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
    }
  }
}

//#include "G4LElastic.hh"
#include "G4HadronElasticProcess.hh"
#include "FullModelHadronicProcess.hh"
#include "ToyModelHadronicProcess.hh"

#include "G4PionMinusInelasticProcess.hh"
//#include "G4LEPionMinusInelastic.hh"
//#include "G4HEPionMinusInelastic.hh"

#include "G4PionPlusInelasticProcess.hh"
//#include "G4LEPionPlusInelastic.hh"
//#include "G4HEPionPlusInelastic.hh"
#include "G4DecayTable.hh"
#include "G4hhIonisation.hh"

void ExRhadPhysicsList::addCustomPhysics()
{
  G4cout << " CustomPhysics: adding CustomPhysics processes  " <<G4endl;
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();

  while((*theParticleIterator)())
    {	
      int i = 0;
      G4ParticleDefinition* particle = theParticleIterator->value();
      //      G4cout<<particle->GetParticleName()<<G4endl;
      CustomParticle* cp = dynamic_cast<CustomParticle*>(particle);
      if(CustomParticleFactory::isCustomParticle(particle))
	{
	  G4cout << particle->GetParticleName()<<", "<<particle->GetPDGEncoding() 
		 << " is Custom. Mass is "
		 <<particle->GetPDGMass()/GeV
		 <<" GeV."<<G4endl;

	  G4DecayTable* table = particle->GetDecayTable();
	  G4cout<<"table: "<<table<<G4endl;
	  if (table!=0) table->DumpInfo();

	  if(cp->GetCloud()!=0)
	    {
	      G4cout<<"Cloud mass is "
		    <<cp->GetCloud()->GetPDGMass()/GeV
		    <<" GeV. Spectator mass is "
		    <<static_cast<CustomParticle*>(particle)->GetSpectator()->GetPDGMass()/GeV
		    <<" GeV." 
		    << G4endl;
	    }
	  G4ProcessManager* pmanager = particle->GetProcessManager();
	  if(pmanager)
	    { 
              if(cp!=0) {
		pmanager->AddDiscreteProcess(new FullModelHadronicProcess()); // Full parametrised model
		//		pmanager->AddDiscreteProcess(new FullModelHadronicProcess()); // Toy-model
	      }
              if(particle->GetPDGCharge()/eplus != 0)
		{ 
		  pmanager->AddProcess(new G4hMultipleScattering,-1, 1,i+1);
		  pmanager->AddProcess(new G4hhIonisation,       -1, 2,i+2);
		}
	      pmanager->DumpInfo();
            }
	  else
	    G4cout << "   No pmanager" << G4endl;
	}


    }
}












#include "G4Decay.hh"

void ExRhadPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadPhysicsList::SetCuts()
{
  //if (verboseLevel >0){
  //  G4cout << "ExRhadPhysicsList::SetCuts:";
  //  G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  //}
  //
  //// set cut values for gamma at first and for e- second and next for e+,
  //// because some processes for e+/e- need cut values for gamma
  ////
  //SetCutValue(defaultCutValue, "gamma");
  //SetCutValue(defaultCutValue, "e-");
  //SetCutValue(defaultCutValue, "e+");
  //
  //if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

