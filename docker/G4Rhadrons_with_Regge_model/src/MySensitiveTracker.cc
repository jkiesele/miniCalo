#include "MySensitiveTracker.hh"
#include "MyTrackerHit.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

using namespace CLHEP;

MySensitiveTracker::MySensitiveTracker(G4String name)
  : G4VSensitiveDetector(name), trackerCollection(0) {
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}


MySensitiveTracker::~MySensitiveTracker() {}


void MySensitiveTracker::Initialize(G4HCofThisEvent* HCE) {
  trackerCollection = new MyTrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
}


G4bool MySensitiveTracker::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4Track* aTrack = aStep->GetTrack();

  const G4ThreeVector* tpos = new G4ThreeVector(aStep->GetPostStepPoint()->GetPosition());
  /*
  if (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="GenericHadronicProcess")
    {
      G4cout<<"Interesting hit. edep: "<<edep<<G4endl;
    }
  */
  if( edep < 0.001*eV ) {delete tpos; return false;}
  /*
  if (aTrack->GetDynamicParticle()->GetDefinition()->GetParticleType()=="rhadron")
    G4cout<<"Making a hit of: "<<aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName()<<G4endl;
  */

  MyTrackerHit* trackerHit = new MyTrackerHit();
  trackerHit->SetEdep( edep );
  trackerHit->SetPos( *tpos  );
  trackerHit->Set4Mom(aTrack->GetDynamicParticle()->Get4Momentum());
  trackerHit->SetParticleDefinition(aTrack->GetDynamicParticle()->GetDefinition());//Er åbenbart et sekundært track...
  trackerCollection->insert( trackerHit );

  delete tpos;

  return true;
}


void MySensitiveTracker::EndOfEvent(G4HCofThisEvent* HCE) {
  static G4int HCID = -1;
  if ( HCID < 0 ) { 
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection( HCID, trackerCollection );
  if ( ! trackerCollection->entries() ) {
    G4cout << " MySensitiveTracker::EndOfEvent  >>> No tracker hit <<< " << G4endl;
  }   
}


void MySensitiveTracker::clear() {} 


void MySensitiveTracker::DrawAll() {} 


void MySensitiveTracker::PrintAll() {} 

