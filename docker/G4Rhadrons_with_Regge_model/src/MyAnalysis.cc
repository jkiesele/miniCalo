#include "MyAnalysis.hh"
#include "G4RunManager.hh" 


MyAnalysis* MyAnalysis::instance = 0;

 
MyAnalysis::MyAnalysis() {


  Nevt=0;
}


MyAnalysis::~MyAnalysis() {


}


MyAnalysis* MyAnalysis::getInstance() {
  if ( instance == 0 ) instance = new MyAnalysis();
  return instance;
}

//void MyAnalysis::fill(const G4int numHits, const G4double energyDeposited, const G4double distance) {
void MyAnalysis::fill(const MyTrackerHitsCollection* MyHits) {

  //This is perhaps a bad way of counting events.
  Nevt++;



  G4int numTrackerHits = MyHits->entries();
  G4cout<<"Number of tracker hits: "<<numTrackerHits<<G4endl;

  for ( G4int i=0; i < numTrackerHits; i++ ) {
    if ( (*MyHits)[i]->GetParticleDefinition()->GetParticleType()=="rhadron" )// Only looking at R-hadron
      {
	G4double xpos = (*MyHits)[i]->GetPosition().x();
	G4double e = (*MyHits)[i]->Get4Mom().e();
	G4double edep = (*MyHits)[i]->GetEdep();
      }
    /*
    if((*MyHits)[i]->GetParticleDefinition()->GetParticleType()=="rhadron")
      G4cout <<"Hit was: "<<(*MyHits)[i]->GetParticleDefinition()->GetParticleName()<<G4endl;
    */
  }

}



