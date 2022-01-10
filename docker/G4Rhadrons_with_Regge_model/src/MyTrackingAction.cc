#include "MyTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "CustomParticle.h"

MyTrackingAction::MyTrackingAction() :
  counter(0) {
}


void MyTrackingAction::PostUserTrackingAction(const G4Track* t) {
  if(dynamic_cast<CustomParticle*>(t->GetDefinition())!=0
     &&
     t->GetDynamicParticle()->GetKineticEnergy()==0.)
  return;
}
