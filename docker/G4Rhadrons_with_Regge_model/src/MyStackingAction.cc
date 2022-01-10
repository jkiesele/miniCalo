#include "MyStackingAction.hh"
#include "G4Track.hh"
#include "CustomParticle.h"

MyStackingAction::MyStackingAction(){}

MyStackingAction::~MyStackingAction(){}

G4ClassificationOfNewTrack MyStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(dynamic_cast<CustomParticle*>(aTrack->GetDefinition())!=0)
    return fUrgent;
  
  return fKill;
}
