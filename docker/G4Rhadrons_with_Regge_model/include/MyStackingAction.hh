#ifndef MyStackingAction_H
#define MyStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"

class MyStackingAction : public G4UserStackingAction
{
public:
  MyStackingAction();
  virtual ~MyStackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);

};


#endif
