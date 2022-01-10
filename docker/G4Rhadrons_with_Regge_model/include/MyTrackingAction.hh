#ifndef USERTRACKINGACTION_HH
#define USERTRACKINGACTION_HH

#include "G4UserTrackingAction.hh"
#include "globals.hh"


class MyTrackingAction : public G4UserTrackingAction {

public:
  MyTrackingAction();
  virtual ~MyTrackingAction() { }
  void PreUserTrackingAction(const G4Track*) { }
  virtual void PostUserTrackingAction(const G4Track*);

private:
  G4int counter;
};

#endif
