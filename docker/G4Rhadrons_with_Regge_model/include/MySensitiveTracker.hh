#ifndef MySensitiveTracker_h
#define MySensitiveTracker_h 1

#include "G4VSensitiveDetector.hh"
#include "MyTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class MySensitiveTracker : public G4VSensitiveDetector {

public:

  MySensitiveTracker(G4String name);
  ~MySensitiveTracker();

  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);

  void clear();
  void DrawAll();
  void PrintAll();

private:

  MyTrackerHitsCollection* trackerCollection;

};


#endif
