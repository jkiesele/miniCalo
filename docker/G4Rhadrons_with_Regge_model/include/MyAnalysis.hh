#ifndef MyAnalysis_h 
#define MyAnalysis_h 1

#include "globals.hh"

#include "G4VHitsCollection.hh"
#include "MySensitiveTracker.hh"
#include "MyTrackerHit.hh"



class MyAnalysis {

public:

  ~MyAnalysis();
  // Close histograms and/or ntuples.
  // This method is called only once.
  
  //  void fill(const G4int numHits, const G4double energyDeposited, const G4double distance); 
  void fill(const MyTrackerHitsCollection*); 
  // This method is called at each event to fill the 
  // histograms and/or ntuples.

  void close();

  static MyAnalysis* getInstance();

private:
  
  MyAnalysis();
  
  static MyAnalysis* instance;

  G4int Nevt;
  
};

#endif

