#include "MyTrackerHit.hh"
#include "G4ios.hh"


MyTrackerHit::MyTrackerHit() : edep(0.0) {}


MyTrackerHit::~MyTrackerHit() {}


MyTrackerHit::MyTrackerHit(const MyTrackerHit &right) :
  G4VHit(),  edep(right.edep), position(right.position) 
{
}


const MyTrackerHit& MyTrackerHit::operator=(const MyTrackerHit &right) {
  edep = right.edep;
  position = right.position;
  fourmomentum = right.fourmomentum;
  particle = right.particle;
  return *this;
}


void MyTrackerHit::Draw() {}


void MyTrackerHit::Print() {}


