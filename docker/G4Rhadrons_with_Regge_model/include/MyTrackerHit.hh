#ifndef MyTrackerHit_h
#define MyTrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Track.hh"

class MyTrackerHit; 
typedef G4THitsCollection<MyTrackerHit> MyTrackerHitsCollection;


class MyTrackerHit : public G4VHit {

public:

  MyTrackerHit();
  ~MyTrackerHit();
  MyTrackerHit(const MyTrackerHit &right);
  // Constructors and Destructor.

  const MyTrackerHit& operator=(const MyTrackerHit &right);
  // Assignment operator.

public:
  
  void Draw();
  void Print();
  
  inline void SetEdep(const G4double de);
  inline void SetPos(const G4ThreeVector tpos);
  inline void AddEdep(const G4double de);
  inline void Set4Mom(const G4LorentzVector p4);
  inline void SetParticleDefinition(G4ParticleDefinition* pdef);

  inline G4double GetEdep() const;
  inline G4ThreeVector GetPosition() const;
  inline G4LorentzVector Get4Mom() const;
  inline G4ParticleDefinition* GetParticleDefinition() const;


private:

  G4double edep;      // Deposited energy.  
  G4ThreeVector position;
  G4LorentzVector fourmomentum;
  G4ParticleDefinition* particle;
};

// Set methods

inline void MyTrackerHit::SetEdep(const G4double de) { 
  edep = de;
}

inline void MyTrackerHit::SetPos(const G4ThreeVector tpos) {
  //  G4cout<<"Treating hit at position: "<< tpos<<G4endl;
  position = tpos;
}

inline void MyTrackerHit::AddEdep(const G4double de) { 
  edep += de;
}

inline void MyTrackerHit::Set4Mom(const G4LorentzVector p4){
  fourmomentum = p4;
}

inline void MyTrackerHit::SetParticleDefinition(G4ParticleDefinition* pdef){
  particle = pdef;
}

// Get methods

inline G4double MyTrackerHit::GetEdep() const { 
  return edep;
}

inline G4ThreeVector MyTrackerHit::GetPosition() const { 
  return position;
}

inline G4LorentzVector MyTrackerHit::Get4Mom() const {
  return fourmomentum;
}

inline G4ParticleDefinition* MyTrackerHit::GetParticleDefinition() const {
  //  G4cout<<"Returning Particle of type: "<<particle->GetParticleType()<<G4endl;
  return particle;
}

#endif


