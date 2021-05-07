/*
 * B4PartGeneratorBase.hh
 *
 *  Created on: 6 Feb 2020
 *      Author: jkiesele
 */

#ifndef DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_
#define DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4INCLRandomSeedVector.hh"

class G4ParticleGun;

class B4PartGeneratorBase : public G4VUserPrimaryGeneratorAction{
public:
    B4PartGeneratorBase(): G4VUserPrimaryGeneratorAction(),energy_(0){}

    virtual ~B4PartGeneratorBase(){};

    virtual bool isJetGenerator()=0;

    enum particles{
        particles_size //leave this one
    };

    G4double getEnergy()const{return energy_;}

    virtual std::vector<G4String> generateAvailableParticles()const=0;

    virtual particles getParticle()const=0;

    virtual int isParticle(int i)const=0;

    virtual  G4ParticleGun* getGun(){return 0;}

    static int seedsoffset_;

    static G4String particle;
    static G4double beta;

    const std::vector<std::pair<G4String,int > > & availParticles()const{
        return availParticles_;
    }

protected:

    std::vector<std::pair<G4String, int>> availParticles_;
    std::vector<int> availPids_;

    G4double energy_;

};



#endif /* DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_ */
