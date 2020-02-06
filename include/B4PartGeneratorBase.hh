/*
 * B4PartGeneratorBase.hh
 *
 *  Created on: 6 Feb 2020
 *      Author: jkiesele
 */

#ifndef DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_
#define DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_

#include "G4VUserPrimaryGeneratorAction.hh"


class B4PartGeneratorBase : public G4VUserPrimaryGeneratorAction{
public:
    B4PartGeneratorBase(): G4VUserPrimaryGeneratorAction(),energy_(0),xorig_(0),yorig_(0){};
    virtual ~B4PartGeneratorBase(){};

    virtual bool isJetGenerator()=0;

    enum particles{
        elec=0,muon,pioncharged,pionneutral,klong,kshort,gamma,
        positron,
        quark,gluon,
        particles_size //leave this one
    };

    G4double getEnergy()const{return energy_;}

    virtual std::vector<G4String> generateAvailableParticles()const=0;

    virtual particles getParticle()const=0;

    virtual int isParticle(int i)const=0;


    G4double getX()const{return xorig_;}
    G4double getY()const{return yorig_;}
    G4double getR()const{return std::sqrt(yorig_*yorig_+xorig_*xorig_);}


protected:

    G4double energy_;
    G4double xorig_,yorig_;
};



#endif /* DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_ */
