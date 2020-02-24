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

class B4PartGeneratorBase : public G4VUserPrimaryGeneratorAction{
public:
    B4PartGeneratorBase(): G4VUserPrimaryGeneratorAction(),energy_(0),xorig_(0),yorig_(0){


    };
    virtual ~B4PartGeneratorBase(){};

    virtual bool isJetGenerator()=0;

    enum particles{
        elec=0,muon,pioncharged,pionneutral,klong,kshort,gamma,
        positron,
        minbias,
        displacedjet,
        particles_size //leave this one
    };

    G4double getEnergy()const{return energy_;}

    virtual std::vector<G4String> generateAvailableParticles()const=0;

    virtual particles getParticle()const=0;

    virtual int isParticle(int i)const=0;


    G4double getX()const{return xorig_;}
    G4double getY()const{return yorig_;}
    G4double getR()const{return std::sqrt(yorig_*yorig_+xorig_*xorig_);}


    G4double getDirX()const{return dirx_;}
    G4double getDirY()const{return diry_;}
    G4double getDirZ()const{return dirz_;}

    G4double getDiffProjTheta()const{return diff_proj_theta_;}
    G4double getDiffProjPhi()const{return diff_proj_phi_;}

    G4double getHowParallel()const{return angle_;}


    static int seedsoffset_;


protected:

    G4double energy_;
    G4double xorig_,yorig_;
    G4double dirx_,diry_,dirz_;
    G4double angle_;
    G4double diff_proj_theta_, diff_proj_phi_;

};



#endif /* DISPLACEDCALO_MINICALO_INCLUDE_B4PARTGENERATORBASE_HH_ */
