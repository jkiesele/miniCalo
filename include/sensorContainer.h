/*
 * sensorContainer.h
 *
 *  Created on: 4 Apr 2018
 *      Author: jkiesele
 */

#ifndef B4A_INCLUDE_SENSORCONTAINER_H_
#define B4A_INCLUDE_SENSORCONTAINER_H_

#include "G4VPhysicalVolume.hh"


class sensorContainer{
public:

	sensorContainer(G4VPhysicalVolume * vol,
	G4double eta,
	G4double phi,
	G4double null,
	G4double posx,
	G4double posy,
	G4double posz,
	int layer, G4int copyno):
		vol_(vol),eta_(eta),phi_(phi),area_(null),
		posx_(posx),posy_(posy),posz_(posz),energyscalefactor_(1),
		layer_(layer),copyno_(copyno)
	{
		global_detid_=global_detid_counter_++;
	}



	const G4double& getArea() const {
		return area_;
	}

	const G4double& getEta() const {
		return eta_;
	}

	const G4double& getPhi() const {
		return phi_;
	}

	const G4VPhysicalVolume* getVol() const {
		return vol_;
	}

	const G4double& getPosx() const {
		return posx_;
	}

	const G4double& getPosy() const {
		return posy_;
	}

	const G4double& getPosz() const {
		return posz_;
	}

	G4double getEnergyscalefactor() const {
		return energyscalefactor_;
	}

	void setEnergyscalefactor(G4double energyscalefactor) {
		energyscalefactor_ = energyscalefactor;
	}

	const int& getLayer() const {
		return layer_;
	}

	G4int getCopyNo()const{
	    return copyno_;
	}

	const int& getGlobalDetID()const{
		return global_detid_;
	}

private:
	sensorContainer(){
			global_detid_=global_detid_counter_++;
		}

	G4VPhysicalVolume * vol_;
	G4double eta_;
	G4double phi_;
	G4double area_;

	G4double posx_;
	G4double posy_;
	G4double posz_;

	G4double energyscalefactor_;

	int layer_;

	int global_detid_;

	G4int  copyno_;

	static int global_detid_counter_;

};


#endif /* B4A_INCLUDE_SENSORCONTAINER_H_ */
