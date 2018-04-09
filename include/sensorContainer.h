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
	sensorContainer():vol_(0),dimxy_(0),area_(0),
	posx_(0),posy_(0),posz_(0),energyscalefactor_(1){

	}
	sensorContainer(G4VPhysicalVolume * vol,
	G4double dimxy,
	G4double area,
	G4double posx,
	G4double posy,
	G4double posz,
	G4double layer):
		vol_(vol),dimxy_(dimxy),area_(area),
		posx_(posx),posy_(posy),posz_(posz),energyscalefactor_(1),
		layer_(layer)
	{}



	const G4double& getArea() const {
		return area_;
	}

	const G4double& getDimxy() const {
		return dimxy_;
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

	const G4double& getLayer() const {
		return layer_;
	}

private:
	G4VPhysicalVolume * vol_;
	G4double dimxy_;
	G4double area_;

	G4double posx_;
	G4double posy_;
	G4double posz_;

	G4double energyscalefactor_;

	G4double layer_;

};


#endif /* B4A_INCLUDE_SENSORCONTAINER_H_ */
