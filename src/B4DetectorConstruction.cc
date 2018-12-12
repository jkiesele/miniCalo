//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "sensorContainer.h"

#include <cstdlib>

static G4double epsilon=0.0*mm;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<class T>
static G4String createString(const T& i){
	std::stringstream ss;
	ss << i;
	std::string number=ss.str();
	return number;
}


B4DetectorConstruction::B4DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(false),
  defaultMaterial(0),
  absorberMaterial(0),
  gapMaterial(0)

{

}

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
	DefineGeometry(homogenous_ecal_only);
	// Define materials
	DefineMaterials();

	// Define volumes
	return DefineVolumes();
}

void  B4DetectorConstruction::DefineGeometry(geometry g){

	if(g == standard){
		calorThickness=2000*mm;

	    layerGranularity.clear();
	    layerSplitGranularity.clear();
		nofEELayers = 10;
		nofHB=15;
		for(int i=0;i<nofEELayers+nofHB;i++){
			G4double granularity=8;
			if(i>3)
				granularity=6;
			if(i>7){
				granularity=4;
			}
			if(i>15)
				granularity=2;
			layerGranularity.push_back(granularity);
			if(granularity<2)
				layerSplitGranularity.push_back(0);
			else
				layerSplitGranularity.push_back(2*granularity);
		}

		layerThicknessEE=15*mm;
		layerThicknessHB=(calorThickness-nofEELayers*layerThicknessEE)/(float)nofHB; //100*mm;

	}
	else if(g == homogenous || g == homogenous_ecal_only){
		if (g == homogenous)
			calorThickness=2000*mm;
		else
			calorThickness=300*mm;

	    layerGranularity.clear();
	    layerSplitGranularity.clear();
		nofEELayers = 10;
		nofHB=20;
		for(int i=0;i<nofEELayers+nofHB;i++){
			G4double granularity=10;
			if(i==27)
				granularity=8;
			if(i==28)
				granularity=2;
			if(i==29)
				granularity=1;
			layerGranularity.push_back(granularity);

			layerSplitGranularity.push_back(-granularity/2);
		}

		layerThicknessEE=80*mm * calorThickness/2000*mm;
		layerThicknessHB=(calorThickness-nofEELayers*layerThicknessEE)/(float)nofHB; //100*mm;

	}
	else if(g == ecal_only){
		calorThickness=200*mm;

		layerGranularity.clear();
		layerSplitGranularity.clear();
		nofEELayers = 10;
		nofHB=0;
		for(int i=0;i<nofEELayers+nofHB;i++){
			G4double granularity=10;

			layerGranularity.push_back(granularity);

			layerSplitGranularity.push_back(int(granularity/2));
		}

		layerThicknessEE=20*mm;
		layerThicknessHB=(calorThickness-nofEELayers*layerThicknessEE)/(float)nofHB; //100*mm;

	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


/*
 * creates a single sandwich tile in a layer
 */
G4VPhysicalVolume* B4DetectorConstruction::createSandwich(G4LogicalVolume* layerLV,
		G4double dx,
		G4double dy,
		G4double dz,
		G4ThreeVector position,
		G4String name,
		G4double absorberfraction,
		G4VPhysicalVolume*& absorber){

	auto absdz=absorberfraction*dz;
	auto gapdz=(1-absorberfraction)*dz;

	//auto layerThickness = absoThickness + gapThickness;

	auto sandwichS   = new G4Box("Sandwich_"+name,           // its name
			dx/2-epsilon, dy/2-epsilon, dz/2-epsilon); // its size

	auto sandwichLV  = new G4LogicalVolume(
			sandwichS,           // its solid
			defaultMaterial,  // its material
			"Sandwich_"+name);         // its name



	//
	// Absorber
	//
	auto absorberS
	= new G4Box("Abso_"+name,            // its name
			dx/2-2*epsilon, dy/2-2*epsilon, absdz/2-2*epsilon); // its size

	auto absorberLV
	= new G4LogicalVolume(
			absorberS,        // its solid
			absorberMaterial, // its material
			"Abso_"+name);          // its name

	absorber
	= new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0., 0., -absdz/2), // its position
			absorberLV,       // its logical volume
			"Abso_"+name,           // its name
			sandwichLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps


	//
	// Gap
	//
	auto gapS
	= new G4Box("Gap_"+name,             // its name
			dx/2-2*epsilon, dy/2-2*epsilon, gapdz/2-2*epsilon); // its size

	auto gapLV
	= new G4LogicalVolume(
			gapS,             // its solid
			gapMaterial,      // its material
			"Gap_"+name);           // its name

	auto activeMaterial
	= new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0., 0., gapdz/2), // its position
			gapLV,            // its logical volume
			"Gap_"+name,            // its name
			sandwichLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	//place the sandwich

	//auto sandwichPV =
			new G4PVPlacement(
				0,                // no rotation
				position, // its position
				sandwichLV,       // its logical volume
				"Sandwich_"+name,           // its name
				layerLV,          // its mother  volume
				false,            // no boolean operation
				0,                // copy number
				fCheckOverlaps);  // checking overlaps

	return activeMaterial;


}

G4VPhysicalVolume* B4DetectorConstruction::createLayer(G4LogicalVolume * caloLV,
		G4double thickness,
		G4int granularity, G4double absfraction,G4ThreeVector position,
		G4String name, int layernumber, G4double calibration, G4int    nsmallsensorsrow){



	auto layerS   = new G4Box("Layer_"+name,           // its name
			calorSizeXY/2, calorSizeXY/2, thickness/2); // its size

	auto layerLV  = new G4LogicalVolume(
			layerS,           // its solid
			defaultMaterial,  // its material
			"Layer_"+name);         // its name

	auto layerPV = new G4PVPlacement(
			0,                // no rotation
			position, // its position
			layerLV,       // its logical volume
			"Layer_"+name,           // its name
			caloLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps


	G4double coarsedivider=(G4double)granularity;
	G4double largesensordxy=calorSizeXY/coarsedivider;
	G4double smallsensordxy=largesensordxy/4;
	if(nsmallsensorsrow<0){
    	nsmallsensorsrow=granularity/2;
    	smallsensordxy=largesensordxy;
	}
	if(granularity<2)
		nsmallsensorsrow=0;


	//divide into 4 areas:
	// LG  HG
	// LG  LG
	//
	// LG: low granularity
	// HG: high granularity


	auto placeSensors = [] (
			G4ThreeVector startcorner,
			bool small,
			G4double sensorsize,
			G4double Thickness,
			int gran,
			G4ThreeVector pos,
			G4String lname,
			std::vector<sensorContainer >* acells,
			G4LogicalVolume* layerlogV,
			B4DetectorConstruction* drec,
			G4ThreeVector patentpos,
			G4double absfractio, int laynum, G4double calib) {
		for(int xi=0;xi<gran;xi++){
			G4double posx=startcorner.x()+sensorsize/2+sensorsize*(G4double)xi;
			for(int yi=0;yi<gran;yi++){
				G4double posy=startcorner.y()+sensorsize/2+sensorsize*(G4double)yi;
				if(!small && posy>pos.y() && posx > pos.x()){
					continue; //here are the small sensors
				}
				auto sandwichposition=G4ThreeVector(posx,posy,pos.z());


				G4VPhysicalVolume * absorber=0;
				auto activesensor=drec->createSandwich(layerlogV,sensorsize,sensorsize,
						Thickness,sandwichposition,
						lname+"_sensor_"+createString(xi)+"_"+createString(yi),
						absfractio,absorber);

				sensorContainer sensordesc(activesensor,
						sensorsize,Thickness,sensorsize*sensorsize,
						patentpos.x()+posx,
						patentpos.y()+posy,
						patentpos.z(),laynum,absorber);
				G4cout << "created sensor with ID "<< sensordesc.getGlobalDetID() << G4endl;
				sensordesc.setEnergyscalefactor(calib);
				acells->push_back(sensordesc);
			}

		}
	};

	G4ThreeVector lowerleftcorner=G4ThreeVector(
			0-calorSizeXY/2,
			0-calorSizeXY/2,
			0);


	//place LG sensors:
	if(nsmallsensorsrow>0){
		placeSensors(lowerleftcorner, false,largesensordxy,thickness,
				granularity,G4ThreeVector(0,0,0),name,&activecells_,layerLV,this,
				position,absfraction,layernumber,calibration);
		placeSensors(G4ThreeVector(0,0,0), true,smallsensordxy,thickness,
				nsmallsensorsrow,G4ThreeVector(0,0,0),name,&activecells_,layerLV,
				this,position,absfraction,layernumber,calibration);
	}
	else{

		placeSensors(lowerleftcorner, true,largesensordxy,thickness,
				granularity,G4ThreeVector(0,0,0),name,&activecells_,layerLV,this,
				position,absfraction,layernumber,calibration);
	}
	G4cout << "layer position="<<position <<G4endl;

	return layerPV;

}

void B4DetectorConstruction::createCalo(G4LogicalVolume * caloLV,G4ThreeVector position,G4String name){

	G4double absorberFractionEE=0.0001;
	G4double absorberFractionHB=absorberFractionEE;
	G4double calibrationEE=1;
	G4double calibrationHB=1;

//define the geometries


	G4double lastzpos=-layerThicknessEE;
	for(int i=0;i<nofEELayers+nofHB;i++){
		int granularity=layerGranularity.at(i);
		int splitgranularity = layerSplitGranularity.at(i);
		G4double absfraction=absorberFractionEE;
		G4double thickness=layerThicknessEE;
		G4double calibration=calibrationEE;

		if(i>=nofEELayers){
			absfraction=absorberFractionHB;
			thickness=layerThicknessHB;
			calibration=calibrationHB;
		}
		G4ThreeVector createatposition=G4ThreeVector(0,0,lastzpos+thickness)+position;
		createLayer(
				caloLV,thickness,
				granularity,
				absfraction,
				createatposition,
				name+"layer"+createString(i),i,1,splitgranularity);//calibration);
		G4cout << "created layer "<<  i<<" at "<< createatposition << G4endl;
		lastzpos+=thickness;
	}

	G4cout << "created " << activecells_.size() << " sensors"<<std::endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
	// Lead material defined using NIST Manager
	auto nistManager = G4NistManager::Instance();
	nistManager->FindOrBuildMaterial("G4_Pb");
	nistManager->FindOrBuildMaterial("G4_PbWO4");

	//nistManager->ListMaterials("all");

	// Liquid argon material
	G4double a;  // mass of a mole;
	G4double z;  // z=mean number of protons;
	G4double density;
	new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
	// The argon by NIST Manager is a gas with a different density

	// Vacuum
	new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;


	// Get materials
	defaultMaterial = G4Material::GetMaterial("Galactic");
	absorberMaterial = G4Material::GetMaterial("G4_Pb");
	gapMaterial = G4Material::GetMaterial("G4_PbWO4");

	if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
		G4ExceptionDescription msg;
		msg << "Cannot retrieve materials already defined.";
		G4Exception("B4DetectorConstruction::DefineVolumes()",
				"MyCode0001", FatalException, msg);
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
	// Geometry parameters
	calorSizeXY  = 30.*cm;






    //auto calorThickness = nofEELayers * layerThicknessEE + nofHB*layerThicknessHB;
	auto worldSizeXY = 1.2 * calorSizeXY;
	G4double worldSizeZ  = 11 * m;



	//
	// World
	//
	auto worldS
	= new G4Box("World",           // its name
			worldSizeXY/2, worldSizeXY/2, worldSizeZ); // its size

	auto worldLV
	= new G4LogicalVolume(
			worldS,           // its solid
			defaultMaterial,  // its material
			"World");         // its name

	auto worldPV
	= new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(),  // at (0,0,0)
			worldLV,          // its logical volume
			"World",          // its name
			0,                // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	//
	// Calorimeter
	//
	createCalo(worldLV,G4ThreeVector(0,0,0),"");



	G4cout << "created in total "<< activecells_.size()<<" sensors" <<G4endl;

	//
	// Visualization attributes
	//
	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,.0,.0));
	simpleBoxVisAtt->SetVisibility(true);
	for(auto& v: activecells_){
		v.getVol()->GetLogicalVolume()->SetVisAttributes(simpleBoxVisAtt);
	}
	//
	// Always return the physical World
	//
	return worldPV;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if
	// the field value is not zero.
	G4ThreeVector fieldValue;
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
