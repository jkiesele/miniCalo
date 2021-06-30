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

#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4UserLimits.hh"
#include "G4PVDivision.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4LogicalVolumeModel.hh"

#include "G4Cons.hh"
#include "G4Tubs.hh"
#include <cmath>

#include <CLHEP/Vector/Rotation.h>
#include "CLHEP/Vector/RotationZ.h"
#include "sensorContainer.h"

#include <cstdlib>
//#include "Math/Vector3D.h"

//#define USEDIVISION

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
  m_vacuum(0),
  m_target(0),
  target_material_str_("G4_Pb")
{

    targetvolume_=0;

    limit_in_calo_time_max_=100.*ms;//this could be more low energy stuff
    limit_in_calo_energy_max_=.01*eV;
    limit_world_time_max_=100.*ms; //either got there or not (30ns should be easily sufficient
    limit_world_energy_max_=.1*eV;

}

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
	// Define materials
	DefineMaterials();

	// Define volumes
	return DefineVolumes();
}

void  B4DetectorConstruction::DefineGeometry(geometry geo){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VPhysicalVolume* B4DetectorConstruction::createBox(
        G4String name,
        G4ThreeVector pos,
        G4ThreeVector dxyz,
        G4Material* m,
        G4LogicalVolume* mother,
        G4LogicalVolume*& LV
){

    auto S = new G4Box(name+"_S",           // its name
            dxyz.x()/2, dxyz.y()/2, dxyz.z()/2); // its size

    LV = new G4LogicalVolume(
            S,           // its solid
            m,  // its material
            name+"_LV");         // its name

    //LV->SetUserLimits(new G4UserLimits(
    //        DBL_MAX, //max step length
    //        DBL_MAX, //max track length
    //        DBL_MAX, //max track time
    //        0.1*eV)); //min track energy


    auto PV = new G4PVPlacement(
            0,                // no rotation
            pos,  // at (0,0,0)
            LV,          // its logical volume
            name+"_P",          // its name
            mother,                // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

    return PV;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
	// Lead material defined using NIST Manager
	auto nistManager = G4NistManager::Instance();
	nistManager->FindOrBuildMaterial(target_material_str_);


	// Vacuum
	new G4Material("Galactic", 1., 1.01*g/mole,universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);


	// Get materials
	m_vacuum = G4Material::GetMaterial("Galactic"); // "Galactic"
	m_target = G4Material::GetMaterial(target_material_str_); // "Galactic"

	if ( ! m_vacuum || ! m_target) {
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
	//calorSizeXY  = 200*cm;






    //auto calorThickness = nofEELayers * layerThicknessEE + nofHB*layerThicknessHB;
    G4double worldSizeXY = 1. * m;
	G4double worldSizeZ  = 1. * m;



	//
	// World
	//
	auto worldS
	= new G4Box("World",           // its name
			worldSizeXY/2, worldSizeXY/2, worldSizeZ); // its size

	auto worldLV
	= new G4LogicalVolume(
			worldS,           // its solid
			m_vacuum,  // its material
			"World");         // its name

	//worldLV->SetUserLimits(new G4UserLimits(
	//        worldSizeZ/10., //max step length
    //        worldSizeZ*50., //max track length
    //        limit_world_time_max_, //max track time
    //        limit_world_energy_max_)); //min track energy


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


	//add the target
	G4LogicalVolume * boxvol=0;
	targetvolume_ = createBox("target",G4ThreeVector(),
	        G4ThreeVector(10.*cm , 10.*cm, 1.*cm),
	        m_target,worldLV,boxvol);


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(.7,0.65,0.26));
	simpleBoxVisAtt->SetVisibility(true);
	simpleBoxVisAtt->SetForceSolid(true);
	targetvolume_->GetLogicalVolume()->SetVisAttributes(simpleBoxVisAtt);

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
	G4ThreeVector fieldValue(0,0,0.);
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(2);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);


	auto* fieldprop = G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
	fieldprop->SetMaxLoopCount(100) ;//check 100, 10 is bad, less is bad, default is 1000, maybe a bit less works, too


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
