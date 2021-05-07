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
  fCheckOverlaps(true),
  m_vacuum(0),
  m_pb(0),
  m_pbtungsten(0),
  m_silicon(0),
  m_cu(0)

{


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


G4double computeTheta(const G4double& eta )
        { return 2. * atan( exp( -eta ) ); }

double etaToR(const G4double& eta, const G4double& z){
    return z * exp(-eta);
}
G4Cons * createCons(
        G4String name,
        G4double start_eta,
        G4double eta_width,
        G4double start_z,
        G4double z_length,
        G4double start_phi,
        G4double end_phi){



  //  G4cout << "created Cons "<< name <<"with size (eta, deta, phi,endphi,z,dz): "
  //          << start_eta <<", "<< eta_width <<", "
  //          << start_phi << ", "<<end_phi<<", "
  //          << start_z << ", "<< z_length  << G4endl;
  //
  //  G4cout << "(pRmin1,pRmax1) " << etaToR(start_eta+eta_width ,start_z) <<", " << etaToR(start_eta ,start_z)<< G4endl;
  //  G4cout << "(pRmin2,pRmax2) " << etaToR(start_eta+eta_width ,start_z+z_length) <<", " << etaToR(start_eta ,start_z+z_length)<< G4endl;


    return new G4Cons(
            name,
            etaToR(start_eta+eta_width ,start_z),
            etaToR(start_eta ,start_z),
            etaToR(start_eta+eta_width ,start_z+z_length),
            etaToR(start_eta ,start_z+z_length),
            z_length/2.,
            start_phi,
            end_phi);
}

G4VPhysicalVolume * B4DetectorConstruction::createTubs(

        const G4String& name,
        G4ThreeVector position,
        G4double inner,
        G4double outer,
        G4double zhalf,
        G4Material* material,
        G4LogicalVolume* mother,
        bool active,
        int layernum,
        int copynum,
        bool halfopen
){

    ///
    /*
     * G4VPhysicalVolume * vol,
    G4double eta,
    G4double phi,
    G4double null,
    G4double posx,
    G4double posy,
    G4double posz,
    int layer, G4int copyno
     */

    G4double wedge=2*pi;
    if(halfopen)
        wedge=pi;


    auto S   = new G4Tubs(name+"_S",           // its name
            inner,
            outer,
            zhalf,
            0,
            wedge); // its size

    auto LV  = new G4LogicalVolume(
            S,
            material,
            name+"_LV");
    CLHEP::HepRotationZ zrot;
    zrot.setDelta(pi/2.);
    auto rot= new G4RotationMatrix (zrot);
    auto PV = new G4PVPlacement(
            rot,                //  rotation
            position, // its position
            LV,       // its logical volume
            name+"_P",           // its name
            mother,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

    if(!active)
        return PV;


    activecells_.push_back(
            sensorContainer(
                    PV,
                    outer-inner,//sensor size   // G4double dimxy
                    0,                              // G4double dimz,
                    0,       // G4double area,
                    inner,                   // G4double posx,
                    position.y(),                   // G4double posy,
                    position.z(),                   // G4double posz,
                    layernum,
                    copynum //copyno
            ));

}

/*
 * everything already in world coordinates. simple enough to do that
 */
G4VPhysicalVolume* B4DetectorConstruction::createCMS(
        G4ThreeVector position,
        G4LogicalVolume* worldLV){

    G4double ECalInner=1*m;
    G4double ECalOuter=1.2*m;
    G4double HCalInner=1.3*m;
    G4double HCalOuter=2.3*m;
    G4double SolenoidInner=2.4*m;
    G4double SolenoidOuter=3*m;
    G4double MuonInner=3.1*m;
    G4double MuonOuter=5*m;
    G4double InnerZ= 6.5*m;


    bool consideractive=false;
#ifdef USECMSACTIVE
    consideractive=true;
#endif
    createTubs("ECal",position,ECalInner,ECalOuter,InnerZ,m_pbtungsten,worldLV,consideractive,0,0);
    createTubs("HCal",position,HCalInner,HCalOuter,InnerZ,m_brass,worldLV,     consideractive,1,0);
    createTubs("Sol",position,SolenoidInner,SolenoidOuter,InnerZ,m_brass,worldLV,consideractive,2,0);
    createTubs("Muons",position,MuonInner,MuonOuter,InnerZ,m_iron,worldLV,consideractive,3,0);




return 0;
}

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



G4VPhysicalVolume* B4DetectorConstruction::createLayer(
        G4String name,
        G4Material* m_scintillator,
        G4Material* m_rods,
        G4ThreeVector pos,
        G4ThreeVector dxyz,
        G4ThreeVector roddxyz,
        G4ThreeVector roddistxyz, //distance between rods, not from centre to centre.
        G4LogicalVolume* mother,
        G4int layerno){

    roddxyz.setY(dxyz.y());

    //directly place rods
  //  G4LogicalVolume* layerLV=0;
   // auto layerPV = createBox(name+"_scintillator",pos,dxyz,m_scintillator,mother,layerLV);

    int nrodsz = dxyz.z()/(roddxyz.z()+roddistxyz.z());
    int nrodsx = 2;//fixed, two layers of offset rods

    auto startpos = pos - dxyz/2.; //
    startpos += roddxyz/2.+roddistxyz/2.;
    auto currentpos=startpos;

    for(int ix=0;ix<nrodsz;ix++){
        G4String sx="";
        sx+=ix;
        for(int i=0;i<nrodsz;i++){
            G4String s=sx+"_";
            s+=std::to_string(i);
            G4LogicalVolume* rodLV=0;
            auto rodPV = createBox(name+"_rod_"+s,currentpos,roddxyz,m_rods,mother,rodLV);


            //just vis attributes

            auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(.7,0.65,0.26));
            simpleBoxVisAtt->SetVisibility(true);
            simpleBoxVisAtt->SetForceSolid(true);
            rodPV->GetLogicalVolume()->SetVisAttributes(simpleBoxVisAtt);

            activecells_.push_back(
                            sensorContainer(
                                    rodPV,
                                    0,//sensor size   // G4double dimxy
                                    0,                              // G4double dimz,
                                    0,       // G4double area,
                                    currentpos.x(),                   // G4double posx,
                                    currentpos.y(),                   // G4double posy,
                                    currentpos.z(),                   // G4double posz,
                                    layerno,
                                    0 //copyno
                            ));

            G4cout << "added active cell " << name+"_rod_"+s << " layer "<< layerno << G4endl;

            currentpos += G4ThreeVector(0,0,roddxyz.z()+roddistxyz.z());

            if(ix){
                if(i==nrodsz-2)
                    break;//one less in second row
            }
        }
        currentpos=startpos + G4ThreeVector(roddxyz.x()+roddistxyz.x(), 0, (roddxyz.z()+roddistxyz.z())/2.);
    }
    //now the second row




    return 0;//layerPV;

}

G4VPhysicalVolume* B4DetectorConstruction::createBottle(
        G4ThreeVector position,
        G4LogicalVolume * worldLV){


    auto currentpos = position;
    int startat=0;
#ifdef USECMSACTIVE
    startat=4;
#endif


    for (int i=startat;i<NACTIVELAYERS;i++){
        G4String s="";
        s+=std::to_string(i);
        createLayer("Layer_"+s,
                m_vacuum,
                m_brass,
                currentpos,
                G4ThreeVector(2.*cm, 6.*m, 6*m),//layer size always two rods in x
                G4ThreeVector(1.0*cm, 6*m, 20*cm),//rod size //1cm
                G4ThreeVector(0.*cm, 0, 0.*cm),//rod distance
                worldLV,
                i);


        currentpos += G4ThreeVector(3.*cm,0,0);

    }

    return 0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
	// Lead material defined using NIST Manager
	auto nistManager = G4NistManager::Instance();
	nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_Cu");
	nistManager->FindOrBuildMaterial("G4_PbWO4");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_BRASS"); //lambda_0 between 16 and 19 cm, radiation length between 1.5 and 2 cm
    nistManager->FindOrBuildMaterial("G4_Fe");
 //   nistManager->FindOrBuildMaterial("G4_lAr");

	//nistManager->ListMaterials("all");
	//exit(1);

	// Liquid argon material
	G4double a;  // mass of a mole;
	G4double z;  // z=mean number of protons;
	G4double density;
	new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
	// The argon by NIST Manager is a gas with a different density

    G4Element *Ar = new G4Element("Argon", "Ar", z= 18., a= 39.948  *g/mole);
    G4Material *LAr = new G4Material("LAr",      1.40*g/cm3, 1,  kStateLiquid, 87* kelvin,   1. * atmosphere);
    LAr->AddElement(Ar, 1);

	// Vacuum
	new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);

	// Print materials
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;
// *m_brass, *m_iron, *m_lar

	// Get materials
	m_vacuum = G4Material::GetMaterial("Galactic"); // "Galactic"
	m_pb = G4Material::GetMaterial("G4_Pb");
	m_pbtungsten = G4Material::GetMaterial("G4_PbWO4");
	m_silicon = G4Material::GetMaterial("G4_Si");
	m_cu = G4Material::GetMaterial("G4_Cu");
	m_brass = G4Material::GetMaterial("G4_BRASS");
	m_iron = G4Material::GetMaterial("G4_Fe");
	m_lar = G4Material::GetMaterial("LAr");

	if ( ! m_vacuum || ! m_pb || ! m_pbtungsten || !m_silicon
	        || ! m_brass || ! m_iron || !m_lar) {
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
    G4double worldSizeXY = 20. * m;
	G4double worldSizeZ  = 20. * m;



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

	//
	// Calorimeter
	//
	createCMS(G4ThreeVector(0,0,0),worldLV);
	createBottle(G4ThreeVector(6*m,0,0),worldLV);



	G4cout << "created in total "<< activecells_.size()<<" sensors" <<G4endl;

	//
	// Visualization attributes
	//


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());


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
