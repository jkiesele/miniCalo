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

#include "G4Cons.hh"
#include <cmath>

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
  m_pb(0),
  m_pbtungsten(0),
  m_silicon(0),
  m_cu(0)

{

}

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
	//DefineGeometry(homogenous_ecal_only);
	//DefineGeometry(hcal_only_irregular);
	DefineGeometry(fullendcap);
	// Define materials
	DefineMaterials();

	// Define volumes
	return DefineVolumes();
}

void  B4DetectorConstruction::DefineGeometry(geometry geo){

	if(geo == fullendcap){
		calorThickness=2200*mm;

        nofEELayers = 14;
        Ncalowedges=120;
        nofHB=0;
        int etasegments=30;

        Calo_start_eta=1.5;
        Calo_end_eta=3.0;
        Calo_start_z=300*cm;

		for(int i=0;i<nofEELayers+nofHB;i++){
		    if(i<nofEELayers)
		        layerThicknesses.push_back(34*cm / (float)nofEELayers);//CMS like about 30 X0
		    else
		        layerThicknesses.push_back((calorThickness - 34*cm) / (float)nofHB);

		    layerGranularity.push_back(etasegments);
		}

		layerThicknessEE=26*cm / (float)nofEELayers;
		layerThicknessHB=(calorThickness-nofEELayers*layerThicknessEE)/(float)nofHB; //100*mm;

	}
	else{
	    G4ExceptionDescription msg;
	            msg << "Geometry not supported in this branch";
	            G4Exception("B4DetectorConstruction::DefineGeometry()",
	                    "MyCode0001", FatalException, msg);
	}
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


    G4cout << "created Cons "<< name <<"with size (eta, deta, phi,endphi,z,dz): "
            << start_eta <<", "<< eta_width <<", "
            << start_phi << ", "<<end_phi<<", "
            << start_z << ", "<< z_length  << G4endl;

    G4cout << "(pRmin1,pRmax1) " << etaToR(start_eta+eta_width ,start_z) <<", " << etaToR(start_eta ,start_z)<< G4endl;
    G4cout << "(pRmin2,pRmax2) " << etaToR(start_eta+eta_width ,start_z+z_length) <<", " << etaToR(start_eta ,start_z+z_length)<< G4endl;


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
/*
 * creates a single sandwich tile in a layer
 */
G4VPhysicalVolume* B4DetectorConstruction::createCellWheel(
        G4LogicalVolume* layerLV,
		G4double start_eta,
        G4double eta_width,
		G4double start_z,
        G4double z_length,
        G4int nphi,
        G4int layernum,
        G4int cellnum,
        G4Material* active_m,
        G4Material* abs_m,
        G4double abs_fraction){



    G4double active_z_length=z_length;
    G4double abs_z_length=0;
    if(abs_fraction>0){
        active_z_length *= 1.-abs_fraction;
        abs_z_length = z_length-active_z_length;
    }

    //create the volume

    G4String name = createString(layernum)+"_"+createString(cellnum);

	auto cellS   = createCons("Cell_"+name,
	        start_eta,
	        eta_width,
	        start_z+abs_z_length,
	        active_z_length,
	        0,

	        2.*M_PI/(double)nphi);

	auto cellLV  = new G4LogicalVolume(
	        cellS,           // its solid
	        active_m,  // its material
			"Cell_LV_"+name);         // its name




    G4double maxStep = z_length/200.;
    G4double maxTime = 2.*s;
    G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);
    cellLV->SetUserLimits(stepLimit);


#ifdef USEDIVISION
    G4VPhysicalVolume* activeMaterial
            = new G4PVDivision("Cell_rep_"+name, cellLV,layerLV, kPhi, nphi,0.);
#else
    G4VPhysicalVolume* activeMaterial
        = new G4PVReplica("Cell_rep_"+name, cellLV,
                layerLV, kPhi, nphi, 2.*M_PI/(double)nphi-1e-7);
#endif

    if(abs_fraction>0 && abs_m){
        auto absS   = createCons("Abs_"+name,
                start_eta,
                eta_width,
                start_z,
                abs_z_length,
                0,

                2.*M_PI/(double)nphi);

        auto absLV  = new G4LogicalVolume(
                absS,           // its solid
                abs_m,  // its material
                "Abs_LV_"+name);         // its name

        G4VPhysicalVolume* absMaterial
        = new G4PVReplica("Abs_rep_"+name, absLV,
                layerLV, kPhi, nphi, 2.*M_PI/(double)nphi-1e-7);


    }

    G4int maxcopies = activeMaterial->GetMultiplicity();

    //only consider active material here

    double width_rad = 2.*M_PI/(double)nphi*rad;
    for(G4int i =0;i<maxcopies;i++){
        double phi= width_rad*(double)i; //offset doesn't really matter

        double x = cos(phi)*etaToR(start_eta+eta_width/2. ,start_z+z_length/2.);
        double y = sin(phi)*etaToR(start_eta+eta_width/2. ,start_z+z_length/2.);

        activecells_.push_back(
                sensorContainer(
                        activeMaterial,
                        start_eta+eta_width/2.,//sensor size   // G4double dimxy
                        phi,                              // G4double dimz,
                        0,       // G4double area,
                        x,                   // G4double posx,
                        y,                   // G4double posy,
                        start_z+z_length/2.,                   // G4double posz,
                        layernum,
                        i //copyno
                )
        );
    }
    return activeMaterial;

}

G4VPhysicalVolume* B4DetectorConstruction::createLayer(
        G4LogicalVolume * caloLV,
        G4double start_eta,
        G4double end_eta,
        G4int n_cells_eta,
        G4int n_cells_phi,
        G4double start_z,
        G4double z_length,
        G4ThreeVector position,
		G4String name, int layernumber){


    //create the volume

    //override:
    z_length = 10.382*mm + 0.300*mm ;//absorber plus silicon, this makes 1.85 X0 per layer

    auto layerS   = createCons("Layer_"+name,
            start_eta,
            end_eta-start_eta,
            start_z,
            z_length,
            0,
            2*M_PI);


    auto layerLV  = new G4LogicalVolume(
            layerS,           // its solid
            m_vacuum,  // its material
			"Layer_"+name);         // its name


    auto layerPV = new G4PVPlacement(
                0,                // no rotation
                position+G4ThreeVector(0,0,z_length/2.), // its position
                layerLV,       // its logical volume
                "Layer_"+name,           // its name
                caloLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps

    G4double maxStep = z_length/200.;
    G4double maxTime = 2.*s;
    G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);
    layerLV->SetUserLimits(stepLimit);


	double cell_etawidth = (end_eta-start_eta) / (double)n_cells_eta;
	double cell_phiwidth = 2*M_PI/(double)n_cells_phi *rad;

	//create cells
	int cellno=0;
	for(int ieta=0; ieta<n_cells_eta; ieta++){
	    double startEta = Calo_start_eta + cell_etawidth* (double)ieta;

	    auto ringSvol = createCons("Ring_"+createString(ieta) +"_"+name,
	            startEta,
	            cell_etawidth,
	            start_z,
	            z_length,
	            0,
	            2*M_PI);
	    auto RingLV  = new G4LogicalVolume(
	            ringSvol,           // its solid
	                m_vacuum,  // its material
	                "RingLV_"+createString(ieta) +"_"+name);         // its name
	    RingLV->SetUserLimits(stepLimit);
	    auto RingPV = new G4PVPlacement(
	                    0,                // no rotation
	                    G4ThreeVector(0,0,0), // its position
	                    RingLV,       // its logical volume
	                    "RingPV_"+createString(ieta) +"_"+name,           // its name
	                    layerLV,          // its mother  volume
	                    false,            // no boolean operation
	                    0,                // copy number
	                    fCheckOverlaps);  // checking overlaps

	   // for(int iphi=0;iphi<n_cells_phi;iphi++){
	     //   double starting_angle_rad = cell_phiwidth * (double)iphi;

	        createCellWheel(RingLV,
	                startEta,
	                cell_etawidth,
	                start_z,
	                z_length,
	                n_cells_phi,
	                layernumber,
	                cellno,
	                m_silicon,
	                m_cu,
	                0.9719153717 //absorber fraction
	                );

	        cellno++;
	  //  }
	}



	G4cout << "layer "<<name<<" position="<<position << " " << layerPV->GetTranslation () <<G4endl;

	return layerPV;

}

void B4DetectorConstruction::createCalo(G4LogicalVolume * worldLV,
        G4ThreeVector position,G4String name){


//define the geometries
    /*
     *
    G4double layerThicknessEE,layerThicknessHB;
    std::vector<G4double> layerThicknesses;
    G4int Ncalowedges;
    std::vector<G4int> layerGranularity;
    std::vector<G4int> layerSplitGranularity;
    layerGranularity -> neta
     *
     */


	G4double lastzpos=Calo_start_z;//;
	for(int i=0;i<nofEELayers+nofHB;i++){

	    auto thickness = layerThicknesses.at(i);
	    auto newheels = layerGranularity.at(i);

		G4ThreeVector createatposition=G4ThreeVector(0,0,lastzpos)+position;

		createLayer(worldLV,
		        Calo_start_eta,
		        Calo_end_eta,
		        newheels, //n_cells_eta
		        Ncalowedges, // n_cells_phi,
		        lastzpos, // start_z,
		        thickness, // z_length,
		        G4ThreeVector(0,0,lastzpos),
		        "Layer_"+createString(i),  i);

		lastzpos+=thickness;
	}


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
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;


	// Get materials
	m_vacuum = G4Material::GetMaterial("Galactic");
	m_pb = G4Material::GetMaterial("G4_Pb");
	m_pbtungsten = G4Material::GetMaterial("G4_PbWO4");
	m_silicon = G4Material::GetMaterial("G4_Si");
	m_cu = G4Material::GetMaterial("G4_Cu");

	if ( ! m_vacuum || ! m_pb || ! m_pbtungsten || !m_silicon ) {
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
    G4double worldSizeXY = 2. * m;
	G4double worldSizeZ  = 8. * m;



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


	G4double maxStep = worldSizeZ/100.;
	G4double maxTime = 20.*s;
	G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);
	worldLV->SetUserLimits(stepLimit);


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
	createCalo(worldLV,G4ThreeVector(0,0,Calo_start_z),"");



	G4cout << "created in total "<< activecells_.size()<<" sensors" <<G4endl;

	//
	// Visualization attributes
	//


	worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

	double g=0;
	for(auto& v: activecells_){
	    auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,g,.0));
	    simpleBoxVisAtt->SetVisibility(true);
		v.getVol()->GetLogicalVolume()->SetVisAttributes(simpleBoxVisAtt);
		g+= 1./(double)activecells_.size();
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
	G4ThreeVector fieldValue(0,0,0);//1.);
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(2);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);


	auto* fieldprop = G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
	fieldprop->SetMaxLoopCount(10) ;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
