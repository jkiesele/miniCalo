//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExRhadDetectorConstruction.cc,v 1.4 2005/12/14 13:51:24 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExRhadDetectorConstruction.hh"
#include "ExRhadDetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "MySensitiveTracker.hh"
#include "G4SDManager.hh"


using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadDetectorConstruction::ExRhadDetectorConstruction()
:solidWorld(0),logicWorld(0),physiWorld(0),
 Box1 (0), logicBox1 (0),  physBox1 (0),
 Box2 (0), logicBox2 (0),  physBox2 (0),
 Box3 (0), logicBox3 (0),  physBox3 (0),
 Box4 (0), logicBox4 (0),  physBox4 (0),
 magField(0)
{
  // materials
  DefineMaterials();

  //Setting defaults values for detector geometry

  Box1Size = 30*cm;
  Box2Size = 1*m;
  Box3Size = 30*cm;
  Box4Size = 40*cm;
  DetectorSpacing = 10*cm;
  WorldSize = Box1Size+Box2Size+Box3Size+Box4Size+3*DetectorSpacing/2;
  WorldWidth = 1*m;

  Box1Pos   = -WorldSize+Box1Size;
  Box2Pos = Box1Pos + Box1Size + Box2Size+DetectorSpacing;
  Box3Pos   = Box2Pos + Box2Size + Box3Size+DetectorSpacing;
  Box4Pos = Box3Pos + Box3Size + Box4Size+DetectorSpacing;

  // create commands for interactive definition of the calorimeter
  detectorMessenger = new ExRhadDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExRhadDetectorConstruction::~ExRhadDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExRhadDetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExRhadDetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  
G4int iz, n;                 //iz=number of protons  in an isotope; 
                             // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element("enriched Uranium",symbol="U",ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//

/*
new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
*/

//
// define a material from elements.   case 1: chemical molecule
//

G4Material* H2O = 
new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

Sci = 
new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

myWater = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
myWater->AddElement(H, natoms=2);
myWater->AddElement(O, natoms=1);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
myWater->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

myIron = new G4Material("Iron"     , z=26., a= 55.85*g/mole , density= 7.870*g/cm3);
myLead = new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
 myCarbon = new G4Material("Carbon", density = 2267*kg/m3, ncomponents=1);
 myCarbon->AddElement(C, natoms=1);

 myHydrogen = new G4Material("Hydrogen", density = 0.0899*kg/m3, ncomponents=1);
 myHydrogen->AddElement(H, natoms=1);

G4Material* Myl = 
new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
Myl->AddElement(C, natoms=10);
Myl->AddElement(H, natoms= 8);
Myl->AddElement(O, natoms= 4);

G4Material* SiO2 = 
new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

G4Material* Air = 
new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

 myAir = Air;

//
// define a material from elements and/or others materials (mixture of mixtures)
//

G4Material* Aerog = 
new G4Material("Aerogel", density= 0.200*g/cm3, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
Aerog->AddElement (C   , fractionmass= 0.1*perCent);

//
// examples of gas in non STP conditions
//

G4Material* CO2 = 
new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                              kStateGas, 325.*kelvin, 50.*atmosphere);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
G4Material* steam = 
new G4Material("WaterSteam", density= 0.3*mg/cm3, ncomponents=1,
                             kStateGas, 500.*kelvin, 2.*atmosphere);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

G4Material* Vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 3.e-18*pascal, 2.73*kelvin);

G4Material* beam = 
new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                       kStateGas, STP_Temperature, 2.e-2*bar);
beam->AddMaterial(Air, fractionmass=1.);

G4cout << "Materials:"<<G4endl;
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the World
defaultMaterial  = Vacuum;


 Box1Material = Sci;
 Box2Material = myIron;
 Box3Material = myAir;
 Box4Material = myIron;
 
 G4cout <<"Materials: "<<Box1Material<<", "<<Box2Material<<", "<<Box3Material<<", "<<Box4Material<<G4endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExRhadDetectorConstruction::ConstructCalorimeter()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //     
  // World
  //

  solidWorld = new G4Box("World",				//its name
			 WorldSize+10*cm,WorldWidth+10*cm,WorldWidth+10*cm);	//its size
                          
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume				 
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Placing our "detector elements"
  //
  
  Box1 = new G4Box("Box1Solid",Box1Size,WorldWidth,WorldWidth);
  logicBox1 = new G4LogicalVolume(Box1,Box1Material,"logicBox1");
  physBox1 =  new G4PVPlacement(0,
				  G4ThreeVector(Box1Pos,0,0),
				  logicBox1,
				  "PhysBox1",
				  logicWorld,
				  false,
				  0);

  Box2 = new G4Box("Box2Solid",Box2Size,WorldWidth,WorldWidth);
  logicBox2 = new G4LogicalVolume(Box2,Box2Material,"logicBox2");
  physBox2 =  new G4PVPlacement(0,
				  G4ThreeVector(Box2Pos,0,0),
				  logicBox2,
				  "PhysBox2",
				  logicWorld,
				  false,
				  0);

  Box3 = new G4Box("Box3Solid",Box3Size,WorldWidth,WorldWidth);
  logicBox3 = new G4LogicalVolume(Box3,Box3Material,"logicBox3");
  physBox3 =  new G4PVPlacement(0,
				  G4ThreeVector(Box3Pos,0,0),
				  logicBox3,
				  "PhysBox3",
				  logicWorld,
				  false,
				  0);

  Box4 = new G4Box("Box4Solid",Box4Size,WorldWidth,WorldWidth);
  logicBox4 = new G4LogicalVolume(Box4,Box4Material,"logicBox4");
  physBox4 =  new G4PVPlacement(0,
				  G4ThreeVector(Box4Pos,0,0),
				  logicBox4,
				  "PhysBox4",
				  logicWorld,
				  false,
				  0);
  
  //Sensitivity

  G4bool warning;
  if (!G4SDManager::GetSDMpointer()->FindSensitiveDetector("TrackerBox",warning))//Only add the tracker if it isn't already there...
    {
      Tracker = new MySensitiveTracker( "TrackerBox" );
      G4SDManager::GetSDMpointer()->AddNewDetector( Tracker  );
    }

  logicBox1->SetSensitiveDetector( Tracker );
  logicBox2->SetSensitiveDetector( Tracker );
  logicBox3->SetSensitiveDetector( Tracker );
  logicBox4->SetSensitiveDetector( Tracker );

  //  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  logicWorld->SetVisAttributes (G4Colour(0.25,0.25,0.25));

  logicBox1->SetVisAttributes (G4Colour(0,1,1));
  logicBox2->SetVisAttributes (G4Colour(0.5,0.5,0.5));
  logicBox3->SetVisAttributes (G4Colour(1,1,0));
  logicBox4->SetVisAttributes (G4Colour(0.5,0.5,0.5));

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);

  //
  //always return the physical World
  //
  return physiWorld;
}


#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void ExRhadDetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field

  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,fieldValue,0.));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void ExRhadDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
