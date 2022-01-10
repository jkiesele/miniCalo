//
// ********************************************************************
// * DISCLAMER                                                       *
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
// $Id: ExRhadDetectorConstruction.hh,v 1.3 2005/12/14 13:51:24 mackepr Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExRhadDetectorConstruction_h
#define ExRhadDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UniformMagField;
class ExRhadDetectorMessenger;
class MySensitiveTracker;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExRhadDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  ExRhadDetectorConstruction();
  ~ExRhadDetectorConstruction();
  
public:
  inline void SetWorldSize(G4double);

  inline void SetBox1Size(G4double);
  inline void SetBox2Size(G4double);
  inline void SetBox3Size(G4double);
  inline void SetBox4Size(G4double);

  inline void SetBox1Pos(G4double);
  inline void SetBox2Pos(G4double);
  inline void SetBox3Pos(G4double);
  inline void SetBox4Pos(G4double);

  inline void SetDetectorSpacing(G4double);

  inline void SetBox1Material(G4String);
  inline void SetBox2Material(G4String);
  inline void SetBox3Material(G4String);
  inline void SetBox4Material(G4String);

  void SetMagField(G4double);
     
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
     
public:
  
  G4double GetWorldSizeX()           {return WorldSize;}; 
  G4double GetWorldSizeYZ()          {return WorldWidth;};
     
private:
  G4Material* defaultMaterial;

  //My stuff
  G4Material* Sci;
  
  G4Material* myIron;
  G4Material* myCarbon;
  G4Material* myHydrogen;
  G4Material* myLead;
  G4Material* myWater;
  G4Material* myAir;

  G4Material* Box1Material;
  G4Material* Box2Material;
  G4Material* Box3Material;
  G4Material* Box4Material;

  G4double WorldSize;
  G4double WorldWidth;

  G4double Box1Size;
  G4double Box2Size;
  G4double Box3Size;
  G4double Box4Size;

  G4double DetectorSpacing;

  G4double Box1Pos;
  G4double Box2Pos;
  G4double Box3Pos;
  G4double Box4Pos;

  MySensitiveTracker* Tracker;

  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World

  G4Box*             solidCalor;    //pointer to the solid Calor 
  G4LogicalVolume*   logicCalor;    //pointer to the logical Calor
  G4VPhysicalVolume* physiCalor;    //pointer to the physical Calor

  G4Box* MyBox;
  G4LogicalVolume* logicMyBox;
  G4VPhysicalVolume* physMyBox;

  G4Box* Box1;
  G4LogicalVolume* logicBox1;
  G4VPhysicalVolume* physBox1;

  G4Box* Box2;
  G4LogicalVolume* logicBox2;
  G4VPhysicalVolume* physBox2;

  G4Box* Box3;
  G4LogicalVolume* logicBox3;
  G4VPhysicalVolume* physBox3;

  G4Box* Box4;
  G4LogicalVolume* logicBox4;
  G4VPhysicalVolume* physBox4;



  G4Box* SciBox;
  G4LogicalVolume* logicSciBox;
  G4VPhysicalVolume* physSciBox;

  G4Box* IronBox1;
  G4LogicalVolume* logicIronBox1;
  G4VPhysicalVolume* physIronBox1;

  G4Box* AirBox;
  G4LogicalVolume* logicAirBox;
  G4VPhysicalVolume* physAirBox;

  G4Box* IronBox2;
  G4LogicalVolume* logicIronBox2;
  G4VPhysicalVolume* physIronBox2;

     
  G4Box*             solidLayer;    //pointer to the solid Layer 
  G4LogicalVolume*   logicLayer;    //pointer to the logical Layer
  G4VPhysicalVolume* physiLayer;    //pointer to the physical Layer
         
  G4Box*             solidAbsorber; //pointer to the solid Absorber
  G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
  G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
  G4Box*             solidGap;      //pointer to the solid Gap
  G4LogicalVolume*   logicGap;      //pointer to the logical Gap
  G4VPhysicalVolume* physiGap;      //pointer to the physical Gap
     
  G4UniformMagField* magField;      //pointer to the magnetic field
     
  ExRhadDetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
private:
    
  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();     
};

inline void ExRhadDetectorConstruction::SetWorldSize(G4double val)
{
  WorldSize = val;
}

inline void ExRhadDetectorConstruction::SetBox1Size(G4double val)
{
  Box1Size = val;
}

inline void ExRhadDetectorConstruction::SetBox2Size(G4double val)
{
  Box2Size = val;
}

inline void ExRhadDetectorConstruction::SetBox3Size(G4double val)
{
  Box3Size = val;
}

inline void ExRhadDetectorConstruction::SetBox4Size(G4double val)
{
  Box4Size = val;
}

inline void ExRhadDetectorConstruction::SetBox1Pos(G4double val)
{
  Box1Pos = val;
}

inline void ExRhadDetectorConstruction::SetBox2Pos(G4double val)
{
  Box2Pos = val;
}

inline void ExRhadDetectorConstruction::SetBox3Pos(G4double val)
{
  Box3Pos = val;
}

inline void ExRhadDetectorConstruction::SetBox4Pos(G4double val)
{
  Box4Pos = val;
}

inline void ExRhadDetectorConstruction::SetDetectorSpacing(G4double val)
{
  DetectorSpacing = val;
  WorldSize = Box1Size+Box2Size+Box3Size+Box4Size+3*DetectorSpacing/2;
  Box1Pos   = -WorldSize+Box1Size;
  Box2Pos = Box1Pos + Box1Size + Box2Size+DetectorSpacing;
  Box3Pos   = Box2Pos + Box2Size + Box3Size+DetectorSpacing;
  Box4Pos = Box3Pos + Box3Size + Box4Size+DetectorSpacing;
}

inline void ExRhadDetectorConstruction::SetBox1Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) Box1Material = pttoMaterial;
}

inline void ExRhadDetectorConstruction::SetBox2Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) Box2Material = pttoMaterial;
}

inline void ExRhadDetectorConstruction::SetBox3Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) Box3Material = pttoMaterial;
}

inline void ExRhadDetectorConstruction::SetBox4Material(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) Box4Material = pttoMaterial;
}



#endif

