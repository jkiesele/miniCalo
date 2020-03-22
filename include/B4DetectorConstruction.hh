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
// $Id: B4DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4DetectorConstruction.hh
/// \brief Definition of the B4DetectorConstruction class

#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "defines.h"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"

#include "sensorContainer.h"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4Material;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class B4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B4DetectorConstruction();
    virtual ~B4DetectorConstruction();

  public:
    enum geometry{
    	standard,
		homogenous,
		homogenous_ecal_only,
		ecal_only,
		ecal_only_hi_granular,
		hcal_only_irregular,
		ecal_only_irregular,
		homogenous_no_tracker,
		fullendcap
    };
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    void  DefineGeometry(geometry g);

    bool isActiveVolume(G4VPhysicalVolume*)const;

    const std::vector<sensorContainer>* getActiveSensors()const;

     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();


    G4VPhysicalVolume* createCellWheel(
            G4ThreeVector position,
            G4LogicalVolume* layerLV,
            G4double start_eta,
            G4double eta_width,
            G4double start_z,
            G4double z_length,
            G4int nphi,
            G4int layernum,
            G4int cellnum,
            G4Material* active_m,
            G4Material* abs_m=0,
            G4double abs_fraction=0);

    G4VPhysicalVolume* createLayer(
            G4LogicalVolume * caloLV,
            G4double start_eta,
            G4double end_eta,
            G4int n_cells_eta,
            G4int n_cells_phi,
            G4double start_z,
            G4double z_length,
            G4ThreeVector position,
            G4String name, int layernumber);
  
    void createCalo(G4LogicalVolume * worldLV,G4ThreeVector position,G4String name);
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger

    std::vector<sensorContainer> activecells_;

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    G4double layerThicknessEE,layerThicknessHB;
    G4double Calo_start_eta;
    G4double Calo_end_eta;
    G4double Calo_start_z;
    std::vector<G4double> layerThicknesses;
    G4int Ncalowedges;
    std::vector<G4int> layerGranularity;
    std::vector<G4int> layerSplitGranularity;
    G4double calorSizeXY;
    G4Material * m_vacuum, *m_pb, *m_pbtungsten, *m_silicon, *m_cu;
    G4int nofEELayers,nofHB, noTrackLayers;
    G4double calorThickness;


    G4double limit_in_calo_time_max_, limit_in_calo_energy_max_;
    G4double limit_world_time_max_,limit_world_energy_max_;


};

// inline functions

inline const std::vector<sensorContainer>* B4DetectorConstruction::getActiveSensors()const{
	return &activecells_;
}

     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

