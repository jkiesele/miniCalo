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
// $Id: exampleB4a.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file exampleB4a.cc
/// \brief Main program of the B4a example

#include "B4DetectorConstruction.hh"
#include "B4aActionInitialization.hh"

#include "defines.h"

#ifdef G4MULTITHREADED
#undef G4MULTITHREADED
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include<string>
#include <unistd.h>

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4RandomTools.hh"

//limiter stuff
#include "G4SystemOfUnits.hh"

#include "G4PhysicsListHelper.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#define USEVIS

#include "B4PartGeneratorBase.hh"

#include "ExRhadPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4a [-m macro ] [-f filenumber=seed] [-o seedoffset (to be added to filenumber)] [-b beta]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}


void enable_physlimits(void)
{
    return;
  // cf. Geant 4 HyperNews, Forum "Physics List", Message 129
  // http://geant4-hn.slac.stanford.edu:5090/HyperNews/public/get/phys-list/129.html

  G4UserSpecialCuts *specialCuts = new G4UserSpecialCuts;
  G4StepLimiter     *stepLimiter = new G4StepLimiter;

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator *particleIterator = particleTable->GetIterator();
  // make sure you have called "G4RunManager::Initialize()" before

  particleIterator->reset();
  while ((*particleIterator)()) {
  // iterate through all known particles

    G4ParticleDefinition *particleDefinition = particleIterator->value();
    G4ProcessManager *processManager = particleDefinition->GetProcessManager();

    if (processManager && !particleDefinition->IsShortLived() && particleDefinition->GetPDGCharge() != 0) {
    // the process manager should exist, but we don't need to limit short-lived particles or neutrals

      processManager->AddDiscreteProcess(stepLimiter);
      processManager->AddDiscreteProcess(specialCuts);
      // these transportation-related processes are always applicable

    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //

  if ( argc > 11 ) {
    PrintUsage();

    return 1;
  }
  long rseed=0;
  G4String macro;
  G4String session="";
  G4String outfile="out";
  G4String betastr="";

  B4PartGeneratorBase::particle = "~g_rho0";
  B4PartGeneratorBase::beta = .1;

  B4DetectorConstruction::LArTankSize = 6*m;


#ifdef G4MULTITHREADED
  G4int nThreads = 0;
#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else if (G4String(argv[i]) == "-f" ) {

    	rseed += atoi(argv[i+1]);
    }
    else if (G4String(argv[i]) == "-o" ) {
        rseed += atoi(argv[i+1]);
        }
    else if (G4String(argv[i]) == "-l" ) {
        B4DetectorConstruction::LArTankSize = atof(argv[i+1])*m;
        betastr=argv[i+1];
        }
    else if (G4String(argv[i]) == "-a" ) {
        outfile += argv[i+1];
    }
    else {
      PrintUsage();
      return 1;
    }
  }  

  if(B4DetectorConstruction::LArTankSize < 2*m){
      G4cout << "LAr tank too small " << B4DetectorConstruction::LArTankSize  << G4endl;
      return -1;
  }

  betastr="_tank_"+betastr+"m_x"+betastr+"m_";

  outfile += betastr + std::to_string(rseed)+".root";
  rseed++;
  // Detect interactive mode (if no macro provided) and define UI session
  //
  //G4UIExecutive* ui = 0;
  //  if ( ! macro.size() ) {
  //      G4cout << "G4UIExecutive...." << G4endl;
  //    ui = new G4UIExecutive(argc, argv, session);
  //    G4cout << "G4UIExecutive done" << G4endl;
  //    G4cout << "G4UIExecutive done" << G4endl;
  //    G4cout << "G4UIExecutive done" << G4endl;
  //    usleep(2000000);
  //  }
  //  G4cout << "G4UIExecutive done" << G4endl;
  //

  // Choose the Random engine
  //

    B4PartGeneratorBase::seedsoffset_ = 800*rseed;

    G4cout << "random seed " << rseed << " particle " << B4PartGeneratorBase::particle<< G4endl;
  
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  auto runManager = new G4MTRunManager;
  if ( nThreads > 0 ) { 
    runManager->SetNumberOfThreads(nThreads);
  }  
#else
  auto runManager = new G4RunManager;
#endif

  G4cout << "B4DetectorConstruction" << G4endl;
  // Set mandatory initialization classes
  //
  auto detConstruction = new B4DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);


  auto physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());

  //noneed for detection bit
  //physicsList->RegisterPhysics(new ExRhadPhysicsList());

   runManager->SetUserInitialization(physicsList);

  //physicsList->RegisterPhysics(new ExRhadPhysicsList);

//  runManager->SetUserInitialization(new ExRhadPhysicsList);
    
  auto actionInitialization = new B4aActionInitialization(detConstruction);
  actionInitialization->setFilename(outfile);
  runManager->SetUserInitialization(actionInitialization);
  

  enable_physlimits();


  //SetThresholdWarningEnergy(500*MeV)

  // Initialize visualization
  //
#ifdef USEVIS
  G4cout << "G4VisExecutive" << G4endl;
  //ATree DAWNFILE HepRepFile HepRepXML OGL OGLI OGLIX OGLS OGLSX RayTracer VRML1FILE VRML2FILE gMocrenFile
   auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
 // G4VisManager* visManager = new G4VisExecutive("Quiet");
 // visManager->Initialize();
  //auto visManager = new G4VisExecutive();
     visManager->Initialize();
#endif
  // Get the pointer to the User Interface manager
     G4cout << "UImanager" << G4endl;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (argc!=1)   // batch mode
      {
        G4String command = "/control/execute ";
        G4String fileName = macro;
        UImanager->ApplyCommand(command+fileName);
      }
    else
      {  // interactive mode : define UI session

        G4cout << "G4UIExecutive ..." << G4endl;
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        G4cout << "G4UIExecutive done" << G4endl;

        G4cout << "init_vis" << G4endl;
        UImanager->ApplyCommand("/control/execute init_vis.mac");

        if (ui->IsGUI()){
            G4cout << "gui" << G4endl;
          UImanager->ApplyCommand("/control/execute gui.mac");
        }
        G4cout << "session" << G4endl;
        ui->SessionStart();
        delete ui;

      }


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

//  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
