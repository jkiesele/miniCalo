#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
//#include "ExRhadVisManager.hh"
#include "G4VisExecutive.hh"
#endif

#include "ExRhadDetectorConstruction.hh"
#include "ExRhadPhysicsList.hh"


#include "ExRhadPrimaryGeneratorAction.hh"
#include "ExRhadRunAction.hh"
#include "ExRhadEventAction.hh"
#include "ExRhadSteppingAction.hh"
#include "ExRhadSteppingVerbose.hh"
#include "MyStackingAction.hh"
#include "MyTrackingAction.hh"

#include "FTFP_BERT.hh"

int main(int argc,char** argv) {

  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new ExRhadSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  ExRhadDetectorConstruction* detector = new ExRhadDetectorConstruction;

  runManager->SetUserInitialization(detector);

  auto physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new ExRhadPhysicsList());
//Use this if you want to use ExRhadPhysicsList
  runManager->SetUserInitialization(physicsList);

 
 G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#if defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
    }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new ExRhadPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new ExRhadRunAction);
  ExRhadEventAction* eventaction = new ExRhadEventAction;
  runManager->SetUserAction(eventaction);
  runManager->SetUserAction(new MyStackingAction());
  runManager->SetUserAction(new MyTrackingAction());
  runManager->SetUserAction(new ExRhadSteppingAction(detector, eventaction));
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute vis.mac");    
#if defined(G4UI_USE_XM) || defined(G4UI_USE_WIN32)
      // Customize the G4UIXm,Win32 menubar with a macro file :
      UI->ApplyCommand("/control/execute visTutor/gui.mac");
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
