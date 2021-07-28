/// G4 HEADER ///

// Compile only if GEANT4_BUILD_MULTITHREADED is set to ON when G4 was built
#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif

// Compile only if options, such as GEANT4_USE_OPENGL_X11 set to ON when G4 was built
#ifdef G4VIS_USE
  #include "G4VisExecutive.hh"
#endif

// Compile only if options, such as GEANT4_USE_QT set to ON when G4 was built
#ifdef G4UI_USE
  #include "G4UIExecutive.hh"
#endif

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

/// CUSTOM HEADER ///

#ifdef G4MULTITHREADED
  #include "InjectorUserActionInitialization.hh"
#else
  #include "InjectorPrimaryGeneratorAction.hh"
  #include "InjectorSteppingAction.hh"
  #include "InjectorTrackingAction.hh"
  #include "InjectorEventAction.hh"
  #include "InjectorRunAction.hh"
#endif

#include "InjectorDetectorConstruction.hh"
#include "InjectorPhysicsList.hh"

// Fot header files
#include "RunParameters.h"
#include "Crystal.h"
#include "mathematics.h"
#include <iostream>
#include <math.h>

// Global variables

#include "GlobalVariable.hh"
#include <iostream>
#include <vector>
using std::vector;

vector<G4String> vNtupName;
vector<G4String> vNtupTitle;
vector<G4String> vVarToFill;


/////// MAIN FUNCTION ///////

int main(int argc, char** argv) {

  if(argc>=5){
     std::string SEED = argv[4];
     int seed = std::atoi(SEED.c_str());
     CLHEP::HepRandom::setTheSeed(seed);
     std::cout<<"INFO:: HepRandom seed = "<<seed<<std::endl;
  }

  mathematics::set_seed(1);

  G4String tree_option = "";
  if(argc>=4){
    tree_option = argv[3];
  }

  // ------------------------------------------------------------------------------
  // define ntuples

  G4String NtupName[] = {
      "primary", "photon_emit", "xtal_leave", "amor_arrive", "amor_leave", "amd_arrive","amd_leave"};
  G4String NtupTitle[] = { 
      "particles at primary position",
      "particles at photon emission position",
      "particles leaving from crystal target",
      "particles arriving at amorphous target",
      "particles leaving from amorphous target",
      "particles arriving at amd",
      "particles leaving from amd"};
  const G4String VarToFill[] = { 
      "x", "y", "z", "px", "py", "pz", "e", "pdgId", "evtId", "t"};

  if(sizeof(NtupName)/sizeof(G4String)!=sizeof(NtupTitle)/sizeof(G4String)){
    G4cout<<"ERROR__:: inconsistence in size of ntuple names and titles!"<<G4endl;
  }
  
  for(unsigned int i_ntup=0;i_ntup<(sizeof(NtupName)/sizeof(G4String));i_ntup++){
    G4String ntup_name  = NtupName[i_ntup];
    G4String ntup_title  = NtupTitle[i_ntup];
    //Skip other trees for optimisation
    if(tree_option!="all" && tree_option!=ntup_name) continue;
    vNtupName.push_back(ntup_name);
    vNtupTitle.push_back(ntup_title);
  }
  
  for(unsigned int i_var=0;i_var<(sizeof(VarToFill)/sizeof(G4String));i_var++){
    G4String var_name  = VarToFill[i_var];
    vVarToFill.push_back(var_name);
  }

  //-------------------------------------------------------------------------------

  // show some flags 
  G4cout<<"INFO__:: Some G4 flags: "<<G4endl;
  #ifdef G4MULTITHREADED
    G4cout<<"INFO__:: G4MULTITHREADED is ON"<<G4endl;
  #else
    G4cout<<"INFO__:: G4MULTITHREADED is OFF"<<G4endl;
  #endif
  #ifdef G4MULTITHREADED
    G4cout<<"INFO__:: G4VIS_USE is ON"<<G4endl;
  #else
    G4cout<<"INFO__:: G4VIS_USE is OFF"<<G4endl;
  #endif
  #ifdef G4UI_USE
    G4cout<<"INFO__:: G4UI_USE is ON"<<G4endl;
  #else
    G4cout<<"INFO__:: G4UI_USE is OFF"<<G4endl;
  #endif

  // construct run manager
  #ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    if(argv[2]){
      if(atoi(argv[2])>0){
        runManager->SetNumberOfThreads(atoi(argv[2]));
      }
    }
    G4cout << "MT MODE ON " << runManager->GetNumberOfThreads() << G4endl;
  #else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "MT MODE OFF" << G4endl;
  #endif
    
  // activate UI-command base scorer
  G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(0);
    
  // choose the Random engine
  #ifndef G4MULTITHREADED
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  #endif
    
  // set mandatory initialization classes
  auto detector = new InjectorDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new InjectorPhysicsList());

  #ifndef G4MULTITHREADED
    runManager->SetUserAction(new InjectorPrimaryGeneratorAction(detector));
    runManager->SetUserAction(new InjectorEventAction());
    runManager->SetUserAction(new InjectorSteppingAction(detector));
    runManager->SetUserAction(new InjectorTrackingAction());
    runManager->SetUserAction(new InjectorRunAction());
  #else
    runManager->SetUserInitialization(new InjectorUserActionInitialization(detector));
  #endif
    
  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
  #endif
    
  if(argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  } else {
    // define visualization and UI terminal for interactive mode
    #ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      UI->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    #endif
        
    #ifdef G4VIS_USE
      delete visManager;
    #endif
  }
    
  // job termination
  delete runManager;
    
  return 0;
} 
/////// END OF MAIN FUNCTION ///////
