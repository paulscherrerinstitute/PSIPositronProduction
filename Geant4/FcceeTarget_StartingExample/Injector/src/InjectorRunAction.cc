#include "InjectorRunAction.hh"

#include "InjectorRunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "InjectorAnalysis.hh"

#include "GlobalVariable.hh"

InjectorRunAction::InjectorRunAction(): G4UserRunAction(),fFileName("Injector"){

    // Create messenger, analysis manager
    fRunActionMessenger = new InjectorRunActionMessenger(this);
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Create directories
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(1);

    // Creating ntuples

    for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
      const G4String ntup_name  = vNtupName[i_ntup];
      const G4String ntup_title = vNtupTitle[i_ntup];
      analysisManager->CreateNtuple(ntup_name, ntup_title);
      for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
        const G4String varname = vVarToFill[ivar];
        if(varname=="pdgId"||varname=="evtId")
          analysisManager->CreateNtupleIColumn(varname);
        else
          analysisManager->CreateNtupleDColumn(varname);
      } // end loop of variables
      analysisManager->FinishNtuple();
    } // end loop of ntuples
}

InjectorRunAction::~InjectorRunAction(){
    delete G4AnalysisManager::Instance();
}

void InjectorRunAction::BeginOfRunAction(const G4Run* /*run*/){
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4String FileName = this->fFileName;
    analysisManager->OpenFile(FileName);
}

void InjectorRunAction::EndOfRunAction(const G4Run* /*run*/)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();
}

void InjectorRunAction::SetFileName(G4String value){
  fFileName = value;
}

