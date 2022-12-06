#include "InjectorEventAction.hh"

#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"

#include "InjectorTrackingAction.hh"
#include "InjectorAnalysis.hh"

InjectorEventAction::InjectorEventAction()
{
    fVerboseLevel = 0;
}

InjectorEventAction::~InjectorEventAction(){
}

void InjectorEventAction::BeginOfEventAction(const G4Event* /* evt */ ){
}

void InjectorEventAction::EndOfEventAction(const G4Event* evt){
}

