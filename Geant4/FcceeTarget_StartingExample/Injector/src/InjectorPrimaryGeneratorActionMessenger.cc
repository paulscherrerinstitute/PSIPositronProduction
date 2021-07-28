#include "InjectorPrimaryGeneratorActionMessenger.hh"
#include "InjectorPrimaryGeneratorAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4SystemOfUnits.hh"

InjectorPrimaryGeneratorActionMessenger::InjectorPrimaryGeneratorActionMessenger(
    InjectorPrimaryGeneratorAction* action):
  G4UImessenger(),fPrimaryGeneratorAction(action),
  //fXtalThickCmd(0),
  fPrimaryEnergyCmd(0),
  fPrimarySigmaXYCmd(0),
  fPrimarySigmaZCmd(0),
  fPrimarySigmaTCmd(0),
  fPrimarySigmaPxyCmd(0)
{

  //fXtalThickCmd = new G4UIcmdWithADoubleAndUnit("/gun/XtalThick",this);
  //fXtalThickCmd->SetParameterName("XtalThick",true);
  //fXtalThickCmd->SetDefaultValue(1.4);
  //fXtalThickCmd->SetDefaultUnit("mm");

  fPrimaryEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/primaryEnergy",this);
  fPrimaryEnergyCmd->SetParameterName("primaryEnergy",true);
  fPrimaryEnergyCmd->SetDefaultValue(3000);
  fPrimaryEnergyCmd->SetDefaultUnit("MeV");

  fPrimarySigmaECmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaUserE",this);
  fPrimarySigmaECmd->SetParameterName("sigmaUserE",true);
  fPrimarySigmaECmd->SetDefaultValue(3.0);
  fPrimarySigmaECmd->SetDefaultUnit("MeV");

  fPrimarySigmaXYCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaUserXY",this);
  fPrimarySigmaXYCmd->SetParameterName("sigmaUserXY",true);
  fPrimarySigmaXYCmd->SetDefaultValue(2.5);
  fPrimarySigmaXYCmd->SetDefaultUnit("mm");

  fPrimarySigmaZCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaUserZ",this);
  fPrimarySigmaZCmd->SetParameterName("sigmaUserZ",true);
  fPrimarySigmaZCmd->SetDefaultValue(0.3);
  fPrimarySigmaZCmd->SetDefaultUnit("mm");

  //fPrimarySigmaTCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaUserT",this);
  fPrimarySigmaTCmd = new G4UIcmdWithADouble("/gun/sigmaUserT",this);
  fPrimarySigmaTCmd->SetParameterName("sigmaUserT",true);
  fPrimarySigmaTCmd->SetDefaultValue(0.3);
  // Define customised units
  //new G4UnitDefinition("millimeter/c", "mm/c", "Time", CLHEP::millimeter/CLHEP::c_light);
  //fPrimarySigmaTCmd->SetDefaultUnit("mm/c");

  fPrimarySigmaPxyCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaUserPxy",this);
  fPrimarySigmaPxyCmd->SetParameterName("sigmaUserPxy",true);
  fPrimarySigmaPxyCmd->SetDefaultValue(1.0e-3);
  fPrimarySigmaPxyCmd->SetDefaultUnit("MeV");
}

InjectorPrimaryGeneratorActionMessenger::~InjectorPrimaryGeneratorActionMessenger()
{
  //delete fXtalThickCmd;
  delete fPrimaryEnergyCmd;
  delete fPrimarySigmaECmd;
  delete fPrimarySigmaXYCmd;
  delete fPrimarySigmaZCmd;
  delete fPrimarySigmaTCmd;
  delete fPrimarySigmaPxyCmd;
}

void InjectorPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand * command,
    G4String newValue){
  //if( command == fXtalThickCmd ){
  //  fPrimaryGeneratorAction->SetXtalThick(fXtalThickCmd->GetNewDoubleValue(newValue));
  //  return;
  //}
  if( command == fPrimaryEnergyCmd ){
    fPrimaryGeneratorAction->SetPrimaryEnergy(fPrimaryEnergyCmd->GetNewDoubleValue(newValue));
    return;
  }
  if( command == fPrimarySigmaECmd ){
    fPrimaryGeneratorAction->SetSigmaE(fPrimarySigmaECmd->GetNewDoubleValue(newValue));
    return;
  }
  if( command == fPrimarySigmaXYCmd ){
    fPrimaryGeneratorAction->SetSigmaXY(fPrimarySigmaXYCmd->GetNewDoubleValue(newValue));
    return;
  }
  if( command == fPrimarySigmaZCmd ){
    fPrimaryGeneratorAction->SetSigmaZ(fPrimarySigmaZCmd->GetNewDoubleValue(newValue));
    return;
  }
  if( command == fPrimarySigmaTCmd ){
    fPrimaryGeneratorAction->SetSigmaT(fPrimarySigmaTCmd->GetNewDoubleValue(newValue));
    return;
  }
  if( command == fPrimarySigmaPxyCmd ){
    fPrimaryGeneratorAction->SetSigmaPxy(fPrimarySigmaPxyCmd->GetNewDoubleValue(newValue));
    return;
  }
}

