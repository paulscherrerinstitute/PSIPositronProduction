#include "InjectorRunActionMessenger.hh"
#include "InjectorRunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

//--------------------------------------------------------------------//
//
InjectorRunActionMessenger::InjectorRunActionMessenger(
    InjectorRunAction* action):
  G4UImessenger(),fRunAction(action),
  fFileNameCmd(0)
{

  fFileNameCmd = new G4UIcmdWithAString("/output/filename",this);
  fFileNameCmd->SetParameterName("filename",true);
  fFileNameCmd->SetDefaultValue("Injector");

}

//--------------------------------------------------------------------//
//

InjectorRunActionMessenger::~InjectorRunActionMessenger()
{
  delete fFileNameCmd;
}

//--------------------------------------------------------------------//
//

void InjectorRunActionMessenger::SetNewValue(G4UIcommand * command,
    G4String newValue){
  if( command == fFileNameCmd ){
    fRunAction->SetFileName(newValue);
    return;
  }
}

//--------------------------------------------------------------------//
//
