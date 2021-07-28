#include "InjectorDetectorConstructionMessenger.hh"
#include "InjectorDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

InjectorDetectorConstructionMessenger::
InjectorDetectorConstructionMessenger(
        InjectorDetectorConstruction* mpga)
:fDet(mpga){
    fMyXtalDirectory = new G4UIdirectory("/xtal/");
    fMyXtalDirectory->SetGuidance("Xtal setup control commands.");

    //fXtalMaterialCmd = new G4UIcmdWithAString("/xtal/setMaterial",this);
    //fXtalMaterialCmd->SetGuidance("Set Xtal material.");
    //fXtalMaterialCmd->SetParameterName("xMat",true);
    //fXtalMaterialCmd->SetDefaultValue("G4_W");
    
    fXtalSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setSize",this);
    fXtalSizeCmd->SetGuidance("Set Xtal size.");
    fXtalSizeCmd->SetParameterName("xtalSizeX",
                                   "xtalSizeY",
                                   "xtalSizeZ",
                                   true);
    fXtalSizeCmd->SetDefaultValue(G4ThreeVector(6.,2.,6.));
    fXtalSizeCmd->SetDefaultUnit("mm");

    fAmorphousSizeCmd = new G4UIcmdWith3VectorAndUnit("/amorphous/setSize",this);
    fAmorphousSizeCmd->SetGuidance("Set Amorphous size.");
    fAmorphousSizeCmd->SetParameterName("amorphousSizeX",
                                   "amorphousSizeY",
                                   "amorphousSizeZ",
                                   true);
    fAmorphousSizeCmd->SetDefaultValue(G4ThreeVector(2.5,2.5,10.));
    fAmorphousSizeCmd->SetDefaultUnit("mm");


    fAmorphousDistanceCmd = new G4UIcmdWith3VectorAndUnit("/amorphous/setDistance",this);
    fAmorphousDistanceCmd->SetGuidance("Set Amorphous size.");
    fAmorphousDistanceCmd->SetParameterName("amorphousDistanceX",
                                   "amorphousDistanceY",
                                   "amorphousDistanceZ",
                                   true);
    fAmorphousDistanceCmd->SetDefaultValue(G4ThreeVector(0.,0.,2.));
    fAmorphousDistanceCmd->SetDefaultUnit("m");

    fDipoleFieldCmd = new G4UIcmdWith3VectorAndUnit("/dipole/setMagField",this);
    fDipoleFieldCmd->SetGuidance("Set dipole magnetic field values.");
    fDipoleFieldCmd->SetParameterName("magFieldValueX",
                                      "magFieldValueY",
                                      "magFieldValueZ",
                                      true);
    fDipoleFieldCmd->SetDefaultValue(G4ThreeVector(0.,1.,0.));
    fDipoleFieldCmd->SetDefaultUnit("tesla");

    fAMDOptionCmd = new G4UIcmdWithAString("/amd/setOption",this);
    fAMDOptionCmd->SetGuidance("Set AMD option.");
    fAMDOptionCmd->SetParameterName("amd_option",true);
    fAMDOptionCmd->SetDefaultValue("Analytic");
    fAMDFrontGapSizeCmd = new G4UIcmdWith3VectorAndUnit("/amd/setFrontGapSize",this);
    fAMDFrontGapSizeCmd->SetGuidance("Set AMD front gap size.");
    fAMDFrontGapSizeCmd->SetParameterName("amdFrontGapSizeX", "amdFrontGapSizeY", "amdFrontGapSizeZ", true);
    fAMDFrontGapSizeCmd->SetDefaultValue(G4ThreeVector(100.,100.,2.));
    fAMDFrontGapSizeCmd->SetDefaultUnit("mm");
    fAMDR1minCmd = new G4UIcmdWithADoubleAndUnit("/amd/setR1min",this);
    fAMDR1minCmd->SetGuidance("Set AMD R1min.");
    fAMDR1minCmd->SetParameterName("amdR1min",true);
    fAMDR1minCmd->SetDefaultValue(6);
    fAMDR1minCmd->SetDefaultUnit("mm");
    fAMDR1maxCmd = new G4UIcmdWithADoubleAndUnit("/amd/setR1max",this);
    fAMDR1maxCmd->SetGuidance("Set AMD R1max.");
    fAMDR1maxCmd->SetParameterName("amdR1max",true);
    fAMDR1maxCmd->SetDefaultValue(50);
    fAMDR1maxCmd->SetDefaultUnit("mm");
    fAMDR2minCmd = new G4UIcmdWithADoubleAndUnit("/amd/setR2min",this);
    fAMDR2minCmd->SetGuidance("Set AMD R2min.");
    fAMDR2minCmd->SetParameterName("amdR2min",true);
    fAMDR2minCmd->SetDefaultValue(20);
    fAMDR2minCmd->SetDefaultUnit("mm");
    fAMDR2maxCmd = new G4UIcmdWithADoubleAndUnit("/amd/setR2max",this);
    fAMDR2maxCmd->SetGuidance("Set AMD R2max.");
    fAMDR2maxCmd->SetParameterName("amdR2max",true);
    fAMDR2maxCmd->SetDefaultValue(50);
    fAMDR2maxCmd->SetDefaultUnit("mm");
    fAMDLengthCmd = new G4UIcmdWithADoubleAndUnit("/amd/setLength",this);
    fAMDLengthCmd->SetGuidance("Set AMD length.");
    fAMDLengthCmd->SetParameterName("amdLength",true);
    fAMDLengthCmd->SetDefaultValue(100);
    fAMDLengthCmd->SetDefaultUnit("mm");

    fAMDLinFringeActCmd = new G4UIcmdWithABool("/amd/LinFringe/isActive",this);
    fAMDLinFringeActCmd->SetGuidance("Activate AMD linear fringe field.");
    fAMDLinFringeActCmd->SetParameterName("amdLinFringeAct",true);
    fAMDLinFringeActCmd->SetDefaultValue(0);
    fAMDLinFringeB0Cmd = new G4UIcmdWithADoubleAndUnit("/amd/LinFringe/setB0",this);
    fAMDLinFringeB0Cmd->SetGuidance("Set B0 for AMD linear fringe field.");
    fAMDLinFringeB0Cmd->SetParameterName("amdLinFringeB0",true);
    fAMDLinFringeB0Cmd->SetDefaultValue(6);
    fAMDLinFringeB0Cmd->SetDefaultUnit("tesla");
    fAMDLinFringeKCmd = new G4UIcmdWithADouble("/amd/LinFringe/setK_T_MM",this);
    fAMDLinFringeKCmd->SetGuidance("Set K for AMD linear fringe field.");
    fAMDLinFringeKCmd->SetParameterName("amdLinFringeK",true);
    fAMDLinFringeKCmd->SetDefaultValue(0.5);
    //fAMDLinFringeKCmd->SetDefaultUnit("tesla/mm");
    fAMDLinFringeB0ZPosCmd = new G4UIcmdWithADoubleAndUnit("/amd/LinFringe/setB0ZPos",this);
    fAMDLinFringeB0ZPosCmd->SetGuidance("Set B0 z position for AMD linear fringe field. Z=0 at AMD entrance.");
    fAMDLinFringeB0ZPosCmd->SetParameterName("amdLinFringeB0ZPos",true);
    fAMDLinFringeB0ZPosCmd->SetDefaultValue(5);
    fAMDLinFringeB0ZPosCmd->SetDefaultUnit("mm");
    fAMDLinFringeTargZPosCmd = new G4UIcmdWithADoubleAndUnit("/amd/LinFringe/setTargZPos",this);
    fAMDLinFringeTargZPosCmd->SetGuidance("Set target z position for AMD linear fringe field. Z=0 at AMD entrance.");
    fAMDLinFringeTargZPosCmd->SetParameterName("amdLinFringeTargZPos",true);
    fAMDLinFringeTargZPosCmd->SetDefaultValue(-2);
    fAMDLinFringeTargZPosCmd->SetDefaultUnit("mm");
}

InjectorDetectorConstructionMessenger::
~InjectorDetectorConstructionMessenger(){
    //delete fXtalMaterialCmd;
    delete fXtalSizeCmd;
    delete fAmorphousSizeCmd;
    delete fAmorphousDistanceCmd;
    delete fDipoleFieldCmd;
    delete fAMDOptionCmd;
    delete fAMDFrontGapSizeCmd;
    delete fAMDR1minCmd;
    delete fAMDR1maxCmd;
    delete fAMDR2minCmd;
    delete fAMDR2maxCmd;
    delete fAMDLengthCmd;
    delete fAMDLinFringeActCmd;
    delete fAMDLinFringeB0Cmd;
    delete fAMDLinFringeKCmd;
    delete fAMDLinFringeB0ZPosCmd;
    delete fAMDLinFringeTargZPosCmd;
}

void InjectorDetectorConstructionMessenger::SetNewValue(
                                                    G4UIcommand *command,
                                                    G4String newValue){
    //if(command==fXtalMaterialCmd ){
    //    fDet->SetXtalMaterial(newValue);
    //}
    if(command==fXtalSizeCmd ){
        fDet->SetXtalSize(fXtalSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fAmorphousSizeCmd ){
        fDet->SetAmorphousSize(fAmorphousSizeCmd->GetNew3VectorValue(newValue));
    }

    if(command==fAmorphousDistanceCmd ){
        fDet->SetAmorphousDistance(fAmorphousDistanceCmd->GetNew3VectorValue(newValue));
    }

    if(command==fDipoleFieldCmd ){
        fDet->SetDipoleFieldValue(fDipoleFieldCmd->GetNew3VectorValue(newValue));
    }

    if(command==fAMDOptionCmd ){
        fDet->SetAMDOption(newValue);
    }
    if(command==fAMDFrontGapSizeCmd ){
        fDet->SetAMDFrontGapSize(fAMDFrontGapSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fAMDR1minCmd ){
        fDet->SetAMDR1min(fAMDR1minCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDR1maxCmd ){
        fDet->SetAMDR1max(fAMDR1maxCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDR2minCmd ){
        fDet->SetAMDR2min(fAMDR2minCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDR2maxCmd ){
        fDet->SetAMDR2max(fAMDR2maxCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDLengthCmd ){
        fDet->SetAMDLength(fAMDLengthCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDLinFringeActCmd ){
        fDet->SetAMDLinFringeAct(fAMDLinFringeActCmd->GetNewBoolValue(newValue));
    }
    if(command==fAMDLinFringeB0Cmd ){
        fDet->SetAMDLinFringeB0(fAMDLinFringeB0Cmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDLinFringeKCmd ){
        fDet->SetAMDLinFringeK(fAMDLinFringeKCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDLinFringeB0ZPosCmd ){
        fDet->SetAMDLinFringeB0ZPos(fAMDLinFringeB0ZPosCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAMDLinFringeTargZPosCmd ){
        fDet->SetAMDLinFringeTargZPos(fAMDLinFringeTargZPosCmd->GetNewDoubleValue(newValue));
    }

}

/*
G4String InjectorDetectorConstructionMessenger::GetCurrentValue(
                                                    G4UIcommand * command){
    G4String cv;
    
    //if( command==fXtalMaterialCmd ){
    //    cv = fDet->GetXtalMaterial();
    //}
    if( command==fXtalSizeCmd ){
        cv = fXtalSizeCmd->ConvertToString(fDet->GetXtalSize(),"mm");
    }
    if( command==fAmorphousSizeCmd ){
        cv = fAmorphousSizeCmd->ConvertToString(fDet->GetAmorphousSize(),"mm");
    }
    if( command==fAmorphousDistanceCmd ){
        cv = fAmorphousDistanceCmd->ConvertToString(fDet->GetAmorphousDistance(),"m");
    }
    if( command==fDipoleFieldCmd ){
        cv = fDipoleFieldCmd->ConvertToString(fDet->GetDipoleFieldValue(),"tesla");
    }
    if( command==fAMDFrontGapSizeCmd ){
        cv = fAMDFrontGapSizeCmd->ConvertToString(fDet->GetAMDFrontGapSize(),"mm");
    }
    if( command==fAMDOptionCmd ){
        cv = fDet->GetAMDOption();
    }
    if( command==fAMDR1minCmd ){
        cv = fAMDR1minCmd->ConvertToString(fDet->GetAMDR1min(),"mm");
    }
    if( command==fAMDR1maxCmd ){
        cv = fAMDR1maxCmd->ConvertToString(fDet->GetAMDR1max(),"mm");
    }
    if( command==fAMDR2minCmd ){
        cv = fAMDR2minCmd->ConvertToString(fDet->GetAMDR2min(),"mm");
    }
    if( command==fAMDR2maxCmd ){
        cv = fAMDR2maxCmd->ConvertToString(fDet->GetAMDR2max(),"mm");
    }
    if( command==fAMDLengthCmd ){
        cv = fAMDLengthCmd->ConvertToString(fDet->GetAMDLength(),"mm");
    }
    if( command==fAMDLinFringeActCmd ){
        cv = fAMDLinFringeActCmd->ConvertToString(fDet->GetAMDLinFringeAct());
    }
    if( command==fAMDLinFringeB0Cmd ){
        cv = fAMDLinFringeB0Cmd->ConvertToString(fDet->GetAMDLinFringeB0(),"tesla");
    }
    if( command==fAMDLinFringeKCmd ){
        //cv = fAMDLinFringeKCmd->ConvertToString(fDet->GetAMDLinFringeK(),"tesla/mm");
        cv = fAMDLinFringeKCmd->ConvertToString(fDet->GetAMDLinFringeK());
    }
    if( command==fAMDLinFringeB0ZPosCmd ){
        cv = fAMDLinFringeB0ZPosCmd->ConvertToString(fDet->GetAMDLinFringeB0ZPos(),"mm");
    }
    if( command==fAMDLinFringeTargZPosCmd ){
        cv = fAMDLinFringeTargZPosCmd->ConvertToString(fDet->GetAMDLinFringeTargZPos(),"mm");
    }
    
    return cv;
}
*/

