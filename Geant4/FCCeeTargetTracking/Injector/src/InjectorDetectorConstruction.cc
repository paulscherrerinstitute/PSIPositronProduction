#include <TFile.h>
#include <TH2F.h>
//#include <TGraph2D.h>

#include "InjectorDetectorConstruction.hh"
#include "InjectorField.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4TransportationManager.hh"

#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4UserLimits.hh"

#include "G4ThreeVector.hh"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class TFile;
class TH2F;
//class TGraph2D;

InjectorDetectorConstruction::InjectorDetectorConstruction():
fWorldLogic(0){
    
    fXtal = true;
    fXtalSize = G4ThreeVector(100 * CLHEP::millimeter,
                              100 * CLHEP::millimeter,
                              1 * CLHEP::millimeter);
    fXtalLocation = G4ThreeVector(0, 0, 0.5*CLHEP::millimeter);

    //SetXtalMaterial("G4_W");

    fAmorphous = true;
    fAmorphousSize = G4ThreeVector(100 * CLHEP::millimeter,
                                   100 * CLHEP::millimeter,
                                    14 * CLHEP::millimeter);
    fAmorphousLocation = G4ThreeVector(0, 0, 8*CLHEP::millimeter);

    fAMDFrontGap = 1;
    fAMDFrontGapSize = G4ThreeVector( 100 * CLHEP::millimeter, 
                                      100 * CLHEP::millimeter, 
				      2 * CLHEP::millimeter);

    fAMD = true;
    fAMDOption = "";
    fAMDR1min = 5.8*mm; //[mm] 
    fAMDR1max = 50*mm; //[mm]
    fAMDR2min = 20*mm; //[mm]
    fAMDR2max = 50*mm; //[mm]
    fAMDLength = 100*mm; // [mm]
    fAMDLinFringeAct = 0;
    fAMDLinFringeB0 = 6*tesla; // [T]
    fAMDLinFringeK = 0.5; // [T/mm]
    fAMDLinFringeB0ZPos = 5*mm; // [mm]
    fAMDLinFringeTargZPos = -2*mm; // [mm]
    fAMDLocation = G4ThreeVector(0, 0, 67.5*CLHEP::millimeter);

    fMessenger = new InjectorDetectorConstructionMessenger(this);
}

InjectorDetectorConstruction::~InjectorDetectorConstruction(){
  //histFieldMapR->Clear();
  //histFieldMapZ->Clear();
  //histFieldMapR->Delete();
  //histFieldMapZ->Delete();
}

void InjectorDetectorConstruction::DefineMaterials(){
}

G4VPhysicalVolume* InjectorDetectorConstruction::Construct(){
    fWorldSize = G4ThreeVector(1  * CLHEP::meter,
                               1  * CLHEP::meter,
                               10 * CLHEP::meter);
    fWorldMaterial = G4NistManager::
        Instance()->FindOrBuildMaterial("G4_Galactic");
    
    fWorldSolid = new G4Box("World",
                            fWorldSize.x()/2.,
                            fWorldSize.y()/2.,
                            fWorldSize.z()/2.);
    
    fWorldLogic = new G4LogicalVolume(fWorldSolid,
                                      fWorldMaterial,
                                      "World");
    
    fWorldPhysical = new G4PVPlacement(0,
                                       G4ThreeVector(),
                                       fWorldLogic,
                                       "World",
                                       0,
                                       false,
                                       0);
    
    //FillMagFieldHist();
    //FillMagFieldGraph();

    if(fXtalSize.z()==0){ 
      fXtal = false;
      fXtalSize = G4ThreeVector( 0,0,0 );
      fXtalLocation = G4ThreeVector(0, 0, 0);
    }  
    if(fXtal){
      ConstructXtalTarget();
    } else {
      std::cout<<"INFO:: Xtal not constructed!"<<std::endl;
    }

    //if(fAmorphousDistance.z()!=0 && fDipoleFieldValue.mag()!=0){
    //  ConstructDipoleField();
    //}else{
    //  std::cout<<"INFO:: Dipole not constructed!"<<std::endl;
    //}
    if(fAmorphousDistance.z()!=0 && fDipoleFieldValue.mag()!=0){
      std::cout<<"INFO:: Dipole field activated."<<std::endl;
    }else{
      std::cout<<"INFO:: Dipole field not activated!"<<std::endl;
    }

    if(fAmorphousSize.z()==0){ 
      fAmorphous = false;
      fAmorphousSize = G4ThreeVector( 0,0,0 );
      fAmorphousLocation = G4ThreeVector(0, 0, 0);
    }
    if(fAmorphous){
      ConstructAmorphousTarget();
    }else{
      std::cout<<"INFO:: Amorphous not constructed!"<<std::endl;
    }

    if(fAMDFrontGapSize.z()==0){
      fAMDFrontGap = 0;
    }
    if(fAMDFrontGap){
      ConstructAMDFrontGap();
    }else{
      std::cout<<"INFO:: AMD front gap not constructed!"<<std::endl;
    }

    if(fAMDLength==0){ 
      fAMD = false;
    }
    if(fAMD){
      ConstructAMD();
    }else{
      std::cout<<"INFO:: AMD not constructed!"<<std::endl;
    }

    if(fAMDLinFringeAct){
      std::cout<<"INFO:: AMD linear fringe field activated."<<std::endl;
    }else{
      std::cout<<"INFO:: AMD fringe field not activated!"<<std::endl;
    }

    G4UserLimits* fUserLimits = new G4UserLimits();
    //fUserLimits->SetMaxAllowedStep(1*mm);
    //fUserLimits->SetUserMaxTrackLength(500.*m);
    //fUserLimits->SetUserMaxTime(10*ms);
    fUserLimits->SetUserMaxTime(1*us);
    //fUserLimits->SetUserMinEkine(0.1*MeV);
    fWorldLogic->SetUserLimits(fUserLimits);

    return fWorldPhysical;
}  

void InjectorDetectorConstruction::ConstructXtalTarget(){
  std::cout<<"INFO:: Constructing Xtal. Target ..."<<std::endl;

  fXtalSolid = new G4Box("Target", 
                         fXtalSize.x()/2, 
			 fXtalSize.y()/2, 
			 fXtalSize.z()/2);

//  fXtalLogic = new G4LogicalVolume(fXtalSolid, 
//                                   fXtalMaterial, 
//				   "Target");
  fXtalLogic = new G4LogicalVolume(fXtalSolid, 
                                   G4NistManager::Instance()->FindOrBuildMaterial("G4_W"),
				   "Target");

  G4double world_opz = 0 * CLHEP::meter; // [mm]
  fXtalLocation = G4ThreeVector(0, 
                                0, 
				fXtalSize.z()/2 + world_opz);

  fXtalPhysical = new G4PVPlacement( 0,
                                     fXtalLocation,
                                     fXtalLogic,
				     "Target",
                                     fWorldLogic,
                                     false,
                                     0);
}

void InjectorDetectorConstruction::ConstructAmorphousTarget(){
  std::cout<<"INFO:: Constructing Amor. Target ..."<<std::endl;
  fAmorphousSolid = new G4Box("Amorphous", fAmorphousSize.x()/2, fAmorphousSize.y()/2, fAmorphousSize.z()/2);

  fAmorphousLogic = new G4LogicalVolume( fAmorphousSolid,
                                         G4NistManager::Instance()->FindOrBuildMaterial("G4_W"),
                                         "Amorphous");

  G4double world_opz = 0 * CLHEP::meter; // [mm]
  fAmorphousLocation = G4ThreeVector( 0,
                                      0,
                                      world_opz + fXtalSize.z() + fAmorphousDistance.z() + fAmorphousSize.z()/2);

  fAmorphousPhysical = new G4PVPlacement( 0,
                                          fAmorphousLocation,
                                          fAmorphousLogic,
                                          "Amorphous",
                                          fWorldLogic,
                                          false,
                                          0);
}

void InjectorDetectorConstruction::ConstructAMDFrontGap(){
  std::cout<<"INFO:: Constructing AMD front gap ..."<<std::endl;

  fAMDFrontGapSolid = new G4Box("AMDFrontGap", fAMDFrontGapSize.x()/2, fAMDFrontGapSize.y()/2, fAMDFrontGapSize.z()/2);

  fAMDFrontGapLogic = new G4LogicalVolume( fAMDFrontGapSolid,
                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
                                      "AMDFrontGap");

  fAMDFrontGapLocation = G4ThreeVector( 0,
                                        0,
                         fXtalLocation.z() + fXtalSize.z()/2 + fAmorphousDistance.z() + fAmorphousSize.z() + fAMDFrontGapSize.z()/2);

  fAMDFrontGapPhysical = new G4PVPlacement( 0,
                                            fAMDFrontGapLocation,
                                            fAMDFrontGapLogic,
                                            "AMDFrontGap",
                                            fWorldLogic,
                                            false,
                                            0);
}

//void InjectorDetectorConstruction::ConstructAMD(){
//  std::cout<<"INFO:: Constructing AMD ..."<<std::endl;
//
//  //cout<<":::::::::::: This is ConstructAMD ::::::::::::::::"<<endl;
//
//  fAMDSolid = new G4Cons("AMD", fAMDR1min, fAMDR1max, fAMDR2min, fAMDR2max, fAMDLength/2, 0*deg, 360*deg);
//
//  fAMDLogic = new G4LogicalVolume( fAMDSolid,
//                                   G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu"),
//                                   "AMD");
//
//  fAMDLocation = G4ThreeVector( 0,
//                                0,
//                                fXtalLocation.z() + fXtalSize.z()/2 + fAmorphousDistance.z() + fAmorphousSize.z() + fAMDFrontGapSize.z() + fAMDLength/2);
//
//  fAMDPhysical = new G4PVPlacement( 0,
//                                    fAMDLocation,
//                                    fAMDLogic,
//                                    "AMD",
//                                    fWorldLogic,
//                                    false,
//                                    0);
//
//  fAMDVacSolid = new G4Cons("AMDVac", 0, fAMDR1min, 0, fAMDR2min, fAMDLength/2, 0*deg, 360*deg);
//
//  fAMDVacLogic = new G4LogicalVolume( fAMDVacSolid,
//                                      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
//                                      "AMDVac");
//
//  fAMDVacPhysical = new G4PVPlacement( 0,
//                                       fAMDLocation,
//                                       fAMDVacLogic,
//                                       "AMDVac",
//                                       fWorldLogic,
//                                       false,
//                                       0);
//}

void InjectorDetectorConstruction::ConstructAMD(){

  std::cout<<"INFO:: Constructing AMD ..."<<std::endl;

  // z position of AMD entrance
  double z0 = fXtalLocation.z() + fXtalSize.z()/2 + fAmorphousDistance.z() + fAmorphousSize.z() + fAMDFrontGapSize.z();

  double W = 8.33; // [mm], turn width
  double G; // [mm], turn gap
  double Rmax; // [mm], outer radius
  double z1 = z0 + W/2.0; // [mm], 1st turn center
  int N; // number of coil turns

  double R[50];
  double R_Lin_Type1_R3p5_I13p0_F25_G0p2_N12[] = {3.50,5.38,7.25,9.13,11.00,12.88,14.75,16.63,18.50,20.38,22.25,24.13};
  double R_Lin_Type1_R3p5_I13p8_F25[] = {3.50,7.25,11.00,14.75,18.50,22.25,26.00,29.75,33.50,37.25,41.00,44.75,48.60,52.45};
  double R_Lin_Type1_R6p5_I13p8_F25[] = {6.50,10.25,14.00,17.75,21.50,25.25,29.00,32.75,36.50,40.25,44.00,47.75,51.60,55.45};
  double R_Exp_Type1_R3p5_I13p8_F25[] = {3.49,12.14,19.33,24.99,29.42,32.96,35.88,38.36,40.50,42.40,44.10,45.66,47.18,48.61}; 
  double R_Exp_Type1_R6p5_I13p8_F25[] = {6.49,15.14,22.33,27.99,32.42,35.96,38.88,41.36,43.50,45.40,47.10,48.66,50.18,51.61};
  double R_Inv_Type1_R2p5_I13p8_F25[] = {2.45,4.27,6.19,8.24,10.43,12.78,15.34,18.16,21.27,24.79,28.84,33.64,39.66,47.58};
  double R_Inv_Type1_R6p5_I13p8_F25[] = {6.45,8.27,10.19,12.24,14.43,16.78,19.34,22.16,25.27,28.79,32.84,37.64,43.66,51.58};
  double R_Inv_Type2_R2p2_I13p8_F25[] = {2.17,3.43,4.76,6.17,7.69,9.33,11.12,13.12,15.35,17.94,21.02,24.87,30.22,38.87};
  double R_Inv_Type2_R6p5_I13p8_F25[] = {6.52,7.78,9.11,10.52,12.04,13.68,15.47,17.47,19.70,22.29,25.37,29.22,34.57,43.22};
  double R_Inv_Type2_R6p5_I25_F25[] = {6.52,7.78,9.11,10.52,12.04,13.68,15.47,17.47,19.70,22.29,25.37,29.22,34.57,43.22};
  double R_Inv_Type3_R2p0_I13p8_F25[] = {1.95,2.84,3.96,5.54,8.20,12.43,17.26,22.16,27.06,31.96,36.87,41.77,46.77,51.77};
  double R_Inv_Type3_R6p5_I13p8_F25[] = {6.45,7.34,8.46,10.04,12.70,16.93,21.76,26.66,31.56,36.46,41.37,46.27,51.27,56.27};

  if(fAMDOption=="Lin_Type1_R3.5_I13.0_F25_G0.2_N12"){
    G = 0.2;
    Rmax = 40; // [mm]
    N = sizeof(R_Lin_Type1_R3p5_I13p0_F25_G0p2_N12)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Lin_Type1_R3p5_I13p0_F25_G0p2_N12[i]; }
  }
  if(fAMDOption=="Lin_Type1_R3.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Lin_Type1_R3p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Lin_Type1_R3p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Lin_Type1_R6.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Lin_Type1_R6p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Lin_Type1_R6p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Exp_Type1_R3.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Exp_Type1_R3p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Exp_Type1_R3p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Exp_Type1_R6.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Exp_Type1_R6p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Exp_Type1_R6p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type1_R2.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type1_R2p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type1_R2p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type1_R6.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type1_R6p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type1_R6p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type2_R2.2_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type2_R2p2_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type2_R2p2_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type2_R6.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type2_R6p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type2_R6p5_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type2_R6.5_I25_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type2_R6p5_I25_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type2_R6p5_I25_F25[i]; }
  }
  if(fAMDOption=="Inv_Type3_R2.0_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type3_R2p0_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type3_R2p0_I13p8_F25[i]; }
  }
  if(fAMDOption=="Inv_Type3_R6.5_I13.8_F25"){
    G = 0.8;
    Rmax = 60; // [mm]
    N = sizeof(R_Inv_Type3_R6p5_I13p8_F25)/sizeof(double);
    for (int i=0;i<N;i++) { R[i] = R_Inv_Type3_R6p5_I13p8_F25[i]; }
  }

  double Z[50];
  for (int i=0;i<N;i++) { Z[i] = z1 + (W+G) * i; }

  // define AMD solid
  G4Tubs* pAMDSolid[50];
  for (int i=0;i<N;i++) { 
    char buffer[50];
    sprintf(buffer,"AMDSolid_%i",i);
    G4String name = buffer;
    pAMDSolid[i] = new G4Tubs(name, R[i], Rmax, W/2.0, 0*deg, 360*deg);
  }

  // define AMD logical
  G4LogicalVolume* pAMDLogic[50];
  for (int i=0;i<N;i++) { 
    char buffer[50];
    sprintf(buffer,"AMDLogical_%i",i);
    G4String name = buffer;
    pAMDLogic[i] = new G4LogicalVolume( pAMDSolid[i],
                                        G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu"),
                                        name);
  }

  // define AMD physical
  G4VPhysicalVolume* pAMDPhysical[50];
  for (int i=0;i<N;i++) { 
    char buffer[50];
    sprintf(buffer,"AMDPhysical_%i",i);
    G4String name = buffer;
    pAMDPhysical[i] = new G4PVPlacement( 0,
                                         G4ThreeVector( 0, 0, Z[i]),
                                         pAMDLogic[i],
                                         name,
                                         fWorldLogic,
                                         false,
                                         0);
  }

  // define AMD Vacuum solid
  G4Tubs* pAMDVacSolidFirst = new G4Tubs("AMDVacSolidFirst", 0, R[0],   W/2.0, 0*deg, 360*deg);
  G4Tubs* pAMDVacSolidLast  = new G4Tubs("AMDVacSolidLast",  0, R[N-1], W/2.0, 0*deg, 360*deg);

  // define AMD Vacuum logical
  G4LogicalVolume* pAMDVacLogicFirst = new G4LogicalVolume( pAMDVacSolidFirst,
                                                            G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
							    "AMDVacLogicFirst");
  G4LogicalVolume* pAMDVacLogicLast  = new G4LogicalVolume( pAMDVacSolidLast,
                                                            G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
							    "AMDVacLogicLast");

  // define AMD Vacuum physical
  G4VPhysicalVolume* pAMDVacPhysicalFirst = new G4PVPlacement( 0,
                                            G4ThreeVector( 0, 0, Z[0]),
                                            pAMDVacLogicFirst,
                                            "AMDVacPhysicalFirst",
                                            fWorldLogic,
                                            false,
                                            0);
  G4VPhysicalVolume* pAMDVacPhysicalLast  = new G4PVPlacement( 0,
                                            G4ThreeVector( 0, 0, Z[N-1]),
                                            pAMDVacLogicLast,
                                            "AMDVacPhysicalLast",
                                            fWorldLogic,
                                            false,
                                            0);

  // reset some AMD parameters
  fAMDR1min = R[0];
  fAMDR1max = Rmax;
  fAMDR2min = R[N-1];
  fAMDR2max = Rmax;
  fAMDLength = N*W + (N-1)*G; // [mm]
  std::cout<<"DEBUG:: AMD length reset here to "<<fAMDLength<<" [mm] in ConstructAMD() function.\n"<<std::endl;
}

//void InjectorDetectorConstruction::ConstructDipoleField(){
//
//  // define magnetic volumes
//
//  fMagSolid = new G4Box( "Magnet",
//                         fAmorphousSize.x()/2.0,
//                         fAmorphousSize.y()/2.0,
//                         fAmorphousDistance.z()*0.98/2 );
//  fMagLogic = new G4LogicalVolume( fMagSolid,
//                                   G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"),
//                                   "Magnet");
//
//  // create a field
//  fDipoleField = new G4UniformMagField( fDipoleFieldValue ); 
//
//  // set it as default field
//  fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager(); 
//  fFieldMgr->SetDetectorField(fDipoleField);
//
//  // create objects that calculate the trajectory
//  fFieldMgr->CreateChordFinder(fDipoleField);
//
//  // attach field manager to a logical volume
//  fMagLogic->SetFieldManager(fFieldMgr, true);
//
//  // create physical volume
//  fMagPhysical = new G4PVPlacement( 0,
//                                    fXtalLocation + G4ThreeVector(0,0,fXtalSize.z()/2 + fAmorphousDistance.z()/2),
//                                    fMagLogic,
//				    "Magnet",
//                                    fWorldLogic,
//                                    false,
//                                    0);
//}

void InjectorDetectorConstruction::ConstructSDandField(){
  //cout<<":::::::::::: This is ConstructSDandField ::::::::::::::::"<<endl;
  //std::cout<<"INFO:: Constructing field ..."<<std::endl;
  InjectorField* myField = new InjectorField(this);
  G4FieldManager* fieldMgr
                  = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(myField);
  fieldMgr->CreateChordFinder(myField);
}

void InjectorDetectorConstruction::SetXtalMaterial(const G4String& name){

  G4Material* vMaterial;
  vMaterial = G4Material::GetMaterial(name, true); // warning turned on
  if(!vMaterial){
    vMaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
  }

  if (vMaterial && vMaterial != fXtalMaterial) {
    fXtalMaterial = vMaterial;
    if(fXtalLogic){
      fXtalLogic->SetMaterial(vMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

G4String InjectorDetectorConstruction::GetXtalMaterial(){
  if(fXtalMaterial) {
    return fXtalMaterial->GetName();
  }
  return "";
}

void InjectorDetectorConstruction::SetXtalSize(G4ThreeVector size) {
  if(fXtalSize != size) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fXtalSize = size;
  }
}

void InjectorDetectorConstruction::SetAmorphousSize(G4ThreeVector size) {
  if(fAmorphousSize != size) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAmorphousSize = size;
  }
}

void InjectorDetectorConstruction::SetAmorphousDistance(G4ThreeVector distance) {
  if(fAmorphousDistance != distance) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAmorphousDistance = distance;
  }
}

void InjectorDetectorConstruction::SetAMDFrontGapSize(G4ThreeVector size) {
  if(fAMDFrontGapSize != size) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDFrontGapSize = size;
  }
}

void InjectorDetectorConstruction::SetAMDOption(const G4String& name){
  fAMDOption = name;
}

G4String InjectorDetectorConstruction::GetAMDOption(){
  if(fAMDOption) {
    return fAMDOption;
  }else{
    return "";
  }
}

void InjectorDetectorConstruction::SetAMDR1min(G4double value) {
  if(fAMDR1min != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDR1min = value;
  }
}
void InjectorDetectorConstruction::SetAMDR1max(G4double value) {
  if(fAMDR1max != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDR1max = value;
  }
}
void InjectorDetectorConstruction::SetAMDR2min(G4double value) {
  if(fAMDR2min != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDR2min = value;
  }
}
void InjectorDetectorConstruction::SetAMDR2max(G4double value) {
  if(fAMDR2max != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDR2max = value;
  }
}
void InjectorDetectorConstruction::SetAMDLength(G4double value) {
  if(fAMDLength != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLength = value;
  }
  std::cout<<"DEBUG:: AMD length set here to "<<fAMDLength<<" [mm] in SetAMDLength() function.\n"<<std::endl;
}
    G4bool   fAMDLinFringeAct;
    G4double fAMDLinFringeB0;
    G4double fAMDLinFringeK;
    G4double fAMDLinFringeB0ZPos;
    G4double fAMDLinFringeTargZPos;
void InjectorDetectorConstruction::SetAMDLinFringeAct(G4bool value) {
  if(fAMDLinFringeAct != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLinFringeAct = value;
  }
}
void InjectorDetectorConstruction::SetAMDLinFringeB0(G4double value) {
  if(fAMDLinFringeB0 != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLinFringeB0 = value;
  }
}
void InjectorDetectorConstruction::SetAMDLinFringeK(G4double value) {
  if(fAMDLinFringeK != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLinFringeK = value;
  }
}
void InjectorDetectorConstruction::SetAMDLinFringeB0ZPos(G4double value) {
  if(fAMDLinFringeB0ZPos != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLinFringeB0ZPos = value;
  }
}
void InjectorDetectorConstruction::SetAMDLinFringeTargZPos(G4double value) {
  if(fAMDLinFringeTargZPos != value) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fAMDLinFringeTargZPos = value;
  }
}

void InjectorDetectorConstruction::SetDipoleFieldValue(G4ThreeVector field) {
  if(fDipoleFieldValue != field) {
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    fDipoleFieldValue = field;
  }
}

void InjectorDetectorConstruction::FillMagFieldHist(){
  TFile* fFieldMap = new TFile("./FieldHists/B_AMD_2D.root");
  histFieldMapR = (TH2F*) fFieldMap->Get("histFieldMapR");
  histFieldMapZ = (TH2F*) fFieldMap->Get("histFieldMapZ");
}

//void InjectorDetectorConstruction::FillMagFieldGraph(){
//  TFile* fFieldMap = new TFile("./FieldGraphs/B_AMD_2D.root");
//  graphFieldMapR = (TGraph2D*) fFieldMap->Get("g_br");
//  graphFieldMapZ = (TGraph2D*) fFieldMap->Get("g_bz");
//}

