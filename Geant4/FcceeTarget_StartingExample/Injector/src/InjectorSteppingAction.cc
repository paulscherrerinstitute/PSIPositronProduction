#include "InjectorSteppingAction.hh"

#include "InjectorRunAction.hh"
#include "InjectorDetectorConstruction.hh"
#include "InjectorTrackingAction.hh"
#include "InjectorAnalysis.hh"

#include "G4SteppingManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4UnitsTable.hh"

#include "GlobalVariable.hh"

#include <cmath>

InjectorSteppingAction::InjectorSteppingAction(InjectorDetectorConstruction* det):
  fDetector(det),
  fWorld(0), fXtal(0), fAmorphous(0)
{ 
  first = true;
}

InjectorSteppingAction::~InjectorSteppingAction()
{ }

void InjectorSteppingAction::UserSteppingAction(const G4Step* step)
{
  // get fDetector pointers
  if (first) {
    fWorld   	= fDetector->GetWorld();
    fXtal   	= fDetector->GetXtal();
    fAmorphous 	= fDetector->GetAmorphous();
    first  	= false;
  }

  G4ThreeVector fXtalSize = fDetector->GetXtalSize();
  G4ThreeVector fAmorphousDistance = fDetector->GetAmorphousDistance();
  G4ThreeVector fAmorphousSize = fDetector->GetAmorphousSize();
  G4ThreeVector fDipoleField = fDetector->GetDipoleFieldValue();
  G4ThreeVector fAMDFrontGapSize = fDetector->GetAMDFrontGapSize();
  G4double fAMDLength = fDetector->GetAMDLength();

  bool stop_track_after_xtal = 0;
  bool stop_track_after_amor = 0;
  bool stop_track_before_amd = 0;
  bool stop_track_after_amd = 0;

  if ( fAmorphousSize.z()==0 && fAMDLength==0 ) stop_track_after_xtal = 1;
  if ( fAMDFrontGapSize.z()==0 && fAMDLength==0 ) stop_track_after_amor = 1;
  if ( fAMDLength==0 ) stop_track_before_amd = 1;
  if ( 1 ) stop_track_after_amd = 1;

  // z position of original piont (0,0,0) in G4 world
  G4double world_opz = 0 * CLHEP::meter; // [mm]

  // z positions of special points [mm]
  G4double z_xtal_leave = fXtalSize.z();
  G4double z_amor_arrive = z_xtal_leave + fAmorphousDistance.z();
  G4double z_amor_leave = z_amor_arrive + fAmorphousSize.z();
  G4double z_amd_arrive = z_amor_leave + fAMDFrontGapSize.z();
  G4double z_amd_leave = z_amd_arrive + fAMDLength;

  const G4Event *evt = G4RunManager::GetRunManager()->GetCurrentEvent();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4Track* track = step->GetTrack();

  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint(); 

  G4VPhysicalVolume* volumeP1 = point1->GetPhysicalVolume(); //GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* volumeP2 = point2->GetPhysicalVolume();
  G4String VolNameP1 = "UnknownVolume";
  G4String VolNameP2 = "UnknownVolume";
  if(volumeP1) VolNameP1 = volumeP1->GetName();
  if(volumeP2) VolNameP2 = volumeP2->GetName();

  G4int PDGID = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  bool is_photon   = (PDGID==22);
  bool is_electron = (PDGID==11);
  bool is_positron = (PDGID==-11);

  G4StepStatus statusP1 = point1->GetStepStatus();
  G4StepStatus statusP2 = point2->GetStepStatus();

  G4double stepSize = step->GetTrack()->GetStepLength();

  G4double xP1 = (point1->GetPosition()).x();
  G4double yP1 = (point1->GetPosition()).y();
  G4double zP1 = (point1->GetPosition()).z() - world_opz;
  G4double xP2 = (point2->GetPosition()).x();
  G4double yP2 = (point2->GetPosition()).y();
  G4double zP2 = (point2->GetPosition()).z() - world_opz;

  G4double pxP1 = (point1->GetMomentum()).x();
  G4double pyP1 = (point1->GetMomentum()).y();
  G4double pzP1 = (point1->GetMomentum()).z();
  G4double pxP2 = (point2->GetMomentum()).x();
  G4double pyP2 = (point2->GetMomentum()).y();
  G4double pzP2 = (point2->GetMomentum()).z();

  // For hybrid with very strong field

  bool stop_track_before_amor = 0;
  //if ( fAmorphousDistance.z()>0 && fDipoleField.mag()>10 && (is_electron||is_positron) ) {
    //G4double pP1 = sqrt(pxP1*pxP1+pyP1*pyP1+pzP1*pzP1); // MeV
    //G4double pP2 = sqrt(pxP2*pxP2+pyP2*pyP2+pzP2*pzP2); // MeV
    //if ( pP1<1.0||pP2<1.0 )
  //if ( fAmorphousDistance.z()>0 && (is_electron||is_positron) ) {
  //    stop_track_before_amor = 1;
  //}

  if ( fAmorphousDistance.z()>0 && fDipoleField.mag()>0 && (is_electron||is_positron) )
    if ( zP2 >= z_xtal_leave  && zP2 <= z_amor_arrive )
      if ( pzP2 <= 0 )
        track->SetTrackStatus(fStopAndKill);

  // Flags

  bool is_xtal_leave  = 0;
  if ( fXtalSize.z()>0 && statusP2 == fGeomBoundary && VolNameP1 == "Target" && pzP2 > 0 ){
    //if ( zP1 == z_xtal_leave || zP2 == z_xtal_leave ) is_xtal_leave = 1;
    //if ( std::abs(zP1 - z_xtal_leave) < z_xtal_leave*1e-6 || std::abs(zP2 - z_xtal_leave) < z_xtal_leave*1e-6 ) is_xtal_leave = 1;
    if ( std::abs(zP2 - z_xtal_leave) < z_xtal_leave*1e-6 ) is_xtal_leave = 1;
  }
    //else { //if(is_electron){
    //if(is_xtal_leave && is_electron)
      //if(std::abs(zP2 - z_xtal_leave) > z_xtal_leave*1e-13){
    //else if (fXtalSize.z()>0 && VolNameP1 == "Target" && VolNameP2 == "World"){
    //  G4cout<<"DEBUG:: "<<G4endl;
    //  G4cout<<"       (statusP1 == fGeomBoundary): "<<(statusP1 == fGeomBoundary)<<G4endl;
    //  G4cout<<"       (statusP2 == fGeomBoundary): "<<(statusP2 == fGeomBoundary)<<G4endl;
    //  G4cout<<"       VolNameP1: "<<VolNameP1<<G4endl;
    //  G4cout<<"       VolNameP2: "<<VolNameP2<<G4endl;
    //  G4cout<<"       (pzP1 > 0): "<<(pzP1 > 0)<<G4endl;
    //  G4cout<<"       (pzP2 > 0): "<<(pzP2 > 0)<<G4endl;
    //  G4cout<<"       zP1: "<<zP1<<G4endl;
    //  G4cout<<"       zP2: "<<zP2<<G4endl;
    //  G4cout<<"       xP1, yP1: "<<xP1<<", "<<yP1<<G4endl;
    //  G4cout<<"       xP2, yP2: "<<xP2<<", "<<yP2<<G4endl;
    //  G4cout<<"       pxP1, pyP1, pzP1: "<<pxP1<<", "<<pyP1<<","<<pzP1<<G4endl;
    //  G4cout<<"       pxP2, pyP2, pzP2: "<<pxP2<<", "<<pyP2<<","<<pzP2<<G4endl;
    //  G4cout<<"       z_xtal_leave: "<<z_xtal_leave<<G4endl;
    //  G4cout<<"       zP1-z_xtal_leave: "<<(zP1-z_xtal_leave)<<G4endl;
    //  G4cout<<"       zP2-z_xtal_leave: "<<(zP2-z_xtal_leave)<<G4endl;
    //  G4cout<<"       is_xtal_leave: "<<is_xtal_leave<<G4endl;
    //  //}
    //}

  bool is_amor_arrive = 0;
  if ( fAmorphousDistance.z()>0 && statusP2 == fGeomBoundary && VolNameP1 == "World" && pzP2 > 0 ){
    //if ( zP1 == z_amor_arrive || zP2 == z_amor_arrive ) is_amor_arrive = 1;
    //if ( std::abs(zP1 - z_amor_arrive) < z_amor_arrive*1e-6 || std::abs(zP2 - z_amor_arrive) < z_amor_arrive*1e-6 ) is_amor_arrive = 1;
    if ( std::abs(zP2 - z_amor_arrive) < z_amor_arrive*1e-6 ) is_amor_arrive = 1;
  }

  bool is_amor_leave  = 0;
  if ( fAmorphousSize.z()>0 && statusP2 == fGeomBoundary && VolNameP1 == "Amorphous" && pzP2 > 0 ){
    //if ( zP1 == z_amor_leave || zP2 == z_amor_leave ) is_amor_leave = 1;
    //if ( std::abs(zP1 - z_amor_leave) < z_amor_leave*1e-6 || std::abs(zP2 - z_amor_leave) < z_amor_leave*1e-6 ) is_amor_leave = 1;
    if ( std::abs(zP2 - z_amor_leave) < z_amor_leave*1e-6 ) is_amor_leave = 1;
    //else if(is_positron){
    //if(is_amor_leave && is_positron)
    // if(std::abs(zP2 - z_amor_leave) > z_amor_leave*1e-13){
    //  G4cout<<"DEBUG:: "<<G4endl;
    //  G4cout<<"       (statusP2 == fGeomBoundary): "<<(statusP2 == fGeomBoundary)<<G4endl;
    //  G4cout<<"       VolNameP1: "<<VolNameP1<<G4endl;
    //  G4cout<<"       VolNameP2: "<<VolNameP2<<G4endl;
    //  G4cout<<"       (pzP1 > 0): "<<(pzP1 > 0)<<G4endl;
    //  G4cout<<"       (pzP2 > 0): "<<(pzP2 > 0)<<G4endl;
    //  G4cout<<"       zP1: "<<zP1<<G4endl;
    //  G4cout<<"       zP2: "<<zP2<<G4endl;
    //  G4cout<<"       xP1, yP1: "<<xP1<<", "<<yP1<<G4endl;
    //  G4cout<<"       xP2, yP2: "<<xP2<<", "<<yP2<<G4endl;
    //  G4cout<<"       pxP1, pyP1, pzP1: "<<pxP1<<", "<<pyP1<<","<<pzP1<<G4endl;
    //  G4cout<<"       pxP2, pyP2, pzP2: "<<pxP2<<", "<<pyP2<<","<<pzP2<<G4endl;
    //  G4cout<<"       z_amor_leave: "<<z_amor_leave<<G4endl;
    //  G4cout<<"       zP1-z_amor_leave: "<<(zP1-z_amor_leave)<<G4endl;
    //  G4cout<<"       zP2-z_amor_leave: "<<(zP2-z_amor_leave)<<G4endl;
    //  G4cout<<"       is_amor_leave: "<<is_amor_leave<<G4endl;
    //}
  }

  bool is_amd_arrive  = 0;
  if ( fAMDFrontGapSize.z()>0 && statusP2 == fGeomBoundary && VolNameP1 == "AMDFrontGap" && pzP2 > 0 ){
    //if ( zP1 == z_amd_arrive || zP2 == z_amd_arrive ) is_amd_arrive = 1;
    //if ( std::abs(zP1 - z_amd_arrive) < z_amd_arrive*1e-6 || std::abs(zP2 - z_amd_arrive) < z_amd_arrive*1e-6 ) is_amd_arrive = 1;
    if ( std::abs(zP2 - z_amd_arrive) < z_amd_arrive*1e-6 ) is_amd_arrive = 1;
  }

  bool is_amd_leave  = 0;
  if ( fAMDLength>0 && statusP2 == fGeomBoundary && VolNameP1 == "AMDVacPhysicalLast" && pzP2 > 0 ){
    //if ( zP1 == z_amd_leave || zP2 == z_amd_leave ) is_amd_leave = 1;
    //if ( std::abs(zP1 - z_amd_leave) < z_amd_leave*1e-6 || std::abs(zP2 - z_amd_leave) < z_amd_leave*1e-6 ) is_amd_leave = 1;
    if ( std::abs(zP2 - z_amd_leave) < z_amd_leave*1e-6 ) is_amd_leave = 1;
  }

  // Debug

  if(0 && statusP2 == fGeomBoundary){
    G4cout<<"INFO:: "<<G4endl;
    G4cout<<"       fXtalSize.z(): "<<(fXtalSize.z())<<G4endl;
    G4cout<<"       fAmorphousDistance.z(): "<<(fAmorphousDistance.z())<<G4endl;
    G4cout<<"       fAmorphousSize.z(): "<<(fAmorphousSize.z())<<G4endl;
    G4cout<<"       (statusP2 == fGeomBoundary): "<<(statusP2 == fGeomBoundary)<<G4endl;
    G4cout<<"       VolNameP1: "<<VolNameP1<<G4endl;
    G4cout<<"       VolNameP2: "<<VolNameP2<<G4endl;
    G4cout<<"       (pzP2 > 0): "<<(pzP2 > 0)<<G4endl;
    G4cout<<"       zP1: "<<zP1<<G4endl;
    G4cout<<"       zP2: "<<zP2<<G4endl;
    G4cout<<"       is_amor_leave: "<<is_amor_leave<<G4endl;
  }
   
  // Fill ntuples

  // particles leaving from xtal target
  //if (statusP2 == fGeomBoundary && VolNameP1 == "Target" && VolNameP2 == "World" && pzP2 > 0){
  if (is_xtal_leave){
    if (is_photon||is_positron||is_electron){

      // primary injected electrons not counted here, but considered in PrimaryGeneration class

      G4ThreeVector X  = point2->GetPosition();
      G4ThreeVector P  = point2->GetMomentum();
      G4double e    = point2->GetTotalEnergy();
      G4double tG = point2->GetGlobalTime() / (CLHEP::millimeter / CLHEP::c_light); // since the event is created
      G4double tL = point2->GetLocalTime() / (CLHEP::millimeter / CLHEP::c_light);  // since the track is created
      G4double x = X.x();
      G4double y = X.y();
      G4double z = X.z() - world_opz;
      G4double px = P.x();
      G4double py = P.y();
      G4double pz = P.z();

      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="xtal_leave") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(varname=="x")  value = x; // position photon leaving from xtal [mm]
          if(varname=="y")  value = y;
          if(varname=="z")  value = z;
          if(varname=="px") value = px; // photon momentum [MeV]
          if(varname=="py") value = py;
          if(varname=="pz") value = pz;
          if(varname=="e")  value = e;  // energy [MeV]
          if(varname=="t")  value = tG; // time photon leaving from xtal [mm/c]
          if(varname=="pdgId"){ 
	    if(is_photon)   value = 22;
	    if(is_positron) value = -11;
	    if(is_electron) value = 11;
	  }

          if(varname=="evtId") value = evt->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples
    } // end if of particles
    // kill all particles to avoid double-counting due to magnetic field
    if(stop_track_after_xtal) track->SetTrackStatus(fStopAndKill);
  } // end if of leaving xtal

  // particles arriving at amorphous target
  //if (statusP2 == fGeomBoundary && VolNameP1 == "World" && VolNameP2 == "Amorphous" && pzP2 > 0){
  if (is_amor_arrive){
    if (is_photon||is_electron||is_positron){

      G4ThreeVector X  = point2->GetPosition();
      G4ThreeVector P  = point2->GetMomentum();
      G4double e    = point2->GetTotalEnergy();
      G4double tG = point2->GetGlobalTime() / (CLHEP::millimeter / CLHEP::c_light); // since the event is created
      G4double tL = point2->GetLocalTime() / (CLHEP::millimeter / CLHEP::c_light);  // since the track is created
      G4double x = X.x();
      G4double y = X.y();
      G4double z = X.z() - world_opz;
      G4double px = P.x();
      G4double py = P.y();
      G4double pz = P.z();

      // fill my ntuples
      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="amor_arrive") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(varname=="x")  value = x; // position photon arriving at amorphous [mm]
          if(varname=="y")  value = y;
          if(varname=="z")  value = z;
          if(varname=="px") value = px; // photon momentum [MeV]
          if(varname=="py") value = py;
          if(varname=="pz") value = pz;
          if(varname=="e")  value = e;  // energy [MeV]
          if(varname=="t")  value = tG; // time photon arriving at amorphous [mm/c]
          if(varname=="pdgId"){ 
	    if(is_photon)   value = 22;
	    if(is_electron) value = 11;
	    if(is_positron) value = -11;
	  }

          if(varname=="evtId") value = evt->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples
    } // end if of particles
    if(stop_track_before_amor) track->SetTrackStatus(fStopAndKill);
  } // end if of arriving amorphous

  // particles leaving from amorphous target
  //if (statusP2 == fGeomBoundary && VolNameP1 == "Amorphous" && VolNameP2 == "World" && pzP2 > 0){
  if (is_amor_leave){
    if ( is_positron || ((is_photon||is_electron) && vNtupName.size()!=1) ) { // fill only positrons when doing optimisation

      G4ThreeVector X  = point2->GetPosition();
      G4ThreeVector P  = point2->GetMomentum();
      G4double e    = point2->GetTotalEnergy();
      G4double tG = point2->GetGlobalTime() / (CLHEP::millimeter / CLHEP::c_light); // since the event is created
      G4double tL = point2->GetLocalTime() / (CLHEP::millimeter / CLHEP::c_light);  // since the track is created
      G4double x = X.x();
      G4double y = X.y();
      G4double z = X.z() - world_opz;
      G4double px = P.x();
      G4double py = P.y();
      G4double pz = P.z();

      // fill my ntuples
      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="amor_leave") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(varname=="x")  value = x; // position positron leaving from amorphous [mm]
          if(varname=="y")  value = y;
          if(varname=="z")  value = z;
          if(varname=="px") value = px; // positron momentum [MeV]
          if(varname=="py") value = py;
          if(varname=="pz") value = pz;
          if(varname=="e")  value = e;  // energy [MeV]
          if(varname=="t")  value = tG; // time positron leaving from amorphous [mm/c]
          if(varname=="pdgId"){ 
	    if(is_photon)   value = 22;
	    if(is_electron) value = 11;
	    if(is_positron) value = -11;
	  }

          if(varname=="evtId") value = evt->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples
    } // end if of particles
    // kill all particles to avoid double-counting due to magnetic field
    //G4ThreeVector v3_amd_frontgap_size = fDetector->GetAMDFrontGapSize();
    if(stop_track_after_amor) track->SetTrackStatus(fStopAndKill);
  } // end if of leaving amorphous

  // particles leaving from gap / arriving at amd
  if (is_amd_arrive){
    if ( is_positron || ((is_photon||is_electron) && vNtupName.size()!=1) ) { // fill only positrons when doing optimisation

      G4ThreeVector X  = point2->GetPosition();
      G4ThreeVector P  = point2->GetMomentum();
      G4double e    = point2->GetTotalEnergy();
      G4double tG = point2->GetGlobalTime() / (CLHEP::millimeter / CLHEP::c_light); // since the event is created
      G4double tL = point2->GetLocalTime() / (CLHEP::millimeter / CLHEP::c_light);  // since the track is created
      G4double x = X.x();
      G4double y = X.y();
      G4double z = X.z() - world_opz;
      G4double px = P.x();
      G4double py = P.y();
      G4double pz = P.z();

      // fill my ntuples
      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="amd_arrive") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(varname=="x")  value = x; // position positron leaving from gap / arriving at amd [mm]
          if(varname=="y")  value = y;
          if(varname=="z")  value = z;
          if(varname=="px") value = px; // positron momentum [MeV]
          if(varname=="py") value = py;
          if(varname=="pz") value = pz;
          if(varname=="e")  value = e;  // energy [MeV]
          if(varname=="t")  value = tG; // time positron leaving from amorphous [mm/c]
          if(varname=="pdgId"){ 
	    if(is_photon)   value = 22;
	    if(is_electron) value = 11;
	    if(is_positron) value = -11;
	  }

          if(varname=="evtId") value = evt->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples
    } // end if of particles
    // kill all particles to avoid double-counting due to magnetic field
    if(stop_track_before_amd) track->SetTrackStatus(fStopAndKill);
  } // end if of leaving gap_targ_amd / arriving at amd

  // particles leaving amd
  if (is_amd_leave){
    if ( is_positron || ((is_photon||is_electron) && vNtupName.size()!=1) ) { // fill only positrons when doing optimisation

      G4ThreeVector X  = point2->GetPosition();
      G4ThreeVector P  = point2->GetMomentum();
      G4double e    = point2->GetTotalEnergy();
      G4double tG = point2->GetGlobalTime() / (CLHEP::millimeter / CLHEP::c_light); // since the event is created
      G4double tL = point2->GetLocalTime() / (CLHEP::millimeter / CLHEP::c_light);  // since the track is created
      G4double x = X.x();
      G4double y = X.y();
      G4double z = X.z() - world_opz;
      G4double px = P.x();
      G4double py = P.y();
      G4double pz = P.z();

      // fill my ntuples
      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="amd_leave") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(varname=="x")  value = x; // position positron leaving from gap / arriving at amd [mm]
          if(varname=="y")  value = y;
          if(varname=="z")  value = z;
          if(varname=="px") value = px; // positron momentum [MeV]
          if(varname=="py") value = py;
          if(varname=="pz") value = pz;
          if(varname=="e")  value = e;  // energy [MeV]
          if(varname=="t")  value = tG; // time positron leaving from amorphous [mm/c]
          if(varname=="pdgId"){ 
	    if(is_photon)   value = 22;
	    if(is_electron) value = 11;
	    if(is_positron) value = -11;
	  }

          if(varname=="evtId") value = evt->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples
    } // end if of particles
    // kill all particles to avoid double-counting due to magnetic field
    if(stop_track_after_amd) track->SetTrackStatus(fStopAndKill);
  } // end if of leaving amd
}

