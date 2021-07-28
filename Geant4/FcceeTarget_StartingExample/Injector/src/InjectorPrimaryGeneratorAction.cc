#include "InjectorPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"

#include "InjectorPrimaryGeneratorActionMessenger.hh"
#include "InjectorAnalysis.hh"
#include "InjectorDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include "Fot.h"
#include "Particle.h"
#include "Crystal.h"
#include "PhotonCollection.h"
#include "RunParameters.h"

#include "GlobalVariable.hh"

#include <iostream>
#include <list>
#include <math.h>
#include <iomanip> // std::setprecision

#include "G4PhysicalConstants.hh"

using namespace std;

//InjectorPrimaryGeneratorAction::InjectorPrimaryGeneratorAction():
InjectorPrimaryGeneratorAction::InjectorPrimaryGeneratorAction(InjectorDetectorConstruction* det) : fDet(det)
{
    fParticleGun = new G4ParticleGun();
    fGunMessenger = new InjectorPrimaryGeneratorActionMessenger(this);

    G4ThreeVector xtalSize = fDet->GetXtalSize();
    XtalThick = xtalSize.z(); // [mm]

    primaryEnergy = 3000; // [MeV]
    sigmaUserE = 3.0;     // [MeV]
    sigmaUserXY = 2.5;    // [mm]
    sigmaUserZ = 0.3;     // [mm]
    sigmaUserT = 0.3;     // [mm/c]
    sigmaUserPxy = 0.01;  // [MeV]
}

InjectorPrimaryGeneratorAction::~InjectorPrimaryGeneratorAction(){
    delete fParticleGun;
}

void InjectorPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

  G4double Energy 	= this->primaryEnergy; // [MeV]
  G4double sigmaE 	= this->sigmaUserE;    // [MeV]
  G4double xtal_thick 	= this->XtalThick;     // [mm]
  G4double sigmaXY   	= this->sigmaUserXY;   // [mm]
  G4double sigmaT 	= this->sigmaUserT;    // [[mm/c]]
  G4double sigmaPxy 	= this->sigmaUserPxy;  // [MeV]

  // z position of original piont (0,0,0) in G4 world
  G4double world_opz = 0 * CLHEP::meter; // [mm]

  G4double e = G4RandGauss::shoot(Energy, sigmaE); // [MeV]

  G4double x = G4RandGauss::shoot(0, sigmaXY); // [mm]
  G4double y = G4RandGauss::shoot(0, sigmaXY); // [mm]
  G4double z = 0;

  G4double px = G4RandGauss::shoot(0, sigmaPxy); // [MeV]
  G4double py = G4RandGauss::shoot(0, sigmaPxy); // [MeV]
  //G4double pz = sqrt( pow(e,2) - pow(CLHEP::electron_mass_c2,2) - pow(px,2) - pow(py,2) ); // [MeV]
  G4double pz2 = pow(e,2) - pow(CLHEP::electron_mass_c2,2) - pow(px,2) - pow(py,2); // [MeV^2]
  G4double pz = 0;
  if( pz2 > 0 ){ 
    pz = sqrt(pz2);
  }else{
    e = sqrt(pow(CLHEP::electron_mass_c2,2) + pow(px,2) + pow(py,2));
  }
  G4double p  = sqrt( pow(px,2) + pow(py,2) + pow(pz,2) ); // [MeV]

  G4double t = G4RandGauss::shoot(0, sigmaT);  // [mm/c]

  // create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill ntuples
  for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
    const G4String ntup_name  = vNtupName[i_ntup];
    if(ntup_name!="primary") continue;
    const G4String ntup_title = vNtupTitle[i_ntup];
    for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
      G4String varname = vVarToFill[ivar];
      G4double value = 0;

      if(ntup_name=="primary"){
        if(varname=="e")  value = e; // [MeV]
        if(varname=="x")  value = x; // [mm]
        if(varname=="y")  value = y;
        if(varname=="z")  value = z;
        if(varname=="px") value = px; // [MeV]
        if(varname=="py") value = py;
        if(varname=="pz") value = pz;
        if(varname=="t")  value = t;  // [mm/c]
        if(varname=="pdgId") value = 11;
      } 

      if(varname=="evtId") value = anEvent->GetEventID();

      if(varname=="pdgId"||varname=="evtId")
        analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
      else
        analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
    } // end loop of variables
    analysisManager->AddNtupleRow(i_ntup);
  } // end loop of ntuples

  if(xtal_thick!=0){

    // channeling process using Fot
    //G4cout<<"DEBUG:: Simulating channeling process with Fot. Xtal_thick = "<<xtal_thick<<" mm."<<G4endl;

    // crystal tungsten target
    Crystal crysW("W",111);
    RunParameters rp(crysW);

    // set exit z position of electron from crystal target
    rp.setZexit(xtal_thick / CLHEP::angstrom);

    // PARAMETRES DU RUN
    G4double phomin	= 2e-3;   // min energy of photon
    G4double etmax 	= 1e+5;
    G4double vtmax 	= 1e-2;
    G4double poimin	= 1;      // weight of crystal
    G4double frekuma 	= 3./259; // testing frequency of KUMA
    rp.set(phomin, etmax, vtmax, poimin, frekuma);

    Fot fot(rp);

    double charge  = -1.0;
    ParticleInCrystal partCrys = fot.makeSingleParticleKumakhov( 
      new Particle( charge, // [e]
                    x  / CLHEP::angstrom, // [angstrom]
                    y  / CLHEP::angstrom,
          	  z  / CLHEP::angstrom,
          	  px / CLHEP::electron_mass_c2, // [me^-1]
          	  py / CLHEP::electron_mass_c2,
          	  e  / CLHEP::electron_mass_c2 ));

    // emit photons provided by Fot using Geant4::ParticleGun

    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particle);

    list<Photon> list_photon_kumakov = fot.getPhotonCollection().getPhotonList();

    G4double RemainedEnergy = e; // electron remained energy after channeling

    list<Photon >::iterator it;
    for(it = list_photon_kumakov.begin() ; it != list_photon_kumakov.end() ; it++)
    {
      G4double e_photon = it->getEnergy() * CLHEP::GeV; // photon energy [MeV]

      // photon momentum
      G4double p_gamma  = e_photon; // [MeV]
      G4double px_gamma = (it->getThetax() * p_gamma); // [MeV]
      G4double py_gamma = (it->getThetay() * p_gamma); // [MeV]
      G4double pz_gamma = sqrt(pow(p_gamma,2) - pow(px_gamma,2) - pow(py_gamma,2)); // [MeV]

      // position and moving direction (unit vector) of photon when emitted
      G4ThreeVector world_position = G4ThreeVector( it->getXemis() * CLHEP::angstrom + x,
                                                    it->getYemis() * CLHEP::angstrom + y,
          			                  it->getZemis() * CLHEP::angstrom + z + world_opz);
      G4ThreeVector direction_emission = G4ThreeVector( px_gamma / p_gamma,
                                                        py_gamma / p_gamma,
          			                      pz_gamma / p_gamma);

      RemainedEnergy -= e_photon;

      // position and time when photons emitted
      G4double x_emit = world_position.x(); // [mm]
      G4double y_emit = world_position.y(); // [mm]
      G4double z_emit = world_position.z() - world_opz; // [mm]
      G4double t_emit = t + z_emit; // [mm/c]

      // NB. photons leaving from crystal target will be filled in SteppingAction

      // fill my ntuples
      for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
        const G4String ntup_name  = vNtupName[i_ntup];
        if(ntup_name!="photon_emit") continue;
        const G4String ntup_title = vNtupTitle[i_ntup];
        for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
          G4String varname = vVarToFill[ivar];
          G4double value = 0;

          if(ntup_name=="photon_emit"){
            if(varname=="x")  value = x_emit; // [mm]
            if(varname=="y")  value = y_emit;
            if(varname=="z")  value = z_emit;
            if(varname=="px") value = px_gamma; // [MeV]
            if(varname=="py") value = py_gamma;
            if(varname=="pz") value = pz_gamma;
            if(varname=="e")  value = e_photon;
            if(varname=="t")  value = t_emit; // [mm/c]
            if(varname=="pdgId") value = 22;
          }

          if(varname=="evtId") value = anEvent->GetEventID();

          if(varname=="pdgId"||varname=="evtId")
            analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
          else
            analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
        } // end loop of variables
        analysisManager->AddNtupleRow(i_ntup);
      } // end loop of ntuples

      // set Geant4::ParticleGun
      fParticleGun->SetParticleEnergy(e_photon);
      fParticleGun->SetParticlePosition(world_position);
      fParticleGun->SetParticleMomentumDirection(direction_emission);
      fParticleGun->SetParticleTime(t_emit * CLHEP::millimeter / CLHEP::c_light);
      fParticleGun->GeneratePrimaryVertex(anEvent);
    } 

    // simulate electrons provided by Fot using Geant4::ParticleGun

    particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particle);

    G4double ze  = rp.getZexit() * CLHEP::angstrom; // [mm]
    G4double zpe = partCrys.getZPosition() * CLHEP::angstrom;  // [mm]
    G4double pxe = CLHEP::electron_mass_c2 * partCrys.getPx(); // [MeV]
    G4double pye = CLHEP::electron_mass_c2 * partCrys.getPy(); // [MeV]
    G4double ge  = partCrys.getGamma(); // [-]
    G4double pze = sqrt( (ge*ge - 1) * pow(CLHEP::electron_mass_c2, 2) - pow(pxe,2) - pow(pye,2) ); // [MeV]
    G4double pe  = sqrt( pow(pxe,2) + pow(pye,2) + pow(pze,2) ); // [MeV]

    if(0 || isnan(pze)){ // DEBUG
     if(isnan(pze)){
      G4cout<<"ERROR__:: pze: "<<pze<<" <0, pxe: "<<pxe<<", pye: "<<pye<<", ge: "<<ge<<G4endl;
      G4cout<<"          pe: "<<pe<<G4endl;
      G4cout<<"          e: "<<RemainedEnergy<<G4endl;
     }
     if(1){
      G4cout<<"xtal_arrive (electron): "<<G4endl;
      G4cout<<"    x = "<<x<<";// [mm]"<<G4endl;
      G4cout<<"    y = "<<y<<";// [mm]"<<G4endl;
      G4cout<<"    z = 0;// [mm]"<<G4endl;
      G4cout<<"    px = "<<px<<";// [MeV]"<<G4endl;
      G4cout<<"    py = "<<py<<";// [MeV]"<<G4endl;
      G4cout<<"    pz = "<<pz<<";// [MeV]"<<G4endl;
      G4cout<<"    e  =  "<<e<<";// [MeV]"<<G4endl;
      G4cout<<"number of photons: "<<(list_photon_kumakov.size())<<G4endl;
      list<Photon >::iterator it_tmp;
      double e_photon_total = 0;
      for(it_tmp = list_photon_kumakov.begin() ; it_tmp != list_photon_kumakov.end() ; it_tmp++){
        double e_photon = it_tmp->getEnergy() * CLHEP::GeV; // photon energy [MeV]
        G4cout<<"    photon energy: "<<e_photon<<" [MeV]"<<endl;
        e_photon_total += e_photon;
      }
      G4cout<<"    total photon energy: "<<e_photon_total<<G4endl;
      G4cout<<"xtal_leave (electron): "<<G4endl;
      G4cout<<"    x = "<<(partCrys.getXPosition()*CLHEP::angstrom)<<";// [mm]"<<G4endl;
      G4cout<<"    y = "<<(partCrys.getYPosition()*CLHEP::angstrom)<<";// [mm]"<<G4endl;
      G4cout<<"    z_p (partCrys.getZPosition) = "<<zpe<<";// [mm]"<<G4endl;
      G4cout<<"    z_e (rp.getZexit)           = "<<ze<<";// [mm]"<<G4endl;
      G4cout<<"    px = "<<pxe<<";// [MeV]"<<G4endl;
      G4cout<<"    py = "<<pye<<";// [MeV]"<<G4endl;
      G4cout<<"    pz = "<<pze<<";// [MeV]"<<G4endl;
      G4cout<<"    e  =  "<<RemainedEnergy<<";// [MeV]"<<G4endl;
     }
    }

    G4double kineticEnergy = RemainedEnergy - CLHEP::electron_mass_c2; // [MeV]

    G4ThreeVector world_position  = G4ThreeVector( partCrys.getXPosition() * CLHEP::angstrom + x,
                                                   partCrys.getYPosition() * CLHEP::angstrom + y,
          					 ze + z + world_opz);
    G4ThreeVector direction_emission = G4ThreeVector( pxe/pe, pye/pe, pze/pe);

    // electrons: position and time when leaving from crystal target
    G4double xtal_leave_t = t + ze; // [mm/c]

    // fill my ntuples
    for(unsigned int i_ntup=0;i_ntup<vNtupName.size();i_ntup++){
      const G4String ntup_name  = vNtupName[i_ntup];
      if(ntup_name!="xtal_leave") continue;
      const G4String ntup_title = vNtupTitle[i_ntup];
      for(unsigned int ivar=0;ivar<vVarToFill.size();ivar++){
        G4String varname = vVarToFill[ivar];
        G4double value = 0;

        if(ntup_name=="xtal_leave"){
          if(varname=="x")  value = world_position.x(); // [mm]
          if(varname=="y")  value = world_position.y(); // [mm]
          if(varname=="z")  value = world_position.z() - world_opz; // [mm]
          if(varname=="px") value = direction_emission.x() * pe; // [MeV]
          if(varname=="py") value = direction_emission.y() * pe; // [MeV]
          if(varname=="pz") value = direction_emission.z() * pe; // [MeV]
          if(varname=="e")  value = RemainedEnergy; // [MeV]
          if(varname=="t")  value = xtal_leave_t;   // [mm/c]
          if(varname=="pdgId") value = 11;
        }

        if(varname=="evtId") value = anEvent->GetEventID();

        if(varname=="pdgId"||varname=="evtId")
          analysisManager->FillNtupleIColumn(i_ntup, ivar, value);
        else
          analysisManager->FillNtupleDColumn(i_ntup, ivar, value);
      } // end loop of variables
      analysisManager->AddNtupleRow(i_ntup);
    } // end loop of ntuples

    fParticleGun->SetParticleEnergy(kineticEnergy);
    fParticleGun->SetParticlePosition(world_position);
    fParticleGun->SetParticleMomentumDirection(direction_emission);
    fParticleGun->SetParticleTime(xtal_leave_t * CLHEP::millimeter / CLHEP::c_light);
    fParticleGun->GeneratePrimaryVertex(anEvent);

  } else {
    G4ThreeVector world_position     = G4ThreeVector(x, y, z + world_opz);
    G4ThreeVector direction_emission = G4ThreeVector(px/p, py/p, pz/p);

    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particle);
    G4double kineticEnergy = e - CLHEP::electron_mass_c2; // [MeV]
    fParticleGun->SetParticleEnergy(kineticEnergy);
    fParticleGun->SetParticlePosition(world_position);
    fParticleGun->SetParticleMomentumDirection(direction_emission);
    //fParticleGun->SetParticleTime(0);
    fParticleGun->SetParticleTime(t * CLHEP::millimeter / CLHEP::c_light);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//void InjectorPrimaryGeneratorAction::SetXtalThick(G4double value){
//  XtalThick = value;
//}

void InjectorPrimaryGeneratorAction::SetPrimaryEnergy(G4double value){
  primaryEnergy = value;
}

void InjectorPrimaryGeneratorAction::SetSigmaE(G4double value){
  sigmaUserE = value;
}

void InjectorPrimaryGeneratorAction::SetSigmaXY(G4double value){
  sigmaUserXY = value;
}

void InjectorPrimaryGeneratorAction::SetSigmaZ(G4double value){
    sigmaUserZ = value;
}

void InjectorPrimaryGeneratorAction::SetSigmaT(G4double value){
    sigmaUserT = value;
}

void InjectorPrimaryGeneratorAction::SetSigmaPxy(G4double value){
    sigmaUserPxy = value;
}

G4double InjectorPrimaryGeneratorAction::GetVelocity(G4double p, G4double m){
  // p in [MeV/c]; m in [MeV/c/c]
  return (p/m)*CLHEP::c_light / sqrt(1+pow(p/m, 2));
}
