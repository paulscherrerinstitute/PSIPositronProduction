#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

void get_root_beam_PLC(TString datfilename = "output/Beam_All.dat", TString varname = "B"){

  TString pathnameDir = datfilename;
  TString filename = datfilename;

  pathnameDir.Remove(pathnameDir.Last('/'), pathnameDir.Length() - pathnameDir.Last('/'));
  filename.Remove(0, filename.Last('/') + 1);

  TString pathnameRoot = pathnameDir;
  //pathnameRoot.ReplaceAll("/Dat","/Root");

  TString rootfilename = filename;
  rootfilename.ReplaceAll(".dat",".root");
  //if(varname=="A_TW1")
  //rootfilename.ReplaceAll("TW_","TW1_");

  cout<<"Converting "<<pathnameDir<<"/"<<filename<<" to "<<pathnameRoot<<"/"<<rootfilename<<" ..."<<endl;

  TTree* tree = new TTree(varname,varname);

  double e_gev = 0;
  double x_um = 0;
  double y_um = 0;
  double z_um = 0;
  double xp_urad = 0;
  double yp_urad = 0;

  tree->Branch("e_gev",&e_gev);
  tree->Branch("x_um",&x_um);
  tree->Branch("y_um",&y_um);
  tree->Branch("z_um",&z_um);
  tree->Branch("xp_urad",&xp_urad);
  tree->Branch("yp_urad",&yp_urad);

  ifstream fin(datfilename);

  string  line;
  TString Line;

  while(getline(fin,line)){
    Line = line;
    // find matrix
    if(Line != "# name: "+varname){
      continue;
    }else{
      getline(fin,line); // type
      getline(fin,line); // rows
      Line = line;
      Line.ReplaceAll("# rows: ","");
      int NR = Line.Atoi();
      getline(fin,line); // columns
      Line = line;
      Line.ReplaceAll("# columns: ","");
      int NC = Line.Atoi();
      for(int i=0;i<NR;i++){
        getline(fin,line); Line = line;
	TObjArray* TOA_Word = Line.Tokenize(" ");
        for(int j=0;j<NC;j++){
	  TString Word = ((TObjString*)(TOA_Word->At(j)))->String();
	  double value = Word.Atof();
	  if(j+1==1) e_gev = value;
	  if(j+1==2) x_um = value;
	  if(j+1==3) y_um = value;
	  if(j+1==4) z_um = value;
	  if(j+1==5) xp_urad = value;
	  if(j+1==6) yp_urad = value;
	}

	tree->Fill();
      }
      break;
    }
  }

  TFile* file = new TFile(pathnameRoot+"/"+rootfilename, "RECREATE");

  tree->Write();
  file->Close();
}

