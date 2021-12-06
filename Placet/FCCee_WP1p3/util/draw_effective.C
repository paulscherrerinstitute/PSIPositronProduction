#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

TTree* get_tree(TString datfilename = "output/Beam.dat", TString varname = "B"){

  //datfilename = "Results/CLIC380GeV_Nov2020/Beam.dat";
  //datfilename = "Results/CLIC3TeV_Nov2020/Beam.dat";
  datfilename = "Results/CLIC380GeV_Nov2020_HugoLinAMD/Beam.dat";
  //datfilename = "Results/CLIC380GeV_Nov2020_HugoAMD/Beam.dat";

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
  fin.close();

  return tree;
}


void draw_effective(){

  TTree* tree = get_tree();

  double e1 = 2.86*(1.0-1.2e-2);
  double e2 = 2.86*(1.0+1.2e-2);

  double z1 = -10; // z1
  double z2 = z1+19.8;

  TString e_cut_str = Form("( e_gev > %f && e_gev < %f )", e1, e2);
  TString z_cut_str = Form("( z_um*1e-3 > %f && z_um*1e-3 < %f )", z1, z2);
  TString cut_str = "(" + e_cut_str + "&&" + z_cut_str + ")";

  TH2F* hez = new TH2F("hez", "", 100,-15,25, 100,2.300,3.000);
  TH2F* hez_cut = new TH2F("hez_cut", "", 100,z1,z2, 100,e1,e2);
  tree->Project("hez","e_gev:z_um*1e-3");
  tree->Project("hez_cut","e_gev:z_um*1e-3");

  hez->GetYaxis()->SetTitle("E [GeV]");
  hez->GetXaxis()->SetTitle("Z [mm]");
  hez->GetYaxis()->CenterTitle();
  hez->GetXaxis()->CenterTitle();

  hez_cut->GetYaxis()->SetTitle("E [GeV]");
  hez_cut->GetXaxis()->SetTitle("Z [mm]");
  hez_cut->GetYaxis()->CenterTitle();
  hez_cut->GetXaxis()->CenterTitle();

  gStyle->SetPalette(1); // ncolors=1: Violet -> Red [recommended]
  gStyle->SetOptStat("emrou");

  //TCanvas* can = new TCanvas("can","",10,10,720,620);
  TCanvas* can = new TCanvas("can","",10,10,1520,620);

  can->Divide(2,1);

  can->cd(1);
  hez->Draw("colz");

  TLine line;
  line.SetLineColor(kRed);
  line.DrawLine(z1,e1, z2,e1);
  line.DrawLine(z2,e1, z2,e2);
  line.DrawLine(z2,e2, z1,e2);
  line.DrawLine(z1,e2, z1,e1);

  can->cd(2);
  hez_cut->Draw("colz");

  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gPad->Update();
  gPad->Modified();

  //can->Print("output/phasespace_ez.png");

}
