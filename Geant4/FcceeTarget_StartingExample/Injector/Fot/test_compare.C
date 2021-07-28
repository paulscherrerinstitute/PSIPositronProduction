#include "TApplication.h"
#include "TEntryList.h"
#include "TEventList.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TCut.h"
#include "TFile.h"
#include "TSystem.h"

TChain *chain1;
TChain *chain2;

void file(){
  chain1 = new TChain("tree");
  chain2 = new TChain("tree");

  chain1->Add("output_Fot.root");
  chain2->Add("output_Fot_old.root");
}

void px(){
  file();
  TPostScript ps("px.eps",113);
  TH1F *histo1 = new TH1F("histo1",";Px;Entries",100,-0.001,0.001);
  TH1F *histo2 = new TH1F("histo2",";Px;Entries",100,-0.001,0.001);
  chain1->Draw("px >> histo1");
  chain2->Draw("px >> histo2");
  histo1->Draw();
  histo2->Draw("esame");
  ps.Close();
}

void py(){
  file();
  TPostScript ps("py.eps",113);
  TH1F *histo1 = new TH1F("histo1",";Py;Entries",100,-0.001,0.001);
  TH1F *histo2 = new TH1F("histo2",";Py;Entries",100,-0.001,0.001);
  chain1->Draw("py >> histo1");
  chain2->Draw("py >> histo2");
  histo1->Draw();
  histo2->Draw("esame");
  ps.Close();
}


void energy(){
  file();
  TPostScript ps("e.eps",113);
  TH1F *histo1 = new TH1F("histo1",";Energy;Entries",100,0,0.2);
  TH1F *histo2 = new TH1F("histo2",";Energy;Entries",100,0,0.2);
  chain1->Draw("e >> histo1");
  chain2->Draw("e >> histo2");
  histo1->Draw();
  histo2->Draw("esame");
  ps.Close();
}

void x(){
  file();
  TPostScript ps("x.eps",113);
  TH1F *histo1 = new TH1F("histo1",";x;Entries",100,-100e6,100e6);
  TH1F *histo2 = new TH1F("histo2",";x;Entries",100,-100e6,100e6);
  chain1->Draw("x >> histo1");
  chain2->Draw("x >> histo2");
  histo1->Draw();
  histo2->Draw("esame");
  ps.Close();

}


void y(){
  file();
  TPostScript ps("y.eps",113);
  TH1F *histo1 = new TH1F("histo1",";y;Entries",100,-100e6,100e6);
  TH1F *histo2 = new TH1F("histo2",";y;Entries",100,-100e6,100e6);
  chain1->Draw("y >> histo1");
  chain2->Draw("y >> histo2");
  histo1->Draw();
  histo2->Draw("esame");
  ps.Close();

}

