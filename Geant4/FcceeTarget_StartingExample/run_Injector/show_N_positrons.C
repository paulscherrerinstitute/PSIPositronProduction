{
  TFile f("FCCeeTargetTracking.root");
  TTree* t = (TTree*) f.Get("amor_leave"); // at target exit
  int np = t->GetEntries("pdgId==-11"); // number of e+
  cout<<"Number of e+ at target exit: "<<np<<endl;
}
