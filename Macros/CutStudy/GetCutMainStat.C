void GetCutMainStat(TString fname="noFile"){
////Usage root 'GetCutStats.C("263180_U30_u370.root")'
  if (fname=="noFile"){
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "Macro used to check cuts after main.C" << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << " " << endl;
    cout << "Usage root:" << endl;
    cout << "\'GetCutStats.C(\"../../output.root\")\'" << endl;
    exit();
  }
  TFile* f = TFile::Open(fname);

  TH1D* h1 = (TH1D*)gROOT->FindObject("hCuts");
  Int_t start_bin = 1;
  Int_t total = h1->GetBinContent(start_bin);

  cout << " " << endl;
  cout << setw(25) << "Total Stats" << setw(10) << total << " " << h1->GetXaxis()->GetBinLabel(start_bin) << endl;
  for(Int_t i=start_bin+10; i<h1->GetXaxis()->GetNbins()+1; i=i+10){
    Double_t percent = h1->GetBinContent(i);
    percent = percent/total;
    percent = 100*percent;
    cout << setw(25) << h1->GetXaxis()->GetBinLabel(i+1) << setw(10) << h1->GetBinContent(i) << " " << percent << endl;
  }

}
