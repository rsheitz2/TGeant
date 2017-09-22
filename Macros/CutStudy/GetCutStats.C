void GetCutStats(TString fname="noFile", TString UserNumber="410"){
////Usage root 'GetCutStats.C("263180_U30_u370.root")'
  if (fname=="noFile"){
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "Macro used to check cuts after UserEvent in phast" << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << " " << endl;
    cout << "Usage root:" << endl;
    cout << "\'GetCutStats.C(\"263180_U30_u370.root\")\'" << endl;
    cout << " " << endl;
    cout << "Can also specify another UserEvent number (default is 410)" <<
      endl;
    cout << "Example:" << endl;
    cout << "\'GetCutStats.C(\"263180_U30_u370.root\", \"410\")\'" << endl;
    exit();
  }
  
  TFile* f = TFile::Open(fname);
  TString UserEvent = "UserEvent" + UserNumber;
  cout << UserEvent << endl;

  f->cd(UserEvent);


  TH1D* h1 = (TH1D*)gROOT->FindObject("DiMuonCuts");
  //Int_t start_bin = 41;
  Int_t start_bin = 1;
  Int_t total = h1->GetBinContent(start_bin);

  cout << " " << endl;
  cout << setw(25) << "Total Stats" << setw(10) << total << " " << h1->GetXaxis()->GetBinLabel(start_bin) << endl;
  for(Int_t i=start_bin+10; i<h1->GetXaxis()->GetNbins()+1; i=i+10){
    Double_t percent = h1->GetBinContent(i);
    percent = percent/total;
    percent = 100*percent;
    cout << setw(25) << h1->GetXaxis()->GetBinLabel(i) << setw(10) << h1->GetBinContent(i) << " " << percent << endl;
  }

  TH1D* hHighM = (TH1D*)gROOT->FindObject("HighMassDiMuonCuts");
  Int_t total_highM = hHighM->GetBinContent(start_bin);
  cout << " " << endl;
  cout << setw(25) << "Total Stats High Mass" << setw(10) << total_highM << " cut = " << hHighM->GetXaxis()->GetBinLabel(start_bin) << endl;
  for(Int_t i=start_bin+10; i<hHighM->GetXaxis()->GetNbins()+1; i=i+10){
    Double_t percent = hHighM->GetBinContent(i);
    percent = percent/total_highM;
    percent = 100*percent;
    cout << setw(25) << hHighM->GetXaxis()->GetBinLabel(i) << setw(10) << hHighM->GetBinContent(i) << " " << percent << endl;
  }
}
