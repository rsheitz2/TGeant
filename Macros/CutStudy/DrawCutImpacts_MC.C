void DrawCutImpacts_MC(TString fname=""){

  if(fname == ""){
    cout << " " << endl;
    
    cout << "Macro draws cut impacts on one canvas" << endl;
    cout << " " << endl;
    cout << "Usage:" << endl;
    cout << "root \'DrawCutImpact.C(\"inputfile.root\")\'" << endl;
    
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  TFile *f_input = TFile::Open(fname);
  if (!f_input->IsOpen()) {
    cout << fname << " did not open" << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nCuts = 7;//Number of impact cut hist
  TString cutNames[nCuts] = {"AllData", "xPion,xN,xF", "0.4<qT<5",
			       "TargetZ-cut", "TargetRadius", "Positive_Pz_In",
			       "Physical_Pxyz_tr1_tr2"};
  TH1D *hCut_VxZ[nCuts];
  for (Int_t i=0; i<nCuts; i++) {
    hCut_VxZ[i]=(TH1D*)f_input->Get(Form("VxZ_CutImpact/%s",cutNames[i].Data()));
  }

  TCanvas* c1 = new TCanvas();
  c1->Divide(2, 4);
  for (Int_t i=0; i<nCuts; i++) {
    c1->cd(i+1);
    hCut_VxZ[i]->Draw();
  }

  TCanvas* c2 = new TCanvas();
  for (Int_t i=0; i<nCuts; i++) {
    if (!i) hCut_VxZ[i]->Draw();
    else {
      hCut_VxZ[i]->Draw("same");
      hCut_VxZ[i]->SetLineColor(i+2);
    }
  }


  
}
