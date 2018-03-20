void DrawSeparate(TCanvas *c1, TH2D **h2, TString *title, Int_t nhist){
  for (Int_t i=0; i<nhist; i++) {
    c1->cd(i+1);
    h2[i]->Draw("colz");
    
    h2[i]->SetTitle(title[i].Data());
  }

}


void DrawCutImpacts_2D_MC(TString fname=""){

  if(fname == ""){
    cout << " " << endl;
    
    cout << "Macro draws cut impacts for 2D histograms" << endl;
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
  TH2D *hCut_MuP_PxPy[nCuts], *hCut_MuM_PxPy[nCuts];
  TH2D *hCut_Beam_PxPy[nCuts];
  for (Int_t i=0; i<nCuts; i++) {

    hCut_MuP_PxPy[i]=
      (TH2D*)f_input->Get(Form("MuP_PxPy_CutImpact/%s",cutNames[i].Data()));
    hCut_MuM_PxPy[i]=
      (TH2D*)f_input->Get(Form("MuM_PxPy_CutImpact/%s",cutNames[i].Data()));
    hCut_Beam_PxPy[i]=
      (TH2D*)f_input->Get(Form("Beam_PxPy_CutImpact/%s",cutNames[i].Data()));
  }

  TCanvas* cMuP_PxPy = new TCanvas();//MuP PxPy
  cMuP_PxPy->Divide(2, 4);
  DrawSeparate(cMuP_PxPy, hCut_MuP_PxPy, cutNames, nCuts);

  TCanvas* cMuM_PxPy = new TCanvas();//MuM PxPy
  cMuM_PxPy->Divide(2, 4);
  DrawSeparate(cMuM_PxPy, hCut_MuM_PxPy, cutNames, nCuts);

  TCanvas* cBeam_PxPy = new TCanvas();//Beam PxPy
  cBeam_PxPy->Divide(2, 4);
  DrawSeparate(cBeam_PxPy, hCut_Beam_PxPy, cutNames, nCuts);
 
}
