void DrawSeparate(TCanvas *c1, TH1D **h1, TString *title, Int_t nhist){
  for (Int_t i=0; i<nhist; i++) {
    c1->cd(i+1);
    h1[i]->Draw();
    
    h1[i]->SetTitle(title[i].Data());
  }

}

void DrawSame(TH1D **h1, Int_t nhist, Bool_t norm=false){
  for (Int_t i=0, icolor=1; i<nhist; i++, icolor++) {
    //Scale hist
    if(norm) h1[i]->Scale(1.0/h1[i]->GetEntries() );

    //Set colors
    if (icolor==5) icolor++;//no Yellow
    else if(icolor==10) icolor=1;
    h1[i]->SetLineColor(icolor);
    
    if (!i) h1[i]->Draw();
    else {
      h1[i]->Draw("same");
    }
  }
  
}

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
  TH1D *hCut_MuPTheta[nCuts],*hCut_MuPPhi[nCuts],*hCut_MuPqP[nCuts];
  TH1D *hCut_MuMTheta[nCuts],*hCut_MuMPhi[nCuts],*hCut_MuMqP[nCuts];
  TH1D *hCut_xN[nCuts], *hCut_xPi[nCuts], *hCut_xF[nCuts];
  TH1D *hCut_qT[nCuts];
  TH1D *hCut_PhiPhoton[nCuts], *hCut_PhiPhoton_gen[nCuts];
  TH1D *hCut_PhiS_simple[nCuts], *hCut_PhiS_simple_gen[nCuts];
  TH1D *hCut_PhiPhoton_clone[nCuts], *hCut_PhiPhoton_gen_clone[nCuts];
  TH1D *hCut_PhiS_simple_clone[nCuts], *hCut_PhiS_simple_gen_clone[nCuts];
  for (Int_t i=0; i<nCuts; i++) {
    hCut_VxZ[i]=
      (TH1D*)f_input->Get(Form("VxZ_CutImpact/%s",cutNames[i].Data()));
    
    hCut_MuPTheta[i]=
      (TH1D*)f_input->Get(Form("MuPTheta_CutImpact/%s",cutNames[i].Data()));
    hCut_MuPPhi[i]=
      (TH1D*)f_input->Get(Form("MuPPhi_CutImpact/%s",cutNames[i].Data()));
    hCut_MuPqP[i]=
      (TH1D*)f_input->Get(Form("MuPqP_CutImpact/%s",cutNames[i].Data()));

    hCut_MuMTheta[i]=
      (TH1D*)f_input->Get(Form("MuMTheta_CutImpact/%s",cutNames[i].Data()));
    hCut_MuMPhi[i]=
      (TH1D*)f_input->Get(Form("MuMPhi_CutImpact/%s",cutNames[i].Data()));
    hCut_MuMqP[i]=
      (TH1D*)f_input->Get(Form("MuMqP_CutImpact/%s",cutNames[i].Data()));

    hCut_xN[i]=
      (TH1D*)f_input->Get(Form("xN_CutImpact/%s",cutNames[i].Data()));
    hCut_xPi[i]=
      (TH1D*)f_input->Get(Form("xPi_CutImpact/%s",cutNames[i].Data()));
    hCut_xF[i]=
      (TH1D*)f_input->Get(Form("xF_CutImpact/%s",cutNames[i].Data()));
    hCut_qT[i]=
      (TH1D*)f_input->Get(Form("qT_CutImpact/%s",cutNames[i].Data()));

    hCut_PhiPhoton[i]=
      (TH1D*)f_input->Get(Form("PhiPhoton_CutImpact/%s",cutNames[i].Data()));
    hCut_PhiPhoton_gen[i]=
      (TH1D*)f_input->Get(Form("PhiPhoton_gen_CutImpact/%s",cutNames[i].Data()));
    hCut_PhiS_simple[i]=
      (TH1D*)f_input->Get(Form("PhiS_simple_CutImpact/%s",cutNames[i].Data()));
    hCut_PhiS_simple_gen[i]=
      (TH1D*)f_input->Get(Form("PhiS_simple_gen_CutImpact/%s",cutNames[i].Data()));

    //Clones
    hCut_PhiPhoton_clone[i] = (TH1D*)hCut_PhiPhoton[i]->Clone();
    hCut_PhiPhoton_gen_clone[i] = (TH1D*)hCut_PhiPhoton_gen[i]->Clone();
    hCut_PhiS_simple_clone[i] = (TH1D*)hCut_PhiS_simple[i]->Clone();
    hCut_PhiS_simple_gen_clone[i] = (TH1D*)hCut_PhiS_simple_gen[i]->Clone();
  }

  TCanvas* cVxZ = new TCanvas();
  cVxZ->Divide(2, 4);
  DrawSeparate(cVxZ, hCut_VxZ, cutNames, nCuts);
  TCanvas* cVxZ_same = new TCanvas();
  DrawSame(hCut_VxZ, nCuts);

  TCanvas* cPhiPhoton = new TCanvas();//Phi Photon
  cPhiPhoton->Divide(2, 4);
  DrawSeparate(cPhiPhoton, hCut_PhiPhoton, cutNames, nCuts);
  TCanvas* cPhiPhoton_same = new TCanvas();
  DrawSame(hCut_PhiPhoton, nCuts);
  TCanvas* cPhiPhoton_same_norm = new TCanvas();
  DrawSame(hCut_PhiPhoton_clone, nCuts, true);

  TCanvas* cPhiPhoton_gen = new TCanvas();//Phi Photon generated
  cPhiPhoton_gen->Divide(2, 4);
  DrawSeparate(cPhiPhoton_gen, hCut_PhiPhoton_gen, cutNames, nCuts);
  TCanvas* cPhiPhoton_gen_same = new TCanvas();
  DrawSame(hCut_PhiPhoton_gen, nCuts);
  TCanvas* cPhiPhoton_gen_same_norm = new TCanvas();
  DrawSame(hCut_PhiPhoton_gen_clone, nCuts, true);

  /*TCanvas* cPhiS_simple = new TCanvas();//PhiS simple
  cPhiS_simple->Divide(2, 4);
  DrawSeparate(cPhiS_simple, hCut_PhiS_simple, cutNames, nCuts);
  TCanvas* cPhiS_simple_same = new TCanvas();
  DrawSame(hCut_PhiS_simple, nCuts);
  TCanvas* cPhiS_simple_same_norm = new TCanvas();
  DrawSame(hCut_PhiS_simple_clone, nCuts, true);

  TCanvas* cPhiS_simple_gen = new TCanvas();//PhiS simple generated
  cPhiS_simple_gen->Divide(2, 4);
  DrawSeparate(cPhiS_simple_gen, hCut_PhiS_simple_gen, cutNames, nCuts);
  TCanvas* cPhiS_simple_gen_same = new TCanvas();
  DrawSame(hCut_PhiS_simple_gen, nCuts);
  TCanvas* cPhiS_simple_gen_same_norm = new TCanvas();
  DrawSame(hCut_PhiS_simple_gen_clone, nCuts, true);*/
}
