#ifndef ENDJOBS_H
#define ENDJOBS_H

inline void NameCutImpacts(const Phast &ph, TString var, TString *cutNames,
		    Int_t nCuts){

  ph.h_file->cd();
  TDirectory* dir = (TDirectory*)gROOT->FindObject(Form("%s_CutImpact",
							var.Data()) );
  if(dir != NULL) dir->cd();
  for (Int_t i=0; i<nCuts; i++) {
    TH1D* h1 = (TH1D*)gROOT->FindObject(Form("hCut%i_%s", i, var.Data() ) );
    if(h1 != NULL) h1->SetName(cutNames[i] );  
  }
}//NameCutImpacts


inline void Name2D_CutImpacts(const Phast &ph, TString var, TString *cutNames,
		    Int_t nCuts){

  ph.h_file->cd();
  TDirectory* dir = (TDirectory*)gROOT->FindObject(Form("%s_CutImpact",
							var.Data()) );
  if(dir != NULL) dir->cd();
  for (Int_t i=0; i<nCuts; i++) {
    TH2D* h2 = (TH2D*)gROOT->FindObject(Form("hCut%i_%s", i, var.Data() ) );
    if(h2 != NULL) h2->SetName(cutNames[i] );  
  }
}//NameCutImpacts

#endif
