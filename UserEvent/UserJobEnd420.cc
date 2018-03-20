#include<iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "Phast.h"
#include "EndJobs.h"


//
// It's a place for some "end-of-job
// actions (like histogram fits, histogram
// normalizations etc)
//
// For every UserEventN, complementary 
// UserJobEndN is called (if exists)
//

void UserJobEnd420()
{
  const Phast& ph = Phast::Ref();

  //if(ph.print) cout<<"[ UserJobEnd420 has been called ]"<<endl;
  cout<<"[ UserJobEnd420 has been called ]"<<endl;
  
  const Int_t nCuts = 15;//Should be the same as maxCuts in dy_variables.h
  TString cutNames[nCuts] = {"All_Di-particles", "CommonVX",
			     "HasCommonPrimaryVertex","Opposite_Q",
			     //"Opposite_Q", "XX0>30", "HighMass","ZfirstZlast",
			     "DimuonTrig", "XX0>30", "HighMass",
			     "tr1Zfirst", "tr1Zlast", "tr2Zfirst", "tr2Zlast",
			     "TrackT_defined", "TrackT_within5",
			     //"TracksX/ndf<10", "BeamDecay", "Trigger Val"};
			     "TracksX/ndf<10", "Trigger Val"};

  ph.h_file->cd();
  TDirectory* dir = (TDirectory*)gROOT->FindObject("UserEvent420");
  if(dir != NULL) dir->cd(); 
  TH1D* h1 = (TH1D*)gROOT->FindObject("DiMuonCuts");

  if(h1 != NULL) {
    for (Int_t i=0, cut_bin=1; i<nCuts; i++, cut_bin+=10) {
      h1->GetXaxis()->SetBinLabel(cut_bin, cutNames[i] );  
    }
  }

  NameCutImpacts(ph, "VxZ", cutNames, nCuts);
  NameCutImpacts(ph, "MuPTheta", cutNames, nCuts);
  NameCutImpacts(ph, "MuPPhi", cutNames, nCuts);
  NameCutImpacts(ph, "MuPqP", cutNames, nCuts);
  NameCutImpacts(ph, "MuMTheta", cutNames, nCuts);
  NameCutImpacts(ph, "MuMPhi", cutNames, nCuts);
  NameCutImpacts(ph, "MuMqP", cutNames, nCuts);
  NameCutImpacts(ph, "xN", cutNames, nCuts);
  NameCutImpacts(ph, "xPi", cutNames, nCuts);
  NameCutImpacts(ph, "xF", cutNames, nCuts);
  NameCutImpacts(ph, "qT", cutNames, nCuts);
  NameCutImpacts(ph, "PhiPhoton", cutNames, nCuts);
  NameCutImpacts(ph, "PhiS", cutNames, nCuts);
  NameCutImpacts(ph, "PhiS_simple", cutNames, nCuts);

  Name2D_CutImpacts(ph, "MuP_PxPy", cutNames, nCuts);
  Name2D_CutImpacts(ph, "MuM_PxPy", cutNames, nCuts);
  Name2D_CutImpacts(ph, "Beam_PxPy", cutNames, nCuts);
  
}
