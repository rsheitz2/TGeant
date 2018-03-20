#ifndef common_h
#define common_h
  
  // C++ 
  #include <iostream>
  #include <fstream>
  #include <cmath>
  #include <cstdio>
  #include <cstring>
  #include <ctime>
  #include <iomanip>
  #include <vector>
  #include <map>
  #include <cstdlib>
  #include <algorithm>
  
  // ROOT
  #include <TROOT.h>
  #include <TChain.h>
  #include <TFile.h>
  #include <TSelector.h>
  #include <TH1.h>
  #include <TH2.h>
  #include <TProfile.h>
  #include <TNtuple.h>
  #include <TRandom.h>
  #include <TTree.h>
  #include <TStyle.h>
  #include <TGraph.h>
  #include <TCanvas.h>
  #include <TGraphErrors.h>
  #include <TMath.h>
  #include <TLorentzVector.h>
  #include <TVector3.h>
  #include <TLorentzRotation.h>
  
  // Phast 
  #include "Phast.h"
  #include "PaSetup.h"
  #include "PaEvent.h"
  #include "PaParticle.h"
  #include "PaAlgo.h"
  #include "PaMCtrack.h"
  #include "PaMCvertex.h"
  #include "PaMCgen.h"
  #include "G3part.h"
  #include "PaTrack.h"
  #include "PaMetaDB.h"
  
  
  //const Int_t    ALL     = 1 ;
  //const Double_t PI      = TMath::Pi();
  //const Double_t RAD2DEG = 180./PI;
  //const Double_t DEG2RAD = PI/180.;
  
  
  //////////
  // Mass                                             
  //////////
  
  // G3partMass[] -> PDG mass (modified; 2013.07.03.)
  const Double_t MASS_MU    = 0.1056583715 ; // muon
  const Double_t MASS_PI    = 0.13957018   ; // pion
  const Double_t MASS_P     = 0.938272046  ; // proton
  const Double_t MASS_N     = 0.939565379  ; // neutron
  const Double_t MASS_K     = 0.493677  ;
  const Double_t MASS_L1116 = 1.115683 ;
  
  //Int_t PointsHodoscopes2    ( const PaTrack& tr, Int_t mode );
  //Int_t PointsHodoscopes2opp ( const PaTrack& tr, Int_t mode );

  Int_t PointsHodoscopes3(   const PaTrack& tr, Bool_t* flagLAS,  Bool_t* flagOUTER,  Bool_t* flagMIDDLE );
  Int_t PointsHodoscopes3opp(const PaTrack& tr, Bool_t* flagLAS,  Bool_t* flagOUTER,  Bool_t* flagMIDDLE );
  Int_t PointsHodoscopes3oppV(const PaTPar& parV, Bool_t* flagLAS,  Bool_t* flagOUTER,  Bool_t* flagMIDDLE );
  
  Bool_t dataSkimLIP(PaEvent& e);
  
#endif // #ifndef common_h


