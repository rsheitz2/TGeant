#include "common.hxx"
#include "dy_functions.h"
#include "dy_variables.h"

//////////////Written by Robert Heitz///////////
////User event: 
////-Gets information for dimuon events
////-To be used with Real Data mDSTs
///Specify -U10 for real data and -U20 for MC data
////-Doesn't quiet get all the information for MC data (but has most of it)
//// Can now track down cut impacts
////
////Testing to determine the correct cuts
///////////////////////////////////////////////

//const static Int_t maxCuts = 15;//Should be the same as nCuts in UserJobEnd420.cc
//const static Int_t nCutHist = 12;
//const static Int_t n2D_cutHist = 3;

void UserEvent420(PaEvent& e){
  //Tree
  static TTree* Particles(NULL);

  //Variables for dimuom tree
  //////////////////////////////////////
  //Setup globally in dy_variables.h
  static Bool_t isMonteCarlo;
  static Double_t SM1_p1x, SM1_p1y, SM1_p2x, SM1_p2y;
  static Double_t SM2_p1x, SM2_p1y, SM2_p2x, SM2_p2y;
  static Double_t HG01_p1x, HG01_p1y, HG02_y1_p1x, HG02_y1_p1y;
  static Double_t HG02_y2_p1x, HG02_y2_p1y;
  static Double_t HG01_p2x, HG01_p2y, HG02_y1_p2x, HG02_y1_p2y;
  static Double_t HG02_y2_p2x, HG02_y2_p2y;
 
  //Histograms
  static TH1D* h1[50];
  //static TH2D* h2[10];
  static TH1D *hCut_VxZ[maxCuts];
  static TH1D*hCut_MuPTheta[maxCuts],*hCut_MuPPhi[maxCuts],*hCut_MuPqP[maxCuts];
  static TH1D*hCut_MuMTheta[maxCuts],*hCut_MuMPhi[maxCuts],*hCut_MuMqP[maxCuts];
  static TH1D *hCut_xN[maxCuts], *hCut_xPi[maxCuts], *hCut_xF[maxCuts];
  static TH1D *hCut_qT[maxCuts], *hCut_PhiPhoton[maxCuts];
  static TH1D *hCut_PhiS[maxCuts], *hCut_PhiS_simple[maxCuts];
  
  static TH2D *hCut_MuP_PxPy[maxCuts], *hCut_MuM_PxPy[maxCuts];
  static TH2D *hCut_Beam_PxPy[maxCuts];
  
  static TH1D *hCuts[nCutHist+2][maxCuts];
  static TH2D *h2DCuts[n2D_cutHist][maxCuts];

  static bool first(true);
  if (first){
    Phast::Ref().HistFileDir("UserEvent420");

    if (Phast::Ref().UserFlag(0) == 20) isMonteCarlo = true;
    else if (Phast::Ref().UserFlag(0) == 10) isMonteCarlo = false;
    else {
      cout << " " << endl;
      cout << "Please specify -U20 for MC or -U real data" << endl;
      cout << " " << endl;
      exit(EXIT_FAILURE);
    }
    
    //////Making Trees////
    Particles = new TTree("Particles", "Particles");
    
    if (isMonteCarlo) SetupMCTreeVariables(Particles);
    else SetupRealTreeVariables(Particles);

    Particles->Branch("SM1_p1x", &SM1_p1x, "SM1_p1x/D");
    Particles->Branch("SM1_p1y", &SM1_p1y, "SM1_p1y/D");
    Particles->Branch("SM1_p2x", &SM1_p2x, "SM1_p2x/D");
    Particles->Branch("SM1_p2y", &SM1_p2y, "SM1_p2y/D");
    Particles->Branch("SM2_p1x", &SM2_p1x, "SM2_p1x/D");
    Particles->Branch("SM2_p1y", &SM2_p1y, "SM2_p1y/D");
    Particles->Branch("SM2_p2x", &SM2_p2x, "SM2_p2x/D");
    Particles->Branch("SM2_p2y", &SM2_p2y, "SM2_p2y/D");
    Particles->Branch("HG01_p1x", &HG01_p1x, "HG01_p1x/D");
    Particles->Branch("HG01_p1y", &HG01_p1y, "HG01_p1y/D");
    Particles->Branch("HG01_p2x", &HG01_p2x, "HG01_p2x/D");
    Particles->Branch("HG01_p2y", &HG01_p2y, "HG01_p2y/D");
    Particles->Branch("HG02_y1_p1x", &HG02_y1_p1x, "HG02_y1_p1x/D");
    Particles->Branch("HG02_y1_p1y", &HG02_y1_p1y, "HG02_y1_p1y/D");
    Particles->Branch("HG02_y1_p2x", &HG02_y1_p2x, "HG02_y1_p2x/D");
    Particles->Branch("HG02_y1_p2y", &HG02_y1_p2y, "HG02_y1_p2y/D");
    Particles->Branch("HG02_y2_p1x", &HG02_y2_p1x, "HG02_y2_p1x/D");
    Particles->Branch("HG02_y2_p1y", &HG02_y2_p1y, "HG02_y2_p1y/D");
    Particles->Branch("HG02_y2_p2x", &HG02_y2_p2x, "HG02_y2_p2x/D");
    Particles->Branch("HG02_y2_p2y", &HG02_y2_p2y, "HG02_y2_p2y/D");
    
    //Statisical changes from cuts
    h1[21] = new TH1D("DiMuonCuts", "Di-Muon cuts", 130, 0, 130);
        
    //Tracking cut impacts
    Int_t ih = 0;
    HistArraySetup(hCut_VxZ, hCuts, 500, -500, 100, ih, "VxZ"); ih++;
    HistArraySetup(hCut_MuPTheta, hCuts, 100, 0, 0.3, ih, "MuPTheta"); ih++;
    HistArraySetup(hCut_MuPPhi, hCuts, 100, -TMath::Pi(), TMath::Pi(), ih,
		   "MuPPhi"); ih++;
    HistArraySetup(hCut_MuPqP, hCuts, 200, 0, 200, ih, "MuPqP"); ih++;
    HistArraySetup(hCut_MuMTheta, hCuts, 100, 0, 0.3, ih, "MuMTheta"); ih++;
    HistArraySetup(hCut_MuMPhi, hCuts, 100, -TMath::Pi(), TMath::Pi(), ih,
		   "MuMPhi"); ih++;
    HistArraySetup(hCut_MuMqP, hCuts, 200, -200, 0, ih, "MuMqP"); ih++;
    HistArraySetup(hCut_xN, hCuts, 100, 0, 1, ih, "xN"); ih++;
    HistArraySetup(hCut_xPi, hCuts, 100, 0, 1, ih, "xPi"); ih++;
    HistArraySetup(hCut_xF, hCuts, 200, -1, 1, ih, "xF"); ih++;
    HistArraySetup(hCut_qT, hCuts, 200, 0, 5, ih, "qT"); ih++;
    HistArraySetup(hCut_PhiPhoton, hCuts, 200, -1.0*TMath::Pi(), TMath::Pi(),
		   ih, "PhiPhoton"); ih++;
    HistArraySetup(hCut_PhiS, hCuts, 200, -1.0*TMath::Pi(), TMath::Pi(),
		   ih, "PhiS"); ih++;
    HistArraySetup(hCut_PhiS_simple, hCuts, 200, -1.0*TMath::Pi(), TMath::Pi(),
		   ih, "PhiS_simple"); ih++;

    ih = 0;
    Hist2D_ArraySetup(hCut_MuP_PxPy, h2DCuts, 100, -5, 5, 100, -5, 5, ih,
		      "MuP_PxPy"); ih++;
    Hist2D_ArraySetup(hCut_MuM_PxPy, h2DCuts, 100, -5, 5, 100, -5, 5, ih,
		      "MuM_PxPy"); ih++;
    Hist2D_ArraySetup(hCut_Beam_PxPy, h2DCuts, 100, -2, 2, 100, -2, 2, ih,
		      "Beam_PxPy"); ih++;
    
    first = false;
  }

  for (Int_t i_p1=0; i_p1<e.NParticle()-1; i_p1++){//p1 loop
    const PaParticle& p1_tmp = e.vParticle(i_p1);
    if (p1_tmp.iTrack() < 0 ) continue;

    for (Int_t i_p2=i_p1+1; i_p2<e.NParticle(); i_p2++){//p2 loop
      const PaParticle& p2_tmp = e.vParticle(i_p2);
      if (p2_tmp.iTrack() < 0 ) continue;

      //All two track combinations 
      Int_t cut_bin = 1, cut_space = 10, icut = 0;
      Int_t i_mHist = 0;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space; icut++;
      
      Int_t vx_p1[50]; //p1 vertices
      for (Int_t i_vx=0; i_vx<p1_tmp.NVertex(); i_vx++) {
	vx_p1[i_vx] = p1_tmp.iVertex(i_vx);
	const PaVertex& tmx_vx = e.vVertex(vx_p1[i_vx]);
	if (tmx_vx.IsPrimary() ) numVxpri_p1++;
      }
      Int_t vx_p2[50]; //p2 vertices
      for (Int_t i_vx=0; i_vx<p2_tmp.NVertex(); i_vx++) {
	vx_p2[i_vx] = p2_tmp.iVertex(i_vx);
	const PaVertex& tmx_vx = e.vVertex(vx_p2[i_vx]);
	if (tmx_vx.IsPrimary() ) numVxpri_p2++;
      }
      vector<Int_t> common_vx;
      common2parVx(common_vx, vx_p1, p1_tmp.NVertex(), vx_p2,p2_tmp.NVertex() );

      if (common_vx.size() == 0) continue; //Has common vertex
      for (vector<Int_t>::iterator ic=common_vx.begin(); ic!=common_vx.end();
	   ic++) {
	if(e.vVertex(*ic).IsPrimary() ) h1[21]->Fill(cut_bin-1 ); 
      }
      cut_bin += cut_space; icut++;
      
      //////////Choose Best Primary vertex
      Int_t best_iVx = -1;
      //Best primary vertex not made for DY
      vector<Int_t>::iterator iv = find(common_vx.begin(), common_vx.end(),
					e.iBestPrimaryVertex() );
      if (iv!=common_vx.end() ){
	best_iVx = *iv;
      }
      else {//Does not have BPV tagged
	for (UInt_t i_vx=0; i_vx<common_vx.size(); i_vx++){//vertex loop
	  const PaVertex& new_vx = e.vVertex(common_vx.at(i_vx) );
	  if (new_vx.IsPrimary() ){//Primary vx
	    if (best_iVx == -1) {
	      best_iVx = common_vx.at(i_vx);
	      continue;
	    }
	    const PaVertex& best_vx = e.vVertex(best_iVx);
	    Double_t best_chi2_red = best_vx.Chi2();
	    Double_t new_chi2_red = new_vx.Chi2();

	    if (new_chi2_red < best_chi2_red) best_iVx = common_vx.at(i_vx);
	  }//Primary vx
	}//vertex loop
      }//Does not have BPV taged
	      
      if (best_iVx == -1) continue;//Vertex is best primary
      const PaVertex& vertex = e.vVertex(best_iVx);

      //Best primary vertex
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      hCut_VxZ[icut]->Fill(vertex.Z() ); icut++;

      if (p1_tmp.Q() == p2_tmp.Q() ) continue;//Opposite signs
      //p1 for positive muon, p2 or negative muon
      const PaParticle &p1 = (p1_tmp.Q() > 0) ?
	e.vParticle(i_p1) : e.vParticle(i_p2);
      const PaParticle &p2 = (p2_tmp.Q() < 0) ?
	e.vParticle(i_p2) : e.vParticle(i_p1);
      const PaTrack& tr1 = e.vTrack(p1.iTrack() );
      const PaTrack& tr2 = e.vTrack(p2.iTrack() );
      const PaTPar& traj_p1 = p1.ParInVtx(best_iVx);
      const PaTPar& traj_p2 = p2.ParInVtx(best_iVx);
      const PaParticle& pIn = e.vParticle(vertex.InParticle() );
      const PaTPar& traj_pIn = pIn.ParInVtx(best_iVx);
      const PaTrack& trIn = e.vTrack(pIn.iTrack() );
      TLorentzVector lv_p1_Mu = traj_p1.LzVec(M_mu);
      TLorentzVector lv_p2_Mu = traj_p2.LzVec(M_mu);
      TLorentzVector lv_diMu = lv_p1_Mu + lv_p2_Mu;
      TLorentzVector lv_pIn = traj_pIn.LzVec(M_pi);
      
      x_beam = lv_diMu.Mag2()/(2*lv_diMu.Dot(lv_pIn) );
      TLorentzVector lv_target(0, 0, 0, M_proton);
      x_target = lv_diMu.Mag2()/(2*lv_diMu.Dot(lv_target) );
      x_feynman = x_beam - x_target;
      q_transverse = lv_diMu.Vect().Cross(lv_pIn.Vect() ).Mag()
	/(lv_pIn.Vect().Mag() );
      
      Double_t SpinDir;
      if (!isMonteCarlo && PaAlgo::GetDYtargetPolarization(e, vertex.Z() > 0) ){
	SpinDir = -1.0;}
      else SpinDir = 1.0;
      TLorentzVector lv_Spin(0, SpinDir, 0, 0);
      TLorentzVector lv_Spin_simple(0, 1.0, 0, 0);
      Bool_t inNH3 = ( (vertex.Z() > -294.5 && vertex.Z() < -239.3) ||
		       (vertex.Z() > -219.5 && vertex.Z() < -164.3) )
	? true : false; 

      //Target frame:
      TLorentzVector lv_pIn_TF(lv_pIn);
      TLorentzVector lv_target_TF(lv_target);
      TLorentzVector lv_Spin_TF(lv_Spin);
      TLorentzVector lv_Spin_simple_TF(lv_Spin_simple);
      TLorentzVector lv_p1_Mu_TF(lv_p1_Mu);
      TLorentzVector lv_p2_Mu_TF(lv_p2_Mu);
      TLorentzVector lv_virtualPhoton_TF(lv_diMu);
      if (inNH3){
	align_wrt_beam_photon(lv_pIn_TF, lv_target_TF, lv_Spin_TF,
			      lv_Spin_simple_TF, lv_virtualPhoton_TF,
			      lv_p1_Mu_TF, lv_p2_Mu_TF);
      }

      Double_t cut_variables[nCutHist] = {vertex.Z(), traj_p1.Theta(),
					  traj_p1.Phi(), traj_p1.qP(),
					  traj_p2.Theta(), traj_p2.Phi(),
					  traj_p2.qP(), x_beam, x_target,
					  x_feynman, q_transverse,
					  lv_diMu.Phi()};
      
      Double_t cut2D_variables[2*n2D_cutHist] = {lv_p1_Mu.X(), lv_p1_Mu.Y(),
						 lv_p2_Mu.X(), lv_p2_Mu.Y(),
						 lv_pIn.X(), lv_pIn.Y()};

      //Opposite Q
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      } 
      FillCuts(hCuts, cut_variables, icut); icut++;


      //Dimuon trigger fired
      if ( !((e.TrigMask() >> 2) & 1) && !((e.TrigMask() >> 8) & 1) ) continue;
      if ( (e.TrigMask() & 1) ) continue;//veto middle
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
                  
      //XX0 > 30
      if (tr1.XX0() < 30.0 || tr2.XX0() < 30.0) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      if (lv_diMu.M() <= 1.0) i_mHist = 30;
      else if (lv_diMu.M() > 1.0 && lv_diMu.M() <= 2.5) i_mHist = 31;
      else if (lv_diMu.M() > 2.5 && lv_diMu.M() <= 4.3) i_mHist = 32;
      else if (lv_diMu.M() > 4.3 && lv_diMu.M() <= 8.5) i_mHist = 33;
      else i_mHist = 34;

      //if (i_mHist != 33) continue;//Mass [4.3, 8.5]
      if (i_mHist != 32) continue;//Mass [4.3, 8.5]
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
            
      //Zfirst/Zlast
      if (tr1.ZFirst() > 300) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      if (tr1.ZLast() < 1500) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut); 
      FillCuts(hCuts, cut_variables, icut); icut++;
            
      if (tr2.ZFirst() > 300) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      if (tr2.ZLast() < 1500) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
            
      /*if (tr1.ZFirst() > 300 || tr1.ZLast() < 1500) continue;
	h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
	hCut_VxZ[icut]->Fill(vertex.Z() ); icut++;
	if (tr2.ZFirst() > 300 || tr2.ZLast() < 1500) continue;
	h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
	Fill2D_Cuts(h2DCuts, cut2D_variables, icut); 
	hCut_VxZ[icut]->Fill(vertex.Z() ); icut++;*/

      //Track time defined
      if (tr1.MeanTime() == 1e+10 || tr2.MeanTime() == 1e+10) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      //Track times with 5ns
      if (TMath::Abs(tr1.MeanTime() - tr2.MeanTime() ) > 5) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;

      //Tracks Chi/ndf < 10
      if (tr1.Chi2tot()/tr1.Ndf() > 10) continue;
      if (tr2.Chi2tot()/tr2.Ndf() > 10) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      Double_t track_chi2 = tr1.Chi2tot() + tr2.Chi2tot() + trIn.Chi2tot();
      track_chi2 = track_chi2/(1.0*(tr1.Ndf() + tr2.Ndf() + trIn.Ndf()) );

      //Beam decay cut
      /*if (TMath::Abs(traj_p1.qP() ) + TMath::Abs(traj_p2.qP() )
	  > 190.0) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      FillCuts(hCuts, cut_variables, icut); icut++;*/

      //Trigger validation
      Bool_t trigValidation = false;
      if ( ((e.TrigMask() >> 8) & 1)//Last_Last
	   && pointsToLAS(tr1) && pointsToLAS(tr2)) trigValidation = true;
      else if ( (e.TrigMask() >> 2) & 1){//Last_Outer
	if (pointsToOuter(tr1) ){
	  if (pointsToLAS(tr2)) trigValidation = true;
	}
        if (pointsToLAS(tr1) ){
	  if (pointsToOuter(tr2) ) trigValidation = true;
	}
      }
      
      if (!trigValidation) continue;
      h1[21]->Fill(cut_bin-1); cut_bin += cut_space;
      Fill2D_Cuts(h2DCuts, cut2D_variables, icut);
      if (inNH3) {
	hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
	hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      }
      FillCuts(hCuts, cut_variables, icut); icut++;
      
      //TLorentzVector lv_pIn = traj_pIn.LzVec(M_pi);

      ///////Numerating values for tree
      //////////////////////////
      AssignPhastParVariables(e, p1, p2, pIn);
      AssignPhastTrackVariables(vertex, tr1, tr2, trIn);
      AssignPhastTrajVariables(traj_p1, traj_p2, traj_pIn);
      
      //numVxpri_p1 determined in vertex loop
      //Positively charged outgoing muon trajectory parameters at vertex
      vP1_X = lv_p1_Mu.X();
      vP1_Y = lv_p1_Mu.Y();
      vP1_Z = lv_p1_Mu.Z();
      vP1_E = lv_p1_Mu.E();
      
      //numVxpri_p2 determined in vertex loop
      //Negatively charged outgoing muon trajectory parameters at vertex
      vP2_X = lv_p2_Mu.X();
      vP2_Y = lv_p2_Mu.Y();
      vP2_Z = lv_p2_Mu.Z();
      vP2_E = lv_p2_Mu.E();

      //Virtual photon specific/Dimuon
      vPhoton_X = lv_diMu.X();
      vPhoton_Y = lv_diMu.Y();
      vPhoton_Z = lv_diMu.Z();
      vPhoton_E = lv_diMu.E();
      vDiMuon_invM = lv_diMu.M();//GeV
      vOpenAngle = lv_p1_Mu.Vect().Angle(lv_p2_Mu.Vect() );
      
      //Beam pion trajectory parameters at vertex
      beam_X = lv_pIn.X();
      beam_Y = lv_pIn.Y();
      beam_Z = lv_pIn.Z();
      beam_E = lv_pIn.E();

      if (isMonteCarlo) {//MC Data
	if (tr1.iMCtrack() >= 0 && tr2.iMCtrack() >= 0 && trIn.iMCtrack() >= 0){

	  const PaMCtrack& MCtr1 = e.vMCtrack(tr1.iMCtrack() );
	  const PaMCtrack& MCtr2 = e.vMCtrack(tr2.iMCtrack() );
	  const PaMCtrack& MCtrIn = e.vMCtrack(trIn.iMCtrack() );
	  TLorentzVector lv_MCtr1 = MCtr1.LzVec();
	  TLorentzVector lv_MCtr2 = MCtr2.LzVec();
	  TLorentzVector lv_MCtrIn( -MCtrIn.P(0), -MCtrIn.P(1), -MCtrIn.P(2),
				    MCtrIn.E() );//generated info saves neg mom

	  TLorentzVector lv_MCdiMu = lv_MCtr1 + lv_MCtr2;
	  MC_x_beam = lv_MCdiMu.Mag2()/(2*lv_MCdiMu.Dot(lv_MCtrIn) );
	  MC_x_target = lv_MCdiMu.Mag2()/(2*lv_MCdiMu.Dot(lv_target) );
	  MC_x_feynman = MC_x_beam - MC_x_target;
	  MC_q_transverse = lv_MCdiMu.Vect().Cross(lv_MCtrIn.Vect() ).Mag()
						/(lv_MCtrIn.Vect().Mag() );

	  AssignMCtrack(MCtr1, MCtr2, MCtrIn);

	  theta_MCtr1 = lv_MCtr1.Theta();
	  phi_MCtr1 = lv_MCtr1.Phi();
	  vMCtr1_X = lv_MCtr1.X();
	  vMCtr1_Y = lv_MCtr1.Y();
	  vMCtr1_Z = lv_MCtr1.Z();
	  vMCtr1_E = lv_MCtr1.E();

	  theta_MCtr2 = lv_MCtr2.Theta();
	  phi_MCtr2 = lv_MCtr2.Phi();
	  vMCtr2_X = lv_MCtr2.X();
	  vMCtr2_Y = lv_MCtr2.Y();
	  vMCtr2_Z = lv_MCtr2.Z();
	  vMCtr2_E = lv_MCtr2.E();

	  theta_MCtrIn = lv_MCtrIn.Theta();
	  phi_MCtrIn = lv_MCtrIn.Phi();
	  vMCtrIn_X = lv_MCtrIn.X();
	  vMCtrIn_Y = lv_MCtrIn.Y();
	  vMCtrIn_Z = lv_MCtrIn.Z();
	  vMCtrIn_E = lv_MCtrIn.E();

	}
	else AssignMCtrackVoid();
      }
      else AssignTarget(e);//Real Data

      PaTPar p1_sm1, p1_sm2;
      Bool_t SM1_p1_reach = tr1.Extrapolate(363.7, p1_sm1);
      Bool_t SM2_p1_reach = tr1.Extrapolate(1825.3, p1_sm2);
      PaTPar p2_sm1, p2_sm2;
      Bool_t SM1_p2_reach = tr2.Extrapolate(363.7, p2_sm1);
      Bool_t SM2_p2_reach = tr2.Extrapolate(1825.3, p2_sm2);

      if (SM1_p1_reach && SM1_p2_reach){
	//h2[0]->Fill(p1_sm1.X(), p1_sm1.Y() );
	//h2[0]->Fill(p2_sm1.X(), p2_sm1.Y() );

	SM1_p1x = p1_sm1.X();
	SM1_p1y = p1_sm1.Y();
	SM1_p2x = p2_sm1.X();
	SM1_p2y = p2_sm1.Y();
      }
      else{
	SM1_p1x = -999.9;
	SM1_p1y = -999.9;
	SM1_p2x = -999.9;
	SM1_p2y = -999.9;
      }
      if (SM2_p1_reach && SM2_p2_reach){
	//h2[1]->Fill(p1_sm2.X(), p1_sm2.Y() );
	//h2[1]->Fill(p2_sm2.X(), p2_sm2.Y() );

	SM2_p1x = p1_sm1.X();
	SM2_p1y = p1_sm1.Y();
	SM2_p2x = p2_sm1.X();
	SM2_p2y = p2_sm1.Y();
      }
      else{
	SM2_p1x = -999.9;
	SM2_p1y = -999.9;
	SM2_p2x = -999.9;
	SM2_p2y = -999.9;
      }

      PaTPar p1_HG1, p1_HG2_y1, p1_HG2_y2;
      PaTPar p2_HG1, p2_HG2_y1, p2_HG2_y2;
      Bool_t HG01_p1_reach = tr1.Extrapolate(583.0, p1_HG1);
      Bool_t HG02_y1_p1_reach = tr1.Extrapolate(1602.0, p1_HG2_y1);
      Bool_t HG02_y2_p1_reach = tr1.Extrapolate(1613.1001, p1_HG2_y2);
      Bool_t HG01_p2_reach = tr1.Extrapolate(583.0, p2_HG1);
      Bool_t HG02_y1_p2_reach = tr1.Extrapolate(1602.0, p2_HG2_y1);
      Bool_t HG02_y2_p2_reach = tr1.Extrapolate(1613.1001, p2_HG2_y2);
      if (HG01_p1_reach && HG01_p2_reach){
	HG01_p1x = p1_sm1.X();
	HG01_p1y = p1_sm1.Y();
	HG01_p2x = p2_sm1.X();
	HG01_p2y = p2_sm1.Y();
      }
      else{
	HG01_p1x = -999.9;
	HG01_p1y = -999.9;
	HG01_p2x = -999.9;
	HG01_p2y = -999.9;
      }
      if (HG02_y1_p1_reach && HG02_y1_p2_reach){
	HG02_y1_p1x = p1_sm1.X();
	HG02_y1_p1y = p1_sm1.Y();
	HG02_y1_p2x = p2_sm1.X();
	HG02_y1_p2y = p2_sm1.Y();
      }
      else{
	HG02_y1_p1x = -999.9;
	HG02_y1_p1y = -999.9;
	HG02_y1_p2x = -999.9;
	HG02_y1_p2y = -999.9;
      }
      if (HG02_y2_p1_reach && HG02_y2_p2_reach){
	HG02_y2_p1x = p1_sm1.X();
	HG02_y2_p1y = p1_sm1.Y();
	HG02_y2_p2x = p2_sm1.X();
	HG02_y2_p2y = p2_sm1.Y();
      }
      else{
	HG02_y2_p1x = -999.9;
	HG02_y2_p1y = -999.9;
	HG02_y2_p2x = -999.9;
	HG02_y2_p2y = -999.9;
      }


      Particles->Fill();
    }//p2 loop
  }//p1 loop

}//End UserEvert function
