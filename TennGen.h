//// Author: Charles Philip Hughes <chughe26@vols.utk.edu>
// Update: 2017-07-20 14:40:27-0500
// Copyright (c) 2017; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// University of Tennessee ALICE collaboration
//
// Header File for TennGen Class
// BE SURE TO INCLUDE THE FOLLOWING LINES IN THE CODE THAT USES THIS:
//
//  class TennGen;
//
//  and include the TennGen.h and TennGen.cxx files 
//  in the same directory you run your macro in


#ifndef ROOT_TennGen
#define ROOT_TennGen

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TennGen                                                              //
//                                                                      //
// Description of the public objects                                    //
//                                                                      // 
// void   SetCentralityBin(Int_t cb){CentralityBin = cb;}               //
// -- set centrality bin                                                //
//                                                                      //
// void InitializeBackground(); -- declare this to start background     //
//                                                                      //
//                                                                      //
// TClonesArray *GetBackground();  -- get background to your macro      //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"	
#include "TLatex.h"
#include "TRandom3.h"

class TennGen : public TObject {

private:
  Int_t CentralityBin;
  Int_t Seed;
  Int_t HarmonicsFlag;
  Double_t EtaRange;
  Bool_t PRINTQA;
  void MyPrivateFunction();
  TClonesArray *BkgdArray;
  //TF1 *pTdist;
  TF1 *pTdist_piPlus;
  TF1 *pTdist_piMinus;
  TF1 *pTdist_kPlus;
  TF1 *pTdist_kMinus;
  TF1 *pTdist_p;
  TF1 *pTdist_pbar;
  TF1 *pTdist_piZero;
  TF1 *v1_pi;
  TF1 *v2_pi;
  TF1 *v3_pi;
  TF1 *v4_pi;
  TF1 *v5_pi;
  TF1 *v1_K;
  TF1 *v2_K;
  TF1 *v3_K;
  TF1 *v4_K;
  TF1 *v5_K;
  TF1 *v1_P;
  TF1 *v2_P;
  TF1 *v3_P;
  TF1 *v4_P;
  TF1 *v5_P;
  //TRandom3 *pTdist;
  TRandom3 *etadist;
  TRandom3 *Harmonics_Phi_Dist_Rand;
  TRandom3 *Psi_1;
  TRandom3 *Psi_3;
  TRandom3 *Psi_5;
  TRandom3 *New_Rand;
  //TF1 *phidist;
  //TF1 *etadist;
  //TF1 *Rando;	
  void SetpiPlusParams();
  void SetpiMinusParams();
  void SetkPlusParams();
  void SetkMinusParams();
  void SetpParams();
  void SetpbarParams();
  void SetpiZeroParams();
  void Setv1_Params();
  void Setv2_Params();
  void Setv3_Params();
  void Setv4_Params();
  void Setv5_Params();
  static Double_t fitf(Double_t *x, Double_t *par);
  static Double_t BlastWavedNdptTimesPt(Double_t *x, Double_t *par);
  static Double_t MyIntegrandBG(Double_t *x, Double_t *par);
  static Double_t MyStaticBGdNdPtTimesPt(Double_t *x, Double_t *par);
 
  
  Double_t dNdPhi(Double_t phiPart, Double_t pT,
     Double_t Psi1 , Double_t Psi2, Double_t Psi3, Double_t Psi4, Double_t Psi5,
    Double_t v1,Double_t v2,Double_t v3,Double_t v4,Double_t v5);
  
  Double_t Psi_1_event;
  Double_t Psi_3_event;
  Double_t Psi_5_event;
  Double_t phi;
  Double_t v2_pi_eval;
  Double_t v3_pi_eval;
  Double_t v4_pi_eval;
  Double_t v5_pi_eval;
  Double_t v2_K_eval;
  Double_t v3_K_eval;
  Double_t v4_K_eval;
  Double_t v5_K_eval;
  Double_t v2_P_eval;
  Double_t v3_P_eval;
  Double_t v4_P_eval;
  Double_t v5_P_eval;  
  
  TFile *vn_QA_file = new TFile("Background_Generator_QA_file.root","RECREATE");
  //TDirectory *vn_QA_plots;
  //TDirectory *Event_Plane_QA_plots;
  TH2D *v1_pi_h = new TH2D("v1_pi_h" , "v_{1} vs p_{T} for #pi^{+}, #pi^{-}, #pi^{0}"  , 200 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v2_pi_h = new TH2D("v2_pi_h" , "v_{2} vs p_{T} for #pi^{+}, #pi^{-}, #pi^{0}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v3_pi_h = new TH2D("v3_pi_h" , "v_{3} vs p_{T} for #pi^{+}, #pi^{-}, #pi^{0}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v4_pi_h = new TH2D("v4_pi_h" , "v_{4} vs p_{T} for #pi^{+}, #pi^{-}, #pi^{0}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v5_pi_h = new TH2D("v5_pi_h" , "v_{5} vs p_{T} for #pi^{+}, #pi^{-}, #pi^{0}"  , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v1_K_h = new TH2D("v1_K_h" , "v_{1} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v2_K_h = new TH2D("v2_K_h" , "v_{2} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v3_K_h = new TH2D("v3_K_h" , "v_{3} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v4_K_h = new TH2D("v4_K_h" , "v_{4} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v5_K_h = new TH2D("v5_K_h" , "v_{5} vs p_{T} for K^{+}, K^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v1_P_h = new TH2D("v1_P_h" , "v_{1} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v2_P_h = new TH2D("v2_P_h" , "v_{2} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3); 
  TH2D *v3_P_h = new TH2D("v3_P_h" , "v_{3} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v4_P_h = new TH2D("v4_P_h" , "v_{4} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);
  TH2D *v5_P_h = new TH2D("v5_P_h" , "v_{5} vs p_{T} for p^{+}, p^{-}" , 100 , 0 , 6 , 100 , -0.1 , 0.3);

  TH1D *v1_sq_h = new TH1D("v1_sq_h" ,"v_{1}^{2} for all particles" , 1000 , 0 , 0.09 );
  TH1D *v2_sq_h = new TH1D("v2_sq_h" ,"v_{2}^{2} for all particles" , 1000 , 0 , 0.09 );
  TH1D *v3_sq_h = new TH1D("v3_sq_h" ,"v_{3}^{2} for all particles" , 1000 , 0 , 0.09 );
  TH1D *v4_sq_h = new TH1D("v4_sq_h" ,"v_{4}^{2} for all particles" , 1000 , 0 , 0.09 );
  TH1D *v5_sq_h = new TH1D("v5_sq_h" ,"v_{5}^{2} for all particles" , 1000 , 0 , 0.09 );

  TH1D *Psi_1_h = new TH1D("Psi_1_h" , "#Psi_{1} for all events thrown" , 100,-(TMath::Pi())/2,(5/2)*TMath::Pi() );
  TH1D *Psi_2_h = new TH1D("Psi_2_h" , "#Psi_{2} for all events thrown" , 100,-(TMath::Pi())/2,(5/2)*TMath::Pi() );
  TH1D *Psi_3_h = new TH1D("Psi_3_h" , "#Psi_{3} for all events thrown" , 100,-(TMath::Pi())/2,(5/2)*TMath::Pi() );
  TH1D *Psi_4_h = new TH1D("Psi_4_h" , "#Psi_{4} for all events thrown" , 100,-(TMath::Pi())/2,(5/2)*TMath::Pi() );
  TH1D *Psi_5_h = new TH1D("Psi_5_h" , "#Psi_{5} for all events thrown" , 100,-(TMath::Pi())/2,(5/2)*TMath::Pi() );

public:
  TennGen();
  virtual ~TennGen();
  void   SetCentralityBin(Int_t cb){CentralityBin = cb;}//inline function - nothing needed in .cxx
  void   SetRandomSeed(UInt_t seed){Seed = seed;}//inline function - nothing needed in .cxx
  void   GetRandomSeed();
  void   SetHarmonics(Int_t hf){HarmonicsFlag = hf;} //inline function - nothing needed in .cxx
  void   SetEtaRange(Double_t etarange){EtaRange = etarange;} //inline function - nothing needed in .cxx
  void   PrintOutQAHistos(Bool_t printQA = kFALSE){PRINTQA = printQA;} //inline funciton -nothing needed in .cxx
  void   InitializeBackground();
  void   CloserFunction(); //close all files and print QA histos to file if desired

  //TF1 *Random1();
  //TClonesArray* GetBackground() = (TClonesArray*) Rando -> GetRandom(Double_t xmin = 0. Double_t xmax = 1.);

   TClonesArray *GetBackground();
   ClassDef(TennGen,1)  //Event structure
};

#endif




















