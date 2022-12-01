
// Copyright (c) 2019; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// MODIFIED BY CHARLES HUGHES <chughe26@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2019-11-06 18:48 EST
//


#ifndef __CINT__
#include "TFile.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TMCParticle.h"
#include "TObject.h"
#include "TennGen.h"	
#include "TLegend.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TTimeStamp.h"
#include <chrono>
#include "TSystem.h"
#include "TROOT.h"
using namespace std::chrono;
#endif
		

class TennGen;

void TennGen_Old_Test(Int_t nEvents, Int_t jobID , Int_t HF , Int_t cent_bin, TString OutputDir ){

  
  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  //_____________________________________________________Creating Plots and Titles___________________________________________//

  char expression1[128];
  char expression2[128];
  char expression3[128];
  char expression4[128];
  char expression5[128];
  char expression6[128];
  char expression7[128];
  char expression8[128];
  char expression9[128];
  char expression10[128];
  char expression11[128];
  char expression12[128];
  char expression13[128];
  char expression14[128];
  char expression15[128];
  char expression16[128];
  char expression17[128];
  char expression18[128];
  char expression19[128];
  char expression20[128];
  char expression21[128];
  char expression22[128];
  char expression23[128];
  char expression24[128];
  char expression25[128];
  char expression26[128];
  char expression27[128];
  char expression28[128];
  char expression29[128];
  char expression30[128];
  char expression31[128];
  char expression32[128];
  char expression33[128];
  char expression34[128];
  char expression35[128];
  char expression36[128];
  char expression37[128];
  char expression38[128];
  char expression39[128];
  char expression40[128];
  char expression41[128];
  char expression42[128];

  sprintf(expression1 , "%s/TennGen_OLD_Test_File_HF=%d_cent_bin=%d.root" , OutputDir.Data() , HF , cent_bin );

  sprintf(expression2 , "p_{T} Distribution for #pi^{+}");
  sprintf(expression3 , "p_{T} Distribution for #pi^{-}");
  sprintf(expression4 , "p_{T} Distribution for #pi^{0}");
  sprintf(expression5 , "p_{T} Distribution for K^{+}");
  sprintf(expression6 , "p_{T} Distribution for K^{-}");
  sprintf(expression7 , "p_{T} Distribution for p^{+}");
  sprintf(expression8 , "p_{T} Distribution for p^{-}");
  sprintf(expression9 , "p_{T} Distribution for charged particles");
  sprintf(expression10 , "p_{T} Distribution for all particles");

  sprintf(expression11 , "#eta Distribution for #pi^{+}");
  sprintf(expression12 , "#eta Distribution for #pi^{-}");
  sprintf(expression13 , "#eta Distribution for #pi^{0}");
  sprintf(expression14 , "#eta Distribution for K^{+}");
  sprintf(expression15 , "#eta Distribution for K^{-}");
  sprintf(expression16 , "#eta Distribution for p^{+}");
  sprintf(expression17 , "#eta Distribution for p^{-}");
  sprintf(expression18 , "#eta Distribution for charged particles");
  sprintf(expression19 , "#eta Distribution for all particles");

  sprintf(expression20 , "#phi Distribution for #pi^{+}");
  sprintf(expression21 , "#phi Distribution for #pi^{-}");
  sprintf(expression22 , "#phi Distribution for #pi^{0}");
  sprintf(expression23 , "#phi Distribution for K^{+}");
  sprintf(expression24 , "#phi Distribution for K^{-}");
  sprintf(expression25 , "#phi Distribution for p^{+}");
  sprintf(expression26 , "#phi Distribution for p^{-}");
  sprintf(expression27 , "#phi Distribution for charged particles");
  sprintf(expression28 , "#phi Distribution for all particles");

  sprintf(expression29 , "#Psi_{EP,1} Distribution for all particles");
  sprintf(expression30 , "#Psi_{EP,2} Distribution for all particles");
  sprintf(expression31 , "#Psi_{EP,3} Distribution for all particles");
  sprintf(expression32 , "#Psi_{EP,4} Distribution for all particles");
  sprintf(expression33 , "#Psi_{EP,5} Distribution for all particles");

  sprintf(expression34 , "#Psi_{EP,1} Distribution for charged particles");
  sprintf(expression35 , "#Psi_{EP,2} Distribution for charged particles");
  sprintf(expression36 , "#Psi_{EP,3} Distribution for charged particles");
  sprintf(expression37 , "#Psi_{EP,4} Distribution for charged particles");
  sprintf(expression38 , "#Psi_{EP,5} Distribution for charged particles");

  sprintf(expression39 , "#phi vs #eta distribution of charged particles");
  sprintf(expression40 , "#phi vs #eta distribution of all particles");

  sprintf(expression41 , "#phi vs #eta distribution of charged particles weighted by p_{T}");
  sprintf(expression42 , "#phi vs #eta distribution of all particles weighted by p_{T}");



//___________________HISTOS_____________________________________________________________________________________//
//___________________START____BACKGROUND PARTICLE LEVEL HISTOS__________________________________________________//

  TH1D *histpT_piPlus = new TH1D("histpT_piPlus", expression2,200,0,10); //for piPlus
  histpT_piPlus -> Sumw2();
  histpT_piPlus->SetXTitle("p_{T} (GeV/c)");
  histpT_piPlus->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_piMinus = new TH1D("histpT_piMinus", expression3,200,0,10); //for piMinus
  histpT_piMinus -> Sumw2();
  histpT_piMinus->SetXTitle("p_{T} (GeV/c)");
  histpT_piMinus->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_piZero = new TH1D("histpT_piZero", expression4,200,0,10); //for piZero
  histpT_piZero -> Sumw2();
  histpT_piZero->SetXTitle("p_{T} (GeV/c)");
  histpT_piZero->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_kPlus = new TH1D("histpT_kPlus", expression5,200,0,10); //for kPlus
  histpT_kPlus -> Sumw2();
  histpT_kPlus->SetXTitle("p_{T} (GeV/c)");
  histpT_kPlus->SetYTitle("dN/dp_{T} ");
 
  TH1D *histpT_kMinus = new TH1D("histpT_kMinus", expression6,200,0,10); //for kMinus
  histpT_kMinus -> Sumw2();
  histpT_kMinus->SetXTitle("p_{T} (GeV/c)");
  histpT_kMinus->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_p = new TH1D("histpT_p", expression7,200,0,10); //for proton
  histpT_p -> Sumw2();
  histpT_p->SetXTitle("p_{T} (GeV/c)");
  histpT_p->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_pbar = new TH1D("histpT_pbar", expression8,200,0,10); //for anti-proton
  histpT_pbar -> Sumw2();
  histpT_pbar->SetXTitle("p_{T} (GeV/c)");
  histpT_pbar->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_ch = new TH1D("histpT_ch", expression9,200,0,10); //for all charged particles
  histpT_ch -> Sumw2();
  histpT_ch->SetXTitle("p_{T} (GeV/c)");
  histpT_ch->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_all = new TH1D("histpT_all", expression10,200,0,10); //for all particles (charged + neutral)
  histpT_all -> Sumw2();
  histpT_all->SetXTitle("p_{T} (GeV/c)");
  histpT_all->SetYTitle("dN/dp_{T} ");


  
  TH1D *histeta_piPlus = new TH1D("histeta_piPlus", expression11,200,-1.5,1.5); //for piPlus
  histeta_piPlus -> Sumw2();
  histeta_piPlus->SetXTitle("eta");
  histeta_piPlus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_piMinus = new TH1D("histeta_piMinus", expression12,200,-1.5,1.5); //for piMinus
  histeta_piMinus -> Sumw2();
  histeta_piMinus->SetXTitle("eta");
  histeta_piMinus->SetYTitle("dN/d#eta");

  TH1D *histeta_piZero = new TH1D("histeta_piZero", expression13,200,-1.5,1.5); //for piZero
  histeta_piZero -> Sumw2();
  histeta_piZero->SetXTitle("eta");
  histeta_piZero->SetYTitle("dN/d#eta");
  
  TH1D *histeta_kPlus = new TH1D("histeta_kPlus", expression14,200,-1.5,1.5); //for kPlus
  histeta_kPlus -> Sumw2();
  histeta_kPlus->SetXTitle("eta");
  histeta_kPlus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_kMinus = new TH1D("histeta_kMinus", expression15,200,-1.5,1.5); //for kMinus
  histeta_kMinus -> Sumw2();
  histeta_kMinus->SetXTitle("eta");
  histeta_kMinus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_p = new TH1D("histeta_p", expression16,200,-1.5,1.5); //for proton
  histeta_p -> Sumw2();
  histeta_p->SetXTitle("eta");
  histeta_p->SetYTitle("dN/d#eta");
  
  TH1D *histeta_pbar = new TH1D("histeta_pbar", expression17,200,-1.5,1.5); //for anti-proton
  histeta_pbar -> Sumw2();
  histeta_pbar->SetXTitle("eta");
  histeta_pbar->SetYTitle("dN/d#eta");

  TH1D *histeta_ch = new TH1D("histeta_ch", expression18,200,-1.5,1.5); //for all charged particle
  histeta_ch -> Sumw2();
  histeta_ch->SetXTitle("eta");
  histeta_ch->SetYTitle("dN/d#eta");

  TH1D *histeta_all = new TH1D("histeta_all", expression19,200,-1.5,1.5); //for all charged particle
  histeta_all -> Sumw2();
  histeta_all->SetXTitle("eta");
  histeta_all->SetYTitle("dN/d#eta");



  TH1D *histphi_piPlus = new TH1D("histphi_piPlus", expression20,200,0,2*TMath::Pi()); //for piPlus
  histphi_piPlus -> Sumw2();
  histphi_piPlus->SetXTitle("phi");
  histphi_piPlus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_piMinus = new TH1D("histphi_piMinus", expression21,200,0,2*TMath::Pi()); //for piMinus
  histphi_piMinus -> Sumw2();
  histphi_piMinus->SetXTitle("phi");
  histphi_piMinus->SetYTitle("dN/d#phi");

  TH1D *histphi_piZero = new TH1D("histphi_piZero", expression22,200,0,2*TMath::Pi()); //for piZero
  histphi_piZero -> Sumw2();
  histphi_piZero->SetXTitle("phi");
  histphi_piZero->SetYTitle("dN/d#phi");
  
  TH1D *histphi_kPlus = new TH1D("histphi_kPlus", expression23,200,0,2*TMath::Pi()); //for kPlus
  histphi_kPlus -> Sumw2();
  histphi_kPlus->SetXTitle("phi");
  histphi_kPlus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_kMinus = new TH1D("histphi_kMinus", expression24,200,0,2*TMath::Pi()); //for kMinus
  histphi_kMinus -> Sumw2();
  histphi_kMinus->SetXTitle("phi");
  histphi_kMinus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_p = new TH1D("histphi_p", expression25,200,0,2*TMath::Pi()); //for proton
  histphi_p -> Sumw2();
  histphi_p->SetXTitle("phi");
  histphi_p->SetYTitle("dN/d#phi");
  
  TH1D *histphi_pbar = new TH1D("histphi_pbar", expression26,200,0,2*TMath::Pi()); //for anti-proton
  histphi_pbar -> Sumw2();
  histphi_pbar->SetXTitle("phi");
  histphi_pbar->SetYTitle("dN/d#phi");

  TH1D *histphi_ch = new TH1D("histphi_ch", expression27,200,0,2*TMath::Pi()); //for all charged particles
  histphi_ch -> Sumw2();
  histphi_ch->SetXTitle("phi");
  histphi_ch->SetYTitle("dN/d#phi");

  TH1D *histphi_all = new TH1D("histphi_all", expression28,200,0,2*TMath::Pi()); //for all charged particles
  histphi_all -> Sumw2();
  histphi_all->SetXTitle("phi");
  histphi_all->SetYTitle("dN/d#phi");



  TH1D *histpsi_1_all = new TH1D("histpsi_1_all", expression29,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_1 for all particles
  histpsi_1_all -> Sumw2();
  histpsi_1_all->SetXTitle("#Psi_{EP,1} (radians)");
  histpsi_1_all->SetYTitle("dN/d#Psi_{EP,1}");

  TH1D *histpsi_2_all = new TH1D("histpsi_2_all", expression30,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_2 for all particles
  histpsi_2_all -> Sumw2();
  histpsi_2_all->SetXTitle("#Psi_{EP,2} (radians)");
  histpsi_2_all->SetYTitle("dN/d#Psi_{EP,2}");

  TH1D *histpsi_3_all = new TH1D("histpsi_3_all", expression31,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_3 for all particles
  histpsi_3_all -> Sumw2();
  histpsi_3_all->SetXTitle("#Psi_{EP,3} (radians)");
  histpsi_3_all->SetYTitle("dN/d#Psi_{EP,3}");

  TH1D *histpsi_4_all = new TH1D("histpsi_4_all", expression32,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_4 for all particles
  histpsi_4_all -> Sumw2();
  histpsi_4_all->SetXTitle("#Psi_{EP,4} (radians)");
  histpsi_4_all->SetYTitle("dN/d#Psi_{EP,4}");

  TH1D *histpsi_5_all = new TH1D("histpsi_5_all", expression33,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_5 for all particles
  histpsi_5_all -> Sumw2();
  histpsi_5_all->SetXTitle("#Psi_{EP,5} (radians)");
  histpsi_5_all->SetYTitle("dN/d#Psi_{EP,5}");



  TH1D *histpsi_1_ch = new TH1D("histpsi_1_ch", expression34,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_1 for charged particles
  histpsi_1_ch -> Sumw2();
  histpsi_1_ch ->SetXTitle("#Psi_{EP,1} (radians)");
  histpsi_1_ch ->SetYTitle("dN/d#Psi_{EP,1}");

  TH1D *histpsi_2_ch = new TH1D("histpsi_2_ch", expression35,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_2 for charged particles
  histpsi_2_ch -> Sumw2();
  histpsi_2_ch->SetXTitle("#Psi_{EP,2} (radians)");
  histpsi_2_ch->SetYTitle("dN/d#Psi_{EP,2}");

  TH1D *histpsi_3_ch = new TH1D("histpsi_3_ch", expression36,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_3 for charged particles
  histpsi_3_ch -> Sumw2();
  histpsi_3_ch ->SetXTitle("#Psi_{EP,3} (radians)");
  histpsi_3_ch ->SetYTitle("dN/d#Psi_{EP,3}");

  TH1D *histpsi_4_ch = new TH1D("histpsi_4_ch", expression37,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_4 for charged particles
  histpsi_4_ch -> Sumw2();
  histpsi_4_ch ->SetXTitle("#Psi_{EP,4} (radians)");
  histpsi_4_ch ->SetYTitle("dN/d#Psi_{EP,4}");

  TH1D *histpsi_5_ch = new TH1D("histpsi_5_ch", expression38,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_5 for charged particles
  histpsi_5_ch -> Sumw2();
  histpsi_5_ch ->SetXTitle("#Psi_{EP,5} (radians)");
  histpsi_5_ch ->SetYTitle("dN/d#Psi_{EP,5}");


  TH2D *histphi_eta_ch = new TH2D("histphi_eta_ch",expression39,200,0,2*TMath::Pi(),200,-1.1,1.1); //for charged particles
  histphi_eta_ch -> Sumw2();
  histphi_eta_ch->SetXTitle("#phi (radians)");
  histphi_eta_ch->SetYTitle("#eta (pseudo-rapidity)");
  histphi_eta_ch->SetZTitle("dN^{ch.}/d#phid#eta");

  TH2D *histphi_eta_all = new TH2D("histphi_eta_all",expression40,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
  histphi_eta_all -> Sumw2();
  histphi_eta_all->SetXTitle("#phi (radians)");
  histphi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
  histphi_eta_all->SetZTitle("dN^{ch.}/d#phid#eta");



  TH2D *hist_pT_phi_eta_ch = new TH2D("hist_pT_phi_eta_ch",expression41,200,0,2*TMath::Pi(),200,-1.1,1.1); //for charged particles
  hist_pT_phi_eta_ch -> Sumw2();
  hist_pT_phi_eta_ch->SetXTitle("#phi (radians)");
  hist_pT_phi_eta_ch->SetYTitle("#eta (pseudo-rapidity)");
  hist_pT_phi_eta_ch->SetZTitle("dp_{T}^{ch.}/d#phid#eta");

  TH2D *hist_pT_phi_eta_all = new TH2D("hist_pT_phi_eta_all",expression42,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
  hist_pT_phi_eta_all -> Sumw2();
  hist_pT_phi_eta_all->SetXTitle("#phi (radians)");
  hist_pT_phi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
  hist_pT_phi_eta_all->SetZTitle("dp_{T}^{ch.}/d#phid#eta");


//___________________END____BACKGROUND PARTICLE LEVEL HISTOS__________________________________________________//

  TCanvas *c1 = new TCanvas("c1" , "Run Time For 1 TennGen Event vs Event Number (Events Run in Sequence)", 200, 10,500,300); 

  Int_t* time_array;
  Int_t* time_interval_array;
  Int_t* event_array;

  time_array = new Int_t[nEvents + 1];
  time_interval_array = new Int_t[nEvents];
  event_array = new Int_t[nEvents];

  ///////////////////////////////////////////////////////////////////////
  
  //initialize counting variables 
  Int_t N_counter = 0 ; //counts events

  //The following block needs to be used on ACF Only

  TTimeStamp *timestamp = new TTimeStamp();
  UInt_t seed5 = ((jobID+9)*97) + (timestamp->GetSec() - static_cast<Int_t> (1446000000));
  //____________________________________________________________________
  if( (seed5>=0) && (seed5<=900000000) ) {
    //pythia->SetMRPY(1, ((rando->Uniform(1.,1000.))*(N_counter + 5 + seed)) );                   // set seed
    //pythia->SetMRPY(2, 0);                      // use new seed
   // TRandom *Rand = new TRandom(seed); //this I will use for the pT efficiency
    cout<<"Random Seed : "<<seed5<<endl;
  } else {cout << "error: time " << seed5 << " is not valid" << endl; exit(2);}
  

  //________Initiliaze Background Generator________________________________________________//
  
  TennGen *bkgd = new TennGen(); //constructor
  
  bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
  bkgd->SetRandomSeed(seed5); //setting the seed
  bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
  bkgd->SetEtaRange(0.9); //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
  bkgd->PrintOutQAHistos(kTRUE);
  bkgd->InitializeBackground();
     //TClonesArray *background =    bkgd->GetBackground();
  TClonesArray *BKGD = bkgd->GetBackground();
  TDatime *time = new TDatime(); //creating new time stamp
  cout<< "\n\nStart Time is " << time->GetTime() <<"\n\n\n" <<endl;
  time_array[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;

  
  //_____________________________________________Starting Event Loop_______________________________________________________________//
  for (Int_t event = 0; event < nEvents; event++) {
    if (event % 10 == 0) cout << "Event # " << event << endl;  
    N_counter ++;


      if( event != 0){
        free(time);
        TDatime *time = new TDatime(); //creating new time stamp
        cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
        time_array[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        //free(bkgd); //reset bkgd generator class variable
        //bkgd->Clear();
        //free(BKGD); //reset BKGD TClonesArray
        BKGD->Clear();
        //TennGen *bkgd = new TennGen();
        //bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
        //Int_t seed_b = (( TMath::CeilNint((rando->Uniform(1.,1000.)))));
        //bkgd->SetRandomSeed( jobID+event ); //setting the seed
        //cout<<"\n\n\n\n\n\njust set new seed\n\n\n\n\n"<<endl;
        //bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
        //bkgd->SetEtaRange(0.9); // set eta range (takes care of scaling in the .cxx) MUST BE UNIFORM ETA RANGE e.g. |eta| < 0.9 !!!
        //bkgd->InitializeBackground();
       //TClonesArray *background =    bkgd->GetBackground();
        TClonesArray *BKGD = bkgd->GetBackground();
      } 
    
      bkgd->GetRandomSeed(); 

      Double_t nParticles= BKGD ->GetEntries();
    

      //The following lines are used in the event plane calculations

      Double_t Sum_Psi_1_numer_all = 0;
      Double_t Sum_Psi_1_denom_all = 0;

      Double_t Sum_Psi_2_numer_all = 0;
      Double_t Sum_Psi_2_denom_all = 0;

      Double_t Sum_Psi_3_numer_all = 0;
      Double_t Sum_Psi_3_denom_all = 0;

      Double_t Sum_Psi_4_numer_all = 0;
      Double_t Sum_Psi_4_denom_all = 0;

      Double_t Sum_Psi_5_numer_all = 0;
      Double_t Sum_Psi_5_denom_all = 0;

      Double_t Sum_Psi_1_numer_ch = 0;
      Double_t Sum_Psi_1_denom_ch = 0;

      Double_t Sum_Psi_2_numer_ch = 0;
      Double_t Sum_Psi_2_denom_ch = 0;

      Double_t Sum_Psi_3_numer_ch = 0;
      Double_t Sum_Psi_3_denom_ch = 0;

      Double_t Sum_Psi_4_numer_ch = 0;
      Double_t Sum_Psi_4_denom_ch = 0;

      Double_t Sum_Psi_5_numer_ch = 0;
      Double_t Sum_Psi_5_denom_ch = 0;

    //____________________________Starting Background Particle Loop___________________________________________________//
      for (Int_t N = 0 ; N < nParticles ; N++){
        TMCParticle* bob = (TMCParticle*)BKGD->At(N); 
        Int_t K_F = bob -> GetKF();
        Double_t p = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) + ((bob->GetPz())*(bob->GetPz())) );
        Double_t pT = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) );
        Double_t eta = 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz())));
        Double_t phi = 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx()));

        //using q-vector method to reconstruct the event planes

        Sum_Psi_1_numer_all = Sum_Psi_1_numer_all + pT*TMath::Sin( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_1_denom_all = Sum_Psi_1_denom_all + pT*TMath::Cos( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_2_numer_all = Sum_Psi_2_numer_all + pT*TMath::Sin( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_2_denom_all = Sum_Psi_2_denom_all + pT*TMath::Cos( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_3_numer_all = Sum_Psi_3_numer_all + pT*TMath::Sin( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_3_denom_all = Sum_Psi_3_denom_all + pT*TMath::Cos( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_4_numer_all = Sum_Psi_4_numer_all + pT*TMath::Sin( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_4_denom_all = Sum_Psi_4_denom_all + pT*TMath::Cos( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_5_numer_all = Sum_Psi_5_numer_all + pT*TMath::Sin( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_5_denom_all = Sum_Psi_5_denom_all + pT*TMath::Cos( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        if(K_F != 111){ // ANYTHING BUT piPlus

          Sum_Psi_1_numer_ch = Sum_Psi_1_numer_ch + pT*TMath::Sin( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_1_denom_ch = Sum_Psi_1_denom_ch + pT*TMath::Cos( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_2_numer_ch = Sum_Psi_2_numer_ch + pT*TMath::Sin( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_2_denom_ch = Sum_Psi_2_denom_ch + pT*TMath::Cos( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_3_numer_ch = Sum_Psi_3_numer_ch + pT*TMath::Sin( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_3_denom_ch = Sum_Psi_3_denom_ch + pT*TMath::Cos( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
 
          Sum_Psi_4_numer_ch = Sum_Psi_4_numer_ch + pT*TMath::Sin( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_4_denom_ch = Sum_Psi_4_denom_ch + pT*TMath::Cos( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_5_numer_ch = Sum_Psi_5_numer_ch + pT*TMath::Sin( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_5_denom_ch = Sum_Psi_5_denom_ch + pT*TMath::Cos( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        }

        histpT_all->Fill(pT); //pT
        histphi_all->Fill( phi ); //azimuthal angle
        histeta_all->Fill(eta); //eta

        histphi_eta_all->Fill(phi,eta);
        hist_pT_phi_eta_all->Fill(phi,eta,pT);
      
        //cout<<"Particle "<<N<<" KF = "<< K_F<<endl;
        if(K_F == 211){ //piPlus
          histpT_piPlus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) ); //transverse momentum
          histeta_piPlus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) ); //pseudo-rapidity
          histphi_piPlus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ); //azimuthal angle 
          histpT_all->Fill(pT); //pT

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        }

        else if(K_F == -211){//piMinus
          histpT_piMinus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_piMinus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_piMinus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );     

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT); 
        } 

        else if(K_F == 111){//piZero
          histpT_piZero->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_piZero->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_piZero->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );
        }

        else if(K_F == 321){//kPlus
          histpT_kPlus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_kPlus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_kPlus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ); 

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        } 

        else if(K_F == -321){//kMinus
          histpT_kMinus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_kMinus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_kMinus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );   

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);   
        }

        else if(K_F == 2212){//p
          histpT_p->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_p->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_p->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );    

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        } 

        else if(K_F == -2212){//pbar
          histpT_pbar->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_pbar->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_pbar->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        }

     
      }//END PARTICLE LOOP BACKGROUND

      histpsi_1_all->Fill( (1.0/1.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_1_numer_all,-Sum_Psi_1_denom_all)));
      histpsi_2_all->Fill( (1.0/2.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_2_numer_all,-Sum_Psi_2_denom_all)));
      histpsi_3_all->Fill( (1.0/3.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_3_numer_all,-Sum_Psi_3_denom_all)));
      histpsi_4_all->Fill( (1.0/4.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_4_numer_all,-Sum_Psi_4_denom_all)));
      histpsi_5_all->Fill( (1.0/5.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_5_numer_all,-Sum_Psi_5_denom_all)));


      histpsi_1_ch->Fill( (1.0/1.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_1_numer_ch,-Sum_Psi_1_denom_ch)));
      histpsi_2_ch->Fill( (1.0/2.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_2_numer_ch,-Sum_Psi_2_denom_ch)));
      histpsi_3_ch->Fill( (1.0/3.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_3_numer_ch,-Sum_Psi_3_denom_ch)));
      histpsi_4_ch->Fill( (1.0/4.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_4_numer_ch,-Sum_Psi_4_denom_ch)));
      histpsi_5_ch->Fill( (1.0/5.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_5_numer_ch,-Sum_Psi_5_denom_ch)));
      

  }//________________________________________________END EVENT LOOP________________________________________

  bkgd->CloserFunction(); //close all the damn files, and print out QA histos to file if desired

  TFile *ff = new TFile(expression1,"RECREATE"); //make output files and directories

  TDirectory *particle_level_histos = ff->mkdir("particle_level_histos");
  TDirectory *all_charged_particles = ff->mkdir("all_charged_particles");
  TDirectory *all_particles = ff->mkdir("all_particles");

   TDatime *time2 = new TDatime(); //creating new time stamp
   cout << "\n\n\n\n\nTime at Event # "<< nEvents << " is " << time2->GetTime() <<"\n\n\n\n" << endl;
   time_array[nEvents] = ((time2->GetYear()*31540000) + (time2->GetMonth()*2628000) + (time2->GetDay()*86400) + (time2->GetHour()*3600) + (time2->GetMinute()*60) + time2->GetSecond()) ;

   //calculating the average time for one event

   Int_t Sum_Time = 0;
   Int_t Total_Time = 0;

   for( Int_t i = 0 ; i < nEvents; i++){
     time_interval_array[i] = time_array[i + 1] - time_array[i];
     Sum_Time = Sum_Time + time_interval_array[i];
     event_array[i] = i + 1;
   }

   cout<< "\n\n\n\n Total Time for = "<<nEvents<<" Events = "<<Sum_Time<<" Seconds \n\n\n\n"<<endl;

   cout<< "\n\n\n\n Average Time For One Event = "<<Sum_Time/(nEvents)<<" Seconds \n\n\n\n"<<endl;

  TH1D *histruntime = new TH1D("histruntime", "Run Times for Background Generator",8,((Sum_Time/nEvents) - 1),((Sum_Time/nEvents) + 1)); //for non-random distribution
  histruntime -> Sumw2();
  histruntime->SetXTitle("Run Time For 1 Event (Seconds)");
  histruntime->SetYTitle("Counts");

  TGraph *rtgraph = new TGraph(nEvents , event_array , time_interval_array);
  rtgraph->SetTitle("Run Time For 1 Background Event vs Event Number (Events Run in Sequence)");
  rtgraph->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
  rtgraph->Write();

  c1->cd(1); 
  rtgraph->Draw("AC");
  c1->Draw();
  c1->Write("");
  //c1->Close();

  for( Int_t j = 0 ; j < nEvents; j++){
    histruntime->Fill(time_interval_array[j]);
  }

  histruntime->Write();
   //________________________________Setting Histogram Styles and Writing to output file_______________________________________//

   histpT_piPlus->SetMarkerColor(2);
   histpT_piPlus->SetMarkerStyle(20);
   histpT_piMinus->SetMarkerColor(3);
   histpT_piMinus->SetMarkerStyle(21);
   histpT_kPlus->SetMarkerColor(4);
   histpT_kPlus->SetMarkerStyle(22);
   histpT_kMinus->SetMarkerColor(5);
   histpT_kMinus->SetMarkerStyle(23);
   histpT_p->SetMarkerColor(6);
   histpT_p->SetMarkerStyle(43);
   histpT_pbar->SetMarkerColor(7);
   histpT_pbar->SetMarkerStyle(47);
   histpT_piZero->SetMarkerColor(13);
   histpT_piZero->SetMarkerStyle(34);

   histeta_piPlus->SetMarkerColor(2);
   histeta_piPlus->SetMarkerStyle(20);
   histeta_piMinus->SetMarkerColor(3);
   histeta_piMinus->SetMarkerStyle(21);
   histeta_kPlus->SetMarkerColor(4);
   histeta_kPlus->SetMarkerStyle(22);
   histeta_kMinus->SetMarkerColor(5);
   histeta_kMinus->SetMarkerStyle(23);
   histeta_p->SetMarkerColor(6);
   histeta_p->SetMarkerStyle(43);
   histeta_pbar->SetMarkerColor(7);
   histeta_pbar->SetMarkerStyle(47);
   histeta_piZero->SetMarkerColor(13);
   histeta_piZero->SetMarkerStyle(34);

   histphi_piPlus->SetMarkerColor(2);
   histphi_piPlus->SetMarkerStyle(20);
   histphi_piMinus->SetMarkerColor(3);
   histphi_piMinus->SetMarkerStyle(21);
   histphi_kPlus->SetMarkerColor(4);
   histphi_kPlus->SetMarkerStyle(22);
   histphi_kMinus->SetMarkerColor(5);
   histphi_kMinus->SetMarkerStyle(23);
   histphi_p->SetMarkerColor(6);
   histphi_p->SetMarkerStyle(43);
   histphi_pbar->SetMarkerColor(7);
   histphi_pbar->SetMarkerStyle(47);
   histphi_piZero->SetMarkerColor(13);
   histphi_piZero->SetMarkerStyle(34);


   particle_level_histos->cd();

   histpT_piPlus->Write();
   histpT_piMinus->Write();
   histpT_kPlus->Write();
   histpT_kMinus->Write();
   histpT_p->Write();
   histpT_pbar->Write();
   histpT_piZero->Write();

   histeta_piPlus->Write();
   histeta_piMinus->Write();
   histeta_kPlus->Write();
   histeta_kMinus->Write();
   histeta_p->Write();
   histeta_pbar->Write();
   histeta_piZero->Write();

   histphi_piPlus->Write();
   histphi_piMinus->Write();
   histphi_kPlus->Write();
   histphi_kMinus->Write();
   histphi_p->Write();
   histphi_pbar->Write();
   histphi_piZero->Write();

   all_charged_particles->cd();

   histpT_ch->Write();
   histeta_ch->Write();
   histphi_ch->Write();

   histpsi_1_ch->Write();
   histpsi_2_ch->Write();
   histpsi_3_ch->Write();
   histpsi_4_ch->Write();
   histpsi_5_ch->Write();

   histphi_eta_ch->Write();
   hist_pT_phi_eta_ch->Write();

   all_particles->cd();

   histpT_all->Write();
   histeta_all->Write();
   histphi_all->Write();
   
   histpsi_1_all->Write();
   histpsi_2_all->Write();
   histpsi_3_all->Write();
   histpsi_4_all->Write();
   histpsi_5_all->Write();

   histphi_eta_all->Write();
   hist_pT_phi_eta_all->Write();

 

//__________________________________________________________________________________________________________________________________//
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe <pT> for pi_plus = "<<histpT_piPlus->GetMean(1)<<" GeV/c, The Sigma for pi_plus = "<<histpT_piPlus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for pi_minus = "<<histpT_piMinus->GetMean(1)<<" GeV/c, The Sigma for pi_minus = "<<histpT_piMinus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for pi_zero = "<<histpT_piZero->GetMean(1)<<" GeV/c, The Sigma for pi_zero = "<<histpT_piZero->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for k_plus = "<<histpT_kPlus->GetMean(1)<<" GeV/c, The Sigma for k_plus = "<<histpT_kPlus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for k_minus = "<<histpT_kMinus->GetMean(1)<<" GeV/c, The Sigma for k_minus = "<<histpT_kMinus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for p = "<<histpT_p->GetMean(1)<<" GeV/c, The Sigma for p = "<<histpT_p->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for p_bar = "<<histpT_pbar->GetMean(1)<<" GeV/c, The Sigma for p_bar = "<<histpT_pbar->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe <pT> for all = "<<histpT_all->GetMean(1)<<" GeV/c, The Sigma for all = "<<histpT_all->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for ch = "<<histpT_ch->GetMean(1)<<" GeV/c, The Sigma for ch = "<<histpT_ch->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe Nall/event = "<<histphi_all->GetEntries()/nEvents<<" particles/event "<<endl;
   cout<<"\nThe Nch/event  = "<<histphi_ch->GetEntries()/nEvents<<" charged particles/event "<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;	
//__________________________________________________________________________________________________________________________________//

 

   
  ff -> Write();
  ff -> Close();
  bkgd->~TennGen();
} //end Background_Test_Macro.C

void Time_Compare_Tests(Int_t nEvents, Int_t jobID , Int_t HF , Int_t cent_bin, TString OutputDir  ){

    
  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  //_____________________________________________________Creating Plots and Titles___________________________________________//

  Int_t* time_array;
  Int_t* time_array_vonNeumann;
  Int_t* time_interval_array;
  Int_t* time_interval_array_vonNeumann;
  Int_t* event_array;
  Int_t* event_array_vonNeumann;

  time_array = new Int_t[nEvents + 1];
  time_interval_array = new Int_t[nEvents];
  event_array = new Int_t[nEvents];

  time_array_vonNeumann = new Int_t[nEvents + 1];
  time_interval_array_vonNeumann = new Int_t[nEvents];
  event_array_vonNeumann = new Int_t[nEvents];



  ///////////////////////////////////////////////////////////////////////
  
  //initialize counting variables 
  Int_t N_counter = 0 ; //counts events
  Int_t N_counter_von = 0; 
  //The following block needs to be used on ACF Only

  TTimeStamp *timestamp = new TTimeStamp();
  UInt_t seed5 = ((jobID+9)*97) + (timestamp->GetSec() - static_cast<Int_t> (1446000000));
  //____________________________________________________________________
  if( (seed5>=0) && (seed5<=900000000) ) {
    //pythia->SetMRPY(1, ((rando->Uniform(1.,1000.))*(N_counter + 5 + seed)) );                   // set seed
    //pythia->SetMRPY(2, 0);                      // use new seed
   // TRandom *Rand = new TRandom(seed); //this I will use for the pT efficiency
    cout<<"Random Seed : "<<seed5<<endl;
  } else {cout << "error: time " << seed5 << " is not valid" << endl; exit(2);}
  

  //________Initiliaze Background Generator________________________________________________//
  
  TennGen *bkgd = new TennGen(); //constructor
  
  bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
  bkgd->SetRandomSeed(seed5); //setting the seed
  bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
  bkgd->SetEtaRange(0.9); //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
  bkgd->PrintOutQAHistos(kTRUE);
  bkgd->InitializeBackground();
     //TClonesArray *background =    bkgd->GetBackground();
  TClonesArray *BKGD = bkgd->GetBackground();
  TDatime *time = new TDatime(); //creating new time stamp
  cout<< "\n\nStart Time is " << time->GetTime() <<"\n\n\n" <<endl;
  time_array[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
   
 // _____________________________________________Starting Event Loop_______________________________________________________________//
  // for (Int_t event = 0; event < nEvents; event++) {
  //   if (event % 10 == 0) cout << "Event # " << event << endl;  
  //   N_counter ++;

  //     if( event != 0){
  //       free(time);
  //       TDatime *time = new TDatime(); //creating new time stamp
  //       cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
  //       time_array[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
  //       //free(bkgd); //reset bkgd generator class variable
  //       //bkgd->Clear();
  //       //free(BKGD); //reset BKGD TClonesArray
        
  //       BKGD->Clear();
  //       //TennGen *bkgd = new TennGen();
  //       //bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
  //       //Int_t seed_b = (( TMath::CeilNint((rando->Uniform(1.,1000.)))));
  //       //bkgd->SetRandomSeed( jobID+event ); //setting the seed
  //       //cout<<"\n\n\n\n\n\njust set new seed\n\n\n\n\n"<<endl;
  //       //bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
  //       //bkgd->SetEtaRange(0.9); // set eta range (takes care of scaling in the .cxx) MUST BE UNIFORM ETA RANGE e.g. |eta| < 0.9 !!!
  //       //bkgd->InitializeBackground();
  //      //TClonesArray *background =    bkgd->GetBackground();
  //       TClonesArray *BKGD = bkgd->GetBackground();
       
  //     } 

  // }

  TDatime *time2 = new TDatime(); //creating new time stamp
  cout << "\n\n\n\n\nTime at Event # "<< nEvents << " is " << time2->GetTime() <<"\n\n\n\n" << endl;
  time_array[nEvents] = ((time2->GetYear()*31540000) + (time2->GetMonth()*2628000) + (time2->GetDay()*86400) + (time2->GetHour()*3600) + (time2->GetMinute()*60) + time2->GetSecond()) ;


  time_array_vonNeumann[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
  for (Int_t event = 0; event < nEvents; event++) {
    if (event % 10 == 0) cout << "Event # " << event << endl;  
    N_counter ++;

      if( event != 0){
        free(time);
        TDatime *time = new TDatime(); //creating new time stamp
        cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
        time_array_vonNeumann[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        //free(bkgd); //reset bkgd generator class variable
        //bkgd->Clear();
        //free(BKGD); //reset BKGD TClonesArray
        auto start = high_resolution_clock::now();
        BKGD->Clear();
        //TennGen *bkgd = new TennGen();
        //bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
        //Int_t seed_b = (( TMath::CeilNint((rando->Uniform(1.,1000.)))));
        //bkgd->SetRandomSeed( jobID+event ); //setting the seed
        //cout<<"\n\n\n\n\n\njust set new seed\n\n\n\n\n"<<endl;
        //bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
        //bkgd->SetEtaRange(0.9); // set eta range (takes care of scaling in the .cxx) MUST BE UNIFORM ETA RANGE e.g. |eta| < 0.9 !!!
        //bkgd->InitializeBackground();
       //TClonesArray *background =    bkgd->GetBackground();
        TClonesArray *BKGD = bkgd->GetBackground_TannerTest();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        time_interval_array_vonNeumann[event] = Int_t(duration.count());
      } 
      
  }

  TDatime *time3 = new TDatime(); //creating new time stamp
  cout << "\n\n\n\n\nTime at Event # "<< nEvents << " is " << time3->GetTime() <<"\n\n\n\n" << endl;
  time_array_vonNeumann[nEvents] = ((time3->GetYear()*31540000) + (time3->GetMonth()*2628000) + (time3->GetDay()*86400) + (time3->GetHour()*3600) + (time3->GetMinute()*60) + time3->GetSecond()) ;

   //calculating the average time for one event

   Int_t Sum_Time = 0;
   Int_t Total_Time_Von = time_array_vonNeumann[nEvents] - time_array_vonNeumann[0];
   Int_t TOTAL_NORM_TIME[nEvents];
   Int_t TOTAL_NORM_VON[nEvents];
   for( Int_t i = 0 ; i < nEvents; i++){
     time_interval_array[i] = time_array[i + 1] - time_array[i];
     //time_interval_array_vonNeumann[i] = time_array_vonNeumann[i + 1] - time_array_vonNeumann[i];
     Sum_Time = Sum_Time + time_interval_array[i];
     Total_Time_Von = Total_Time_Von +time_interval_array_vonNeumann[i];
     cout << "Time interval for event # " << i << " is " << time_interval_array_vonNeumann[i] << endl;
     event_array[i] = i + 1;
     event_array_vonNeumann[i] = i+1;

     TOTAL_NORM_TIME[i] = Sum_Time;
     TOTAL_NORM_VON[i] = Total_Time_Von;
   }

   cout<< "\n\n\n\n Total Time for = "<<nEvents<<" Events = "<<Total_Time_Von<<" Seconds \n\n\n\n"<<endl;

   cout<< "\n\n\n\n Average Time For One Event = "<<Total_Time_Von/(nEvents)<<" Seconds \n\n\n\n"<<endl;
  Double_t average_norm = (Sum_Time/nEvents);
  Double_t average_von = (Total_Time_Von/nEvents);
  Double_t minX_norm =10000.0;
  Double_t minX_von = 10000.0; 
  Double_t maxX_norm = 0.0;
  Double_t maxX_von = 0.0;
  for( Int_t i = 0 ; i < nEvents; i++){
    if(time_interval_array[i] < minX_norm) minX_norm = time_interval_array[i];
    if(time_interval_array[i] > maxX_norm) maxX_norm = time_interval_array[i];
    if(time_interval_array_vonNeumann[i] < minX_von) minX_von = time_interval_array_vonNeumann[i];
    if(time_interval_array_vonNeumann[i] > maxX_von) maxX_von = time_interval_array_vonNeumann[i];
  }

  minX_norm*=0.2;
  minX_von*= 0.2;
  maxX_von*=1.2;
  maxX_norm*=1.2;


  TH1D *histruntime_norm = new TH1D("histruntime_norm", "Run Times for Background Generator",100, minX_norm ,maxX_norm); //for non-random distribution
  histruntime_norm -> Sumw2();
  histruntime_norm->SetXTitle("Run Time For 1 Event (Seconds)");
  histruntime_norm->SetYTitle("Counts");

  TH1D *histruntime_von = new TH1D("histruntime_von", "Run Times for Background Generator",100, minX_von ,maxX_von); //for non-random distribution
  histruntime_von -> Sumw2();
  histruntime_von->SetXTitle("Run Time For 1 Event with Von Nuemann Sampling (Seconds)");
  histruntime_von->SetYTitle("Counts");

  for( Int_t j = 0 ; j < nEvents; j++){
    histruntime_norm->Fill(time_interval_array[j]);
    histruntime_von->Fill(time_interval_array_vonNeumann[j]);
  }

  TCanvas *c1 = new TCanvas("c1" , "Run Time For 1 Event (Seconds)", 200, 10,500,300); 
  c1->cd(1); 
  
  histruntime_norm->SetLineColor(kRed);
  histruntime_norm->SetLineWidth(2);
  histruntime_norm->Draw("HIST");


  histruntime_von->SetLineColor(kBlue);
  histruntime_von->SetLineWidth(2);
  histruntime_von->Draw("HIST SAME");

  TLegend *leg = new TLegend(0.6,0.65,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  leg->SetTextFont(42);
  leg->SetNColumns(1);
  
  leg->AddEntry(histruntime_von,"Von Neumann Sampling (Microseconds)", "l");
  leg->AddEntry(histruntime_norm,"TF1 Sampling", "l");
  leg->Draw();  
  c1->SaveAs(Form("%s/RunTimeFor1Event.png",OutputDir.Data()));


  


  TGraph *rtgraph = new TGraph(nEvents , event_array , time_interval_array);
  rtgraph->SetTitle("Run Time For 1 Background Event vs Event Number (Events Run in Sequence)");
  rtgraph->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
  rtgraph->Write();

  TGraph *rtgraph_von = new TGraph(nEvents , event_array_vonNeumann , time_interval_array_vonNeumann);
  rtgraph_von->SetTitle("Run Time For 1 Background Event with Von Nuemann Sampling vs Event Number (Events Run in Sequence)");
  rtgraph_von->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph_von->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
  rtgraph_von->Write();


  TCanvas *c2 = new TCanvas("c2" , "Run Time For 1 Background Event vs Event Number (Events Run in Sequence)", 200, 10,500,300); 
  c2->cd(1); 

  rtgraph->SetLineColor(kRed);
  rtgraph->SetLineWidth(2);
  rtgraph->Draw("AC");

  rtgraph_von->SetLineColor(kBlue);
  rtgraph_von->SetLineWidth(2);
  rtgraph_von->Draw("C SAME");

  TLegend *leg1 = new TLegend(0.6,0.65,0.8,0.8);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.06);
  leg1->SetTextFont(42);
  leg1->SetNColumns(1);
  
  leg1->AddEntry(rtgraph_von,"Von Neumann Sampling (Microseconds)", "l");
  leg1->AddEntry(rtgraph,"TF1 Sampling", "l");
  leg1->Draw();  
  c2->SaveAs(Form("%s/RunTimeFor1EventvsEventNumber.png",OutputDir.Data()));

  TGraph *rtgraph_tot = new TGraph(nEvents , event_array , TOTAL_NORM_TIME);
  rtgraph_tot->SetTitle("Total Run Time vs Event Number (Events Run in Sequence)");
  rtgraph_tot->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph_tot->GetYaxis()->SetTitle("Total Time");
  rtgraph_tot->Write();

  TGraph *rtgraph_tot_von = new TGraph(nEvents , event_array_vonNeumann , TOTAL_NORM_VON);
  rtgraph_tot_von->SetTitle("Total Run Time with Von Nuemann Sampling vs Event Number (Events Run in Sequence)");
  rtgraph_tot_von->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph_tot_von->GetYaxis()->SetTitle("Total Time");
  rtgraph_tot_von->Write();


  TCanvas *c3 = new TCanvas("c3" , "Total Run Time vs Event Number (Events Run in Sequence)", 200, 10,500,300); 
  c3->cd(1); 

  rtgraph_tot->SetLineColor(kRed);
  rtgraph_tot->SetLineWidth(2);
  rtgraph_tot->Draw("AC");

  rtgraph_tot_von->SetLineColor(kBlue);
  rtgraph_tot_von->SetLineWidth(2);
  rtgraph_tot_von->Draw("C SAME");

  TLegend *leg2 = new TLegend(0.6,0.65,0.8,0.8);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.06);
  leg2->SetTextFont(42);
  leg2->SetNColumns(1);
  
  leg2->AddEntry(rtgraph_von,"Von Neumann Sampling (Microseconds)", "l");
  leg2->AddEntry(rtgraph,"TF1 Sampling", "l");
  leg2->Draw();  
  c3->SaveAs(Form("%s/TotalRunTimevsEventNumber.png",OutputDir.Data()));

  // c1->cd(1); 
  // rtgraph->Draw("AC");
  // c1->Draw();
  // c1->Write("");
  //c1->Close();

  c1->Close();
  c2->Close();
  c3->Close();


}

// nEvents is how many events we want. This is dynamic user input to be run from the command prompt and user specified

void TennGen_New_Test(Int_t nEvents, Int_t jobID , Int_t HF , Int_t cent_bin, TString OutputDir ){

  
  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  //_____________________________________________________Creating Plots and Titles___________________________________________//

  char expression1[128];
  char expression2[128];
  char expression3[128];
  char expression4[128];
  char expression5[128];
  char expression6[128];
  char expression7[128];
  char expression8[128];
  char expression9[128];
  char expression10[128];
  char expression11[128];
  char expression12[128];
  char expression13[128];
  char expression14[128];
  char expression15[128];
  char expression16[128];
  char expression17[128];
  char expression18[128];
  char expression19[128];
  char expression20[128];
  char expression21[128];
  char expression22[128];
  char expression23[128];
  char expression24[128];
  char expression25[128];
  char expression26[128];
  char expression27[128];
  char expression28[128];
  char expression29[128];
  char expression30[128];
  char expression31[128];
  char expression32[128];
  char expression33[128];
  char expression34[128];
  char expression35[128];
  char expression36[128];
  char expression37[128];
  char expression38[128];
  char expression39[128];
  char expression40[128];
  char expression41[128];
  char expression42[128];

  sprintf(expression1 , "%s/TennGen_VonNeumann_Test_File_HF=%d_cent_bin=%d.root" , OutputDir.Data() , HF , cent_bin);

  sprintf(expression2 , "p_{T} Distribution for #pi^{+}");
  sprintf(expression3 , "p_{T} Distribution for #pi^{-}");
  sprintf(expression4 , "p_{T} Distribution for #pi^{0}");
  sprintf(expression5 , "p_{T} Distribution for K^{+}");
  sprintf(expression6 , "p_{T} Distribution for K^{-}");
  sprintf(expression7 , "p_{T} Distribution for p^{+}");
  sprintf(expression8 , "p_{T} Distribution for p^{-}");
  sprintf(expression9 , "p_{T} Distribution for charged particles");
  sprintf(expression10 , "p_{T} Distribution for all particles");

  sprintf(expression11 , "#eta Distribution for #pi^{+}");
  sprintf(expression12 , "#eta Distribution for #pi^{-}");
  sprintf(expression13 , "#eta Distribution for #pi^{0}");
  sprintf(expression14 , "#eta Distribution for K^{+}");
  sprintf(expression15 , "#eta Distribution for K^{-}");
  sprintf(expression16 , "#eta Distribution for p^{+}");
  sprintf(expression17 , "#eta Distribution for p^{-}");
  sprintf(expression18 , "#eta Distribution for charged particles");
  sprintf(expression19 , "#eta Distribution for all particles");

  sprintf(expression20 , "#phi Distribution for #pi^{+}");
  sprintf(expression21 , "#phi Distribution for #pi^{-}");
  sprintf(expression22 , "#phi Distribution for #pi^{0}");
  sprintf(expression23 , "#phi Distribution for K^{+}");
  sprintf(expression24 , "#phi Distribution for K^{-}");
  sprintf(expression25 , "#phi Distribution for p^{+}");
  sprintf(expression26 , "#phi Distribution for p^{-}");
  sprintf(expression27 , "#phi Distribution for charged particles");
  sprintf(expression28 , "#phi Distribution for all particles");

  sprintf(expression29 , "#Psi_{EP,1} Distribution for all particles");
  sprintf(expression30 , "#Psi_{EP,2} Distribution for all particles");
  sprintf(expression31 , "#Psi_{EP,3} Distribution for all particles");
  sprintf(expression32 , "#Psi_{EP,4} Distribution for all particles");
  sprintf(expression33 , "#Psi_{EP,5} Distribution for all particles");

  sprintf(expression34 , "#Psi_{EP,1} Distribution for charged particles");
  sprintf(expression35 , "#Psi_{EP,2} Distribution for charged particles");
  sprintf(expression36 , "#Psi_{EP,3} Distribution for charged particles");
  sprintf(expression37 , "#Psi_{EP,4} Distribution for charged particles");
  sprintf(expression38 , "#Psi_{EP,5} Distribution for charged particles");

  sprintf(expression39 , "#phi vs #eta distribution of charged particles");
  sprintf(expression40 , "#phi vs #eta distribution of all particles");

  sprintf(expression41 , "#phi vs #eta distribution of charged particles weighted by p_{T}");
  sprintf(expression42 , "#phi vs #eta distribution of all particles weighted by p_{T}");



//___________________HISTOS_____________________________________________________________________________________//
//___________________START____BACKGROUND PARTICLE LEVEL HISTOS__________________________________________________//

  TH1D *histpT_piPlus = new TH1D("histpT_piPlus", expression2,200,0,10); //for piPlus
  histpT_piPlus -> Sumw2();
  histpT_piPlus->SetXTitle("p_{T} (GeV/c)");
  histpT_piPlus->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_piMinus = new TH1D("histpT_piMinus", expression3,200,0,10); //for piMinus
  histpT_piMinus -> Sumw2();
  histpT_piMinus->SetXTitle("p_{T} (GeV/c)");
  histpT_piMinus->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_piZero = new TH1D("histpT_piZero", expression4,200,0,10); //for piZero
  histpT_piZero -> Sumw2();
  histpT_piZero->SetXTitle("p_{T} (GeV/c)");
  histpT_piZero->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_kPlus = new TH1D("histpT_kPlus", expression5,200,0,10); //for kPlus
  histpT_kPlus -> Sumw2();
  histpT_kPlus->SetXTitle("p_{T} (GeV/c)");
  histpT_kPlus->SetYTitle("dN/dp_{T} ");
 
  TH1D *histpT_kMinus = new TH1D("histpT_kMinus", expression6,200,0,10); //for kMinus
  histpT_kMinus -> Sumw2();
  histpT_kMinus->SetXTitle("p_{T} (GeV/c)");
  histpT_kMinus->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_p = new TH1D("histpT_p", expression7,200,0,10); //for proton
  histpT_p -> Sumw2();
  histpT_p->SetXTitle("p_{T} (GeV/c)");
  histpT_p->SetYTitle("dN/dp_{T} ");
  
  TH1D *histpT_pbar = new TH1D("histpT_pbar", expression8,200,0,10); //for anti-proton
  histpT_pbar -> Sumw2();
  histpT_pbar->SetXTitle("p_{T} (GeV/c)");
  histpT_pbar->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_ch = new TH1D("histpT_ch", expression9,200,0,10); //for all charged particles
  histpT_ch -> Sumw2();
  histpT_ch->SetXTitle("p_{T} (GeV/c)");
  histpT_ch->SetYTitle("dN/dp_{T} ");

  TH1D *histpT_all = new TH1D("histpT_all", expression10,200,0,10); //for all particles (charged + neutral)
  histpT_all -> Sumw2();
  histpT_all->SetXTitle("p_{T} (GeV/c)");
  histpT_all->SetYTitle("dN/dp_{T} ");


  
  TH1D *histeta_piPlus = new TH1D("histeta_piPlus", expression11,200,-1.5,1.5); //for piPlus
  histeta_piPlus -> Sumw2();
  histeta_piPlus->SetXTitle("eta");
  histeta_piPlus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_piMinus = new TH1D("histeta_piMinus", expression12,200,-1.5,1.5); //for piMinus
  histeta_piMinus -> Sumw2();
  histeta_piMinus->SetXTitle("eta");
  histeta_piMinus->SetYTitle("dN/d#eta");

  TH1D *histeta_piZero = new TH1D("histeta_piZero", expression13,200,-1.5,1.5); //for piZero
  histeta_piZero -> Sumw2();
  histeta_piZero->SetXTitle("eta");
  histeta_piZero->SetYTitle("dN/d#eta");
  
  TH1D *histeta_kPlus = new TH1D("histeta_kPlus", expression14,200,-1.5,1.5); //for kPlus
  histeta_kPlus -> Sumw2();
  histeta_kPlus->SetXTitle("eta");
  histeta_kPlus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_kMinus = new TH1D("histeta_kMinus", expression15,200,-1.5,1.5); //for kMinus
  histeta_kMinus -> Sumw2();
  histeta_kMinus->SetXTitle("eta");
  histeta_kMinus->SetYTitle("dN/d#eta");
  
  TH1D *histeta_p = new TH1D("histeta_p", expression16,200,-1.5,1.5); //for proton
  histeta_p -> Sumw2();
  histeta_p->SetXTitle("eta");
  histeta_p->SetYTitle("dN/d#eta");
  
  TH1D *histeta_pbar = new TH1D("histeta_pbar", expression17,200,-1.5,1.5); //for anti-proton
  histeta_pbar -> Sumw2();
  histeta_pbar->SetXTitle("eta");
  histeta_pbar->SetYTitle("dN/d#eta");

  TH1D *histeta_ch = new TH1D("histeta_ch", expression18,200,-1.5,1.5); //for all charged particle
  histeta_ch -> Sumw2();
  histeta_ch->SetXTitle("eta");
  histeta_ch->SetYTitle("dN/d#eta");

  TH1D *histeta_all = new TH1D("histeta_all", expression19,200,-1.5,1.5); //for all charged particle
  histeta_all -> Sumw2();
  histeta_all->SetXTitle("eta");
  histeta_all->SetYTitle("dN/d#eta");



  TH1D *histphi_piPlus = new TH1D("histphi_piPlus", expression20,200,0,2*TMath::Pi()); //for piPlus
  histphi_piPlus -> Sumw2();
  histphi_piPlus->SetXTitle("phi");
  histphi_piPlus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_piMinus = new TH1D("histphi_piMinus", expression21,200,0,2*TMath::Pi()); //for piMinus
  histphi_piMinus -> Sumw2();
  histphi_piMinus->SetXTitle("phi");
  histphi_piMinus->SetYTitle("dN/d#phi");

  TH1D *histphi_piZero = new TH1D("histphi_piZero", expression22,200,0,2*TMath::Pi()); //for piZero
  histphi_piZero -> Sumw2();
  histphi_piZero->SetXTitle("phi");
  histphi_piZero->SetYTitle("dN/d#phi");
  
  TH1D *histphi_kPlus = new TH1D("histphi_kPlus", expression23,200,0,2*TMath::Pi()); //for kPlus
  histphi_kPlus -> Sumw2();
  histphi_kPlus->SetXTitle("phi");
  histphi_kPlus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_kMinus = new TH1D("histphi_kMinus", expression24,200,0,2*TMath::Pi()); //for kMinus
  histphi_kMinus -> Sumw2();
  histphi_kMinus->SetXTitle("phi");
  histphi_kMinus->SetYTitle("dN/d#phi");
  
  TH1D *histphi_p = new TH1D("histphi_p", expression25,200,0,2*TMath::Pi()); //for proton
  histphi_p -> Sumw2();
  histphi_p->SetXTitle("phi");
  histphi_p->SetYTitle("dN/d#phi");
  
  TH1D *histphi_pbar = new TH1D("histphi_pbar", expression26,200,0,2*TMath::Pi()); //for anti-proton
  histphi_pbar -> Sumw2();
  histphi_pbar->SetXTitle("phi");
  histphi_pbar->SetYTitle("dN/d#phi");

  TH1D *histphi_ch = new TH1D("histphi_ch", expression27,200,0,2*TMath::Pi()); //for all charged particles
  histphi_ch -> Sumw2();
  histphi_ch->SetXTitle("phi");
  histphi_ch->SetYTitle("dN/d#phi");

  TH1D *histphi_all = new TH1D("histphi_all", expression28,200,0,2*TMath::Pi()); //for all charged particles
  histphi_all -> Sumw2();
  histphi_all->SetXTitle("phi");
  histphi_all->SetYTitle("dN/d#phi");



  TH1D *histpsi_1_all = new TH1D("histpsi_1_all", expression29,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_1 for all particles
  histpsi_1_all -> Sumw2();
  histpsi_1_all->SetXTitle("#Psi_{EP,1} (radians)");
  histpsi_1_all->SetYTitle("dN/d#Psi_{EP,1}");

  TH1D *histpsi_2_all = new TH1D("histpsi_2_all", expression30,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_2 for all particles
  histpsi_2_all -> Sumw2();
  histpsi_2_all->SetXTitle("#Psi_{EP,2} (radians)");
  histpsi_2_all->SetYTitle("dN/d#Psi_{EP,2}");

  TH1D *histpsi_3_all = new TH1D("histpsi_3_all", expression31,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_3 for all particles
  histpsi_3_all -> Sumw2();
  histpsi_3_all->SetXTitle("#Psi_{EP,3} (radians)");
  histpsi_3_all->SetYTitle("dN/d#Psi_{EP,3}");

  TH1D *histpsi_4_all = new TH1D("histpsi_4_all", expression32,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_4 for all particles
  histpsi_4_all -> Sumw2();
  histpsi_4_all->SetXTitle("#Psi_{EP,4} (radians)");
  histpsi_4_all->SetYTitle("dN/d#Psi_{EP,4}");

  TH1D *histpsi_5_all = new TH1D("histpsi_5_all", expression33,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_5 for all particles
  histpsi_5_all -> Sumw2();
  histpsi_5_all->SetXTitle("#Psi_{EP,5} (radians)");
  histpsi_5_all->SetYTitle("dN/d#Psi_{EP,5}");



  TH1D *histpsi_1_ch = new TH1D("histpsi_1_ch", expression34,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_1 for charged particles
  histpsi_1_ch -> Sumw2();
  histpsi_1_ch ->SetXTitle("#Psi_{EP,1} (radians)");
  histpsi_1_ch ->SetYTitle("dN/d#Psi_{EP,1}");

  TH1D *histpsi_2_ch = new TH1D("histpsi_2_ch", expression35,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_2 for charged particles
  histpsi_2_ch -> Sumw2();
  histpsi_2_ch->SetXTitle("#Psi_{EP,2} (radians)");
  histpsi_2_ch->SetYTitle("dN/d#Psi_{EP,2}");

  TH1D *histpsi_3_ch = new TH1D("histpsi_3_ch", expression36,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_3 for charged particles
  histpsi_3_ch -> Sumw2();
  histpsi_3_ch ->SetXTitle("#Psi_{EP,3} (radians)");
  histpsi_3_ch ->SetYTitle("dN/d#Psi_{EP,3}");

  TH1D *histpsi_4_ch = new TH1D("histpsi_4_ch", expression37,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_4 for charged particles
  histpsi_4_ch -> Sumw2();
  histpsi_4_ch ->SetXTitle("#Psi_{EP,4} (radians)");
  histpsi_4_ch ->SetYTitle("dN/d#Psi_{EP,4}");

  TH1D *histpsi_5_ch = new TH1D("histpsi_5_ch", expression38,100,-(TMath::Pi())/2,(5/2)*TMath::Pi()); //psi_5 for charged particles
  histpsi_5_ch -> Sumw2();
  histpsi_5_ch ->SetXTitle("#Psi_{EP,5} (radians)");
  histpsi_5_ch ->SetYTitle("dN/d#Psi_{EP,5}");


  TH2D *histphi_eta_ch = new TH2D("histphi_eta_ch",expression39,200,0,2*TMath::Pi(),200,-1.1,1.1); //for charged particles
  histphi_eta_ch -> Sumw2();
  histphi_eta_ch->SetXTitle("#phi (radians)");
  histphi_eta_ch->SetYTitle("#eta (pseudo-rapidity)");
  histphi_eta_ch->SetZTitle("dN^{ch.}/d#phid#eta");

  TH2D *histphi_eta_all = new TH2D("histphi_eta_all",expression40,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
  histphi_eta_all -> Sumw2();
  histphi_eta_all->SetXTitle("#phi (radians)");
  histphi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
  histphi_eta_all->SetZTitle("dN^{ch.}/d#phid#eta");



  TH2D *hist_pT_phi_eta_ch = new TH2D("hist_pT_phi_eta_ch",expression41,200,0,2*TMath::Pi(),200,-1.1,1.1); //for charged particles
  hist_pT_phi_eta_ch -> Sumw2();
  hist_pT_phi_eta_ch->SetXTitle("#phi (radians)");
  hist_pT_phi_eta_ch->SetYTitle("#eta (pseudo-rapidity)");
  hist_pT_phi_eta_ch->SetZTitle("dp_{T}^{ch.}/d#phid#eta");

  TH2D *hist_pT_phi_eta_all = new TH2D("hist_pT_phi_eta_all",expression42,200,0,2*TMath::Pi(),200,-1.1,1.1); //for all particles
  hist_pT_phi_eta_all -> Sumw2();
  hist_pT_phi_eta_all->SetXTitle("#phi (radians)");
  hist_pT_phi_eta_all->SetYTitle("#eta (pseudo-rapidity)");
  hist_pT_phi_eta_all->SetZTitle("dp_{T}^{ch.}/d#phid#eta");


//___________________END____BACKGROUND PARTICLE LEVEL HISTOS__________________________________________________//

  TCanvas *c1 = new TCanvas("c1" , "Run Time For 1 TennGen Event vs Event Number (Events Run in Sequence)", 200, 10,500,300); 

  Int_t* time_array;
  Int_t* time_interval_array;
  Int_t* event_array;

  time_array = new Int_t[nEvents + 1];
  time_interval_array = new Int_t[nEvents];
  event_array = new Int_t[nEvents];

  ///////////////////////////////////////////////////////////////////////
  
  //initialize counting variables 
  Int_t N_counter = 0 ; //counts events

  //The following block needs to be used on ACF Only

  TTimeStamp *timestamp = new TTimeStamp();
  UInt_t seed5 = ((jobID+9)*97) + (timestamp->GetSec() - static_cast<Int_t> (1446000000));
  //____________________________________________________________________
  if( (seed5>=0) && (seed5<=900000000) ) {
    //pythia->SetMRPY(1, ((rando->Uniform(1.,1000.))*(N_counter + 5 + seed)) );                   // set seed
    //pythia->SetMRPY(2, 0);                      // use new seed
   // TRandom *Rand = new TRandom(seed); //this I will use for the pT efficiency
    cout<<"Random Seed : "<<seed5<<endl;
  } else {cout << "error: time " << seed5 << " is not valid" << endl; exit(2);}
  

  //________Initiliaze Background Generator________________________________________________//
  
  TennGen *bkgd = new TennGen(); //constructor
  
  bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
  bkgd->SetRandomSeed(seed5); //setting the seed
  bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
  bkgd->SetEtaRange(0.9); //set eta range must be uniform (e.g. |eta| < 0.9, which is ALICE TPC fiducial acceptance)
  bkgd->PrintOutQAHistos(kTRUE);
  bkgd->InitializeBackground();
     //TClonesArray *background =    bkgd->GetBackground();
  TClonesArray *BKGD = bkgd->GetBackground();
  TDatime *time = new TDatime(); //creating new time stamp
  cout<< "\n\nStart Time is " << time->GetTime() <<"\n\n\n" <<endl;
  time_array[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;

  
  //_____________________________________________Starting Event Loop_______________________________________________________________//
  for (Int_t event = 0; event < nEvents; event++) {
    if (event % 10 == 0) cout << "Event # " << event << endl;  
    N_counter ++;


      if( event != 0){
        free(time);
        TDatime *time = new TDatime(); //creating new time stamp
        cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
        time_array[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        //free(bkgd); //reset bkgd generator class variable
        //bkgd->Clear();
        //free(BKGD); //reset BKGD TClonesArray
        BKGD->Clear();
        //TennGen *bkgd = new TennGen();
        //bkgd->SetCentralityBin(cent_bin); //centraility bin 0 is the  0-5 % most central bin
        //Int_t seed_b = (( TMath::CeilNint((rando->Uniform(1.,1000.)))));
        //bkgd->SetRandomSeed( jobID+event ); //setting the seed
        //cout<<"\n\n\n\n\n\njust set new seed\n\n\n\n\n"<<endl;
        //bkgd->SetHarmonics(HF); // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: uniform dN/dphi no harmonics) , (6 : v1 - v2 , v4 - v5) , (7 : v1 - v3 , v5) , (8 : v1 , v3 - v5) , (9 : v1 only) , (10 : v2 only) , (11 : v3 only) , (12 : v4 only) , (13 : v5 only)
        //bkgd->SetEtaRange(0.9); // set eta range (takes care of scaling in the .cxx) MUST BE UNIFORM ETA RANGE e.g. |eta| < 0.9 !!!
        //bkgd->InitializeBackground();
       //TClonesArray *background =    bkgd->GetBackground();
        TClonesArray *BKGD = bkgd->GetBackground_TannerTest();
      } 
    
      bkgd->GetRandomSeed(); 

      Double_t nParticles= BKGD ->GetEntries();
    

      //The following lines are used in the event plane calculations

      Double_t Sum_Psi_1_numer_all = 0;
      Double_t Sum_Psi_1_denom_all = 0;

      Double_t Sum_Psi_2_numer_all = 0;
      Double_t Sum_Psi_2_denom_all = 0;

      Double_t Sum_Psi_3_numer_all = 0;
      Double_t Sum_Psi_3_denom_all = 0;

      Double_t Sum_Psi_4_numer_all = 0;
      Double_t Sum_Psi_4_denom_all = 0;

      Double_t Sum_Psi_5_numer_all = 0;
      Double_t Sum_Psi_5_denom_all = 0;

      Double_t Sum_Psi_1_numer_ch = 0;
      Double_t Sum_Psi_1_denom_ch = 0;

      Double_t Sum_Psi_2_numer_ch = 0;
      Double_t Sum_Psi_2_denom_ch = 0;

      Double_t Sum_Psi_3_numer_ch = 0;
      Double_t Sum_Psi_3_denom_ch = 0;

      Double_t Sum_Psi_4_numer_ch = 0;
      Double_t Sum_Psi_4_denom_ch = 0;

      Double_t Sum_Psi_5_numer_ch = 0;
      Double_t Sum_Psi_5_denom_ch = 0;

    //____________________________Starting Background Particle Loop___________________________________________________//
      for (Int_t N = 0 ; N < nParticles ; N++){
        TMCParticle* bob = (TMCParticle*)BKGD->At(N); 
        Int_t K_F = bob -> GetKF();
        Double_t p = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) + ((bob->GetPz())*(bob->GetPz())) );
        Double_t pT = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) );
        Double_t eta = 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz())));
        Double_t phi = 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx()));

        //using q-vector method to reconstruct the event planes

        Sum_Psi_1_numer_all = Sum_Psi_1_numer_all + pT*TMath::Sin( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_1_denom_all = Sum_Psi_1_denom_all + pT*TMath::Cos( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_2_numer_all = Sum_Psi_2_numer_all + pT*TMath::Sin( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_2_denom_all = Sum_Psi_2_denom_all + pT*TMath::Cos( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_3_numer_all = Sum_Psi_3_numer_all + pT*TMath::Sin( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_3_denom_all = Sum_Psi_3_denom_all + pT*TMath::Cos( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_4_numer_all = Sum_Psi_4_numer_all + pT*TMath::Sin( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_4_denom_all = Sum_Psi_4_denom_all + pT*TMath::Cos( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        Sum_Psi_5_numer_all = Sum_Psi_5_numer_all + pT*TMath::Sin( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        Sum_Psi_5_denom_all = Sum_Psi_5_denom_all + pT*TMath::Cos( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

        if(K_F != 111){ // ANYTHING BUT piPlus

          Sum_Psi_1_numer_ch = Sum_Psi_1_numer_ch + pT*TMath::Sin( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_1_denom_ch = Sum_Psi_1_denom_ch + pT*TMath::Cos( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_2_numer_ch = Sum_Psi_2_numer_ch + pT*TMath::Sin( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_2_denom_ch = Sum_Psi_2_denom_ch + pT*TMath::Cos( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_3_numer_ch = Sum_Psi_3_numer_ch + pT*TMath::Sin( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_3_denom_ch = Sum_Psi_3_denom_ch + pT*TMath::Cos( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
 
          Sum_Psi_4_numer_ch = Sum_Psi_4_numer_ch + pT*TMath::Sin( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_4_denom_ch = Sum_Psi_4_denom_ch + pT*TMath::Cos( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_5_numer_ch = Sum_Psi_5_numer_ch + pT*TMath::Sin( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_5_denom_ch = Sum_Psi_5_denom_ch + pT*TMath::Cos( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
        }

        histpT_all->Fill(pT); //pT
        histphi_all->Fill( phi ); //azimuthal angle
        histeta_all->Fill(eta); //eta

        histphi_eta_all->Fill(phi,eta);
        hist_pT_phi_eta_all->Fill(phi,eta,pT);
      
        //cout<<"Particle "<<N<<" KF = "<< K_F<<endl;
        if(K_F == 211){ //piPlus
          histpT_piPlus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) ); //transverse momentum
          histeta_piPlus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) ); //pseudo-rapidity
          histphi_piPlus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ); //azimuthal angle 
          histpT_all->Fill(pT); //pT

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        }

        else if(K_F == -211){//piMinus
          histpT_piMinus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_piMinus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_piMinus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );     

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT); 
        } 

        else if(K_F == 111){//piZero
          histpT_piZero->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_piZero->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_piZero->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );
        }

        else if(K_F == 321){//kPlus
          histpT_kPlus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_kPlus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_kPlus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ); 

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        } 

        else if(K_F == -321){//kMinus
          histpT_kMinus->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_kMinus->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_kMinus->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );   

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);   
        }

        else if(K_F == 2212){//p
          histpT_p->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_p->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_p->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );    

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        } 

        else if(K_F == -2212){//pbar
          histpT_pbar->Fill( TMath::Sqrt(bob->GetPx() * bob->GetPx() + bob->GetPy() * bob->GetPy()) );  
          histeta_pbar->Fill( 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz()))) );
          histphi_pbar->Fill( 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) );

          histphi_ch->Fill( phi );
          histpT_ch->Fill(pT);
          histeta_ch->Fill(eta);
          histphi_eta_ch->Fill(phi,eta);
          hist_pT_phi_eta_ch->Fill(phi,eta,pT);
        }

     
      }//END PARTICLE LOOP BACKGROUND

      histpsi_1_all->Fill( (1.0/1.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_1_numer_all,-Sum_Psi_1_denom_all)));
      histpsi_2_all->Fill( (1.0/2.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_2_numer_all,-Sum_Psi_2_denom_all)));
      histpsi_3_all->Fill( (1.0/3.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_3_numer_all,-Sum_Psi_3_denom_all)));
      histpsi_4_all->Fill( (1.0/4.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_4_numer_all,-Sum_Psi_4_denom_all)));
      histpsi_5_all->Fill( (1.0/5.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_5_numer_all,-Sum_Psi_5_denom_all)));


      histpsi_1_ch->Fill( (1.0/1.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_1_numer_ch,-Sum_Psi_1_denom_ch)));
      histpsi_2_ch->Fill( (1.0/2.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_2_numer_ch,-Sum_Psi_2_denom_ch)));
      histpsi_3_ch->Fill( (1.0/3.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_3_numer_ch,-Sum_Psi_3_denom_ch)));
      histpsi_4_ch->Fill( (1.0/4.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_4_numer_ch,-Sum_Psi_4_denom_ch)));
      histpsi_5_ch->Fill( (1.0/5.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_5_numer_ch,-Sum_Psi_5_denom_ch)));
      

  }//________________________________________________END EVENT LOOP________________________________________

  bkgd->CloserFunction(); //close all the damn files, and print out QA histos to file if desired

  TFile *ff = new TFile(expression1,"RECREATE"); //make output files and directories

  TDirectory *particle_level_histos = ff->mkdir("particle_level_histos");
  TDirectory *all_charged_particles = ff->mkdir("all_charged_particles");
  TDirectory *all_particles = ff->mkdir("all_particles");

   TDatime *time2 = new TDatime(); //creating new time stamp
   cout << "\n\n\n\n\nTime at Event # "<< nEvents << " is " << time2->GetTime() <<"\n\n\n\n" << endl;
   time_array[nEvents] = ((time2->GetYear()*31540000) + (time2->GetMonth()*2628000) + (time2->GetDay()*86400) + (time2->GetHour()*3600) + (time2->GetMinute()*60) + time2->GetSecond()) ;

   //calculating the average time for one event

   Int_t Sum_Time = 0;
   Int_t Total_Time = 0;

   for( Int_t i = 0 ; i < nEvents; i++){
     time_interval_array[i] = time_array[i + 1] - time_array[i];
     Sum_Time = Sum_Time + time_interval_array[i];
     event_array[i] = i + 1;
   }

   cout<< "\n\n\n\n Total Time for = "<<nEvents<<" Events = "<<Sum_Time<<" Seconds \n\n\n\n"<<endl;

   cout<< "\n\n\n\n Average Time For One Event = "<<Sum_Time/(nEvents)<<" Seconds \n\n\n\n"<<endl;

  TH1D *histruntime = new TH1D("histruntime", "Run Times for Background Generator",8,((Sum_Time/nEvents) - 1),((Sum_Time/nEvents) + 1)); //for non-random distribution
  histruntime -> Sumw2();
  histruntime->SetXTitle("Run Time For 1 Event (Seconds)");
  histruntime->SetYTitle("Counts");

  TGraph *rtgraph = new TGraph(nEvents , event_array , time_interval_array);
  rtgraph->SetTitle("Run Time For 1 Background Event vs Event Number (Events Run in Sequence)");
  rtgraph->GetXaxis()->SetTitle("Event Iteration Number");
  rtgraph->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
  rtgraph->Write();

  c1->cd(1); 
  rtgraph->Draw("AC");
  c1->Draw();
  c1->Write("");
  //c1->Close();

  for( Int_t j = 0 ; j < nEvents; j++){
    histruntime->Fill(time_interval_array[j]);
  }

  histruntime->Write();
   //________________________________Setting Histogram Styles and Writing to output file_______________________________________//

   histpT_piPlus->SetMarkerColor(2);
   histpT_piPlus->SetMarkerStyle(20);
   histpT_piMinus->SetMarkerColor(3);
   histpT_piMinus->SetMarkerStyle(21);
   histpT_kPlus->SetMarkerColor(4);
   histpT_kPlus->SetMarkerStyle(22);
   histpT_kMinus->SetMarkerColor(5);
   histpT_kMinus->SetMarkerStyle(23);
   histpT_p->SetMarkerColor(6);
   histpT_p->SetMarkerStyle(43);
   histpT_pbar->SetMarkerColor(7);
   histpT_pbar->SetMarkerStyle(47);
   histpT_piZero->SetMarkerColor(13);
   histpT_piZero->SetMarkerStyle(34);

   histeta_piPlus->SetMarkerColor(2);
   histeta_piPlus->SetMarkerStyle(20);
   histeta_piMinus->SetMarkerColor(3);
   histeta_piMinus->SetMarkerStyle(21);
   histeta_kPlus->SetMarkerColor(4);
   histeta_kPlus->SetMarkerStyle(22);
   histeta_kMinus->SetMarkerColor(5);
   histeta_kMinus->SetMarkerStyle(23);
   histeta_p->SetMarkerColor(6);
   histeta_p->SetMarkerStyle(43);
   histeta_pbar->SetMarkerColor(7);
   histeta_pbar->SetMarkerStyle(47);
   histeta_piZero->SetMarkerColor(13);
   histeta_piZero->SetMarkerStyle(34);

   histphi_piPlus->SetMarkerColor(2);
   histphi_piPlus->SetMarkerStyle(20);
   histphi_piMinus->SetMarkerColor(3);
   histphi_piMinus->SetMarkerStyle(21);
   histphi_kPlus->SetMarkerColor(4);
   histphi_kPlus->SetMarkerStyle(22);
   histphi_kMinus->SetMarkerColor(5);
   histphi_kMinus->SetMarkerStyle(23);
   histphi_p->SetMarkerColor(6);
   histphi_p->SetMarkerStyle(43);
   histphi_pbar->SetMarkerColor(7);
   histphi_pbar->SetMarkerStyle(47);
   histphi_piZero->SetMarkerColor(13);
   histphi_piZero->SetMarkerStyle(34);


   particle_level_histos->cd();

   histpT_piPlus->Write();
   histpT_piMinus->Write();
   histpT_kPlus->Write();
   histpT_kMinus->Write();
   histpT_p->Write();
   histpT_pbar->Write();
   histpT_piZero->Write();

   histeta_piPlus->Write();
   histeta_piMinus->Write();
   histeta_kPlus->Write();
   histeta_kMinus->Write();
   histeta_p->Write();
   histeta_pbar->Write();
   histeta_piZero->Write();

   histphi_piPlus->Write();
   histphi_piMinus->Write();
   histphi_kPlus->Write();
   histphi_kMinus->Write();
   histphi_p->Write();
   histphi_pbar->Write();
   histphi_piZero->Write();

   all_charged_particles->cd();

   histpT_ch->Write();
   histeta_ch->Write();
   histphi_ch->Write();

   histpsi_1_ch->Write();
   histpsi_2_ch->Write();
   histpsi_3_ch->Write();
   histpsi_4_ch->Write();
   histpsi_5_ch->Write();

   histphi_eta_ch->Write();
   hist_pT_phi_eta_ch->Write();

   all_particles->cd();

   histpT_all->Write();
   histeta_all->Write();
   histphi_all->Write();
   
   histpsi_1_all->Write();
   histpsi_2_all->Write();
   histpsi_3_all->Write();
   histpsi_4_all->Write();
   histpsi_5_all->Write();

   histphi_eta_all->Write();
   hist_pT_phi_eta_all->Write();

 

//__________________________________________________________________________________________________________________________________//
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe <pT> for pi_plus = "<<histpT_piPlus->GetMean(1)<<" GeV/c, The Sigma for pi_plus = "<<histpT_piPlus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for pi_minus = "<<histpT_piMinus->GetMean(1)<<" GeV/c, The Sigma for pi_minus = "<<histpT_piMinus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for pi_zero = "<<histpT_piZero->GetMean(1)<<" GeV/c, The Sigma for pi_zero = "<<histpT_piZero->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for k_plus = "<<histpT_kPlus->GetMean(1)<<" GeV/c, The Sigma for k_plus = "<<histpT_kPlus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for k_minus = "<<histpT_kMinus->GetMean(1)<<" GeV/c, The Sigma for k_minus = "<<histpT_kMinus->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for p = "<<histpT_p->GetMean(1)<<" GeV/c, The Sigma for p = "<<histpT_p->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for p_bar = "<<histpT_pbar->GetMean(1)<<" GeV/c, The Sigma for p_bar = "<<histpT_pbar->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe <pT> for all = "<<histpT_all->GetMean(1)<<" GeV/c, The Sigma for all = "<<histpT_all->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"\nThe <pT> for ch = "<<histpT_ch->GetMean(1)<<" GeV/c, The Sigma for ch = "<<histpT_ch->GetStdDev(1)<<" GeV/c"<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;
   cout<<"\nThe Nall/event = "<<histphi_all->GetEntries()/nEvents<<" particles/event "<<endl;
   cout<<"\nThe Nch/event  = "<<histphi_ch->GetEntries()/nEvents<<" charged particles/event "<<endl;
   cout<<"//____________________________________________________________________________________________________________________________//\n"<<endl;	
//__________________________________________________________________________________________________________________________________//

 

   
  ff -> Write();
  ff -> Close();
  bkgd->~TennGen();
} //end Background_Test_Macro.C

void VonNeumannTest(Int_t nEvents, Int_t jobID , Int_t HF , Int_t cent_bin, TString OutputDir ){
  //check output directory exists
  if (gSystem->AccessPathName(OutputDir)) {
    cout << "Output directory does not exist. Creating directory " << OutputDir << endl;
    gSystem->mkdir(OutputDir);
  }

  cout <<"Starting old tests" << endl;
 // TennGen_Old_Test(nEvents,jobID ,HF ,cent_bin,OutputDir);
  cout << "Starting new tests" << endl;
  TennGen_New_Test(nEvents,jobID ,HF ,cent_bin,OutputDir);
  cout << "Starting background tests" << endl;
  if(HF==0 || cent_bin ==0) Time_Compare_Tests(nEvents,jobID ,HF ,cent_bin,OutputDir);
} //end Background_Test_Macro.C

