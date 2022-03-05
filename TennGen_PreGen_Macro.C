
// Copyright (c) 2019; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// MODIFIED BY CHARLES HUGHES <chughe26@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2022-03-04 19:31 EST
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
#include "BkgrLoad.h"
#include "TCanvas.h"
#include "TTimeStamp.h"
#endif
		

class TennGen;

//#include "TPythia8.h"
//class TPythia8;

// nEvents is how many events we want. This is dynamic user input to be run from the command prompt and user specified
void TennGen_PreGen_Macro(Int_t nEvents, Int_t jobID , Int_t HF, Int_t cent_bin, Int_t Num_BKGD_Files, Bool_t RT_Stats = kFALSE , Bool_t GRID = kFALSE){


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

  sprintf(expression1 , "TennGen_PreGen_Test_File_HF=%d_cent_bin=%d.root" , HF , cent_bin);

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

  TDatime *time = new TDatime(); //creating new time stamp
  cout<< "\n\nStart Time is " << time->GetTime() <<"\n\n\n" <<endl;

  time_array[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;

  char run_title_str[256];
  sprintf( run_title_str, "BG_Gen_Cent%d_Harm%d" , cent_bin , HF);
 
  Float_t px; 
  Float_t py;
  Float_t pz;
  //Float_t m;
  //Float_t E;
  Short_t KF;
	  
  UInt_t centrality = cent_bin;
  UInt_t harmonics = HF;
  UInt_t version = 0;
  UInt_t events = nEvents;
                                                       
  BkgrLoad *bkgd2 = new BkgrLoad();

  UInt_t seed6 = ((jobID+834)*33) + (timestamp->GetSec() - static_cast<Int_t> (1446000000)); // for file shuffling

//____________________________________________________________________________________________________________________//

  char filepathstr[512];

  if( !GRID ){

    //PUT THE PATH TO THE BKGD_ROOT_FILES directory here !!!!!!!//
    //sprintf(filepathstr,"/home/alidock/ML/BKGD_ROOT_FILES");
    sprintf(filepathstr,"/home/charles/Documents/research/Background_Research/forcharles/Newest_Background_Code_05_15_2018/Harmonic_Code_for_copying/DiJet_Asymmetry/Updated_Code/Latest_most_up_to_date_Background_Code_08_31_2018/Latest_Version_meant_for_anti_kT/Frag_Func_Code/Heavy_Ion_BGLoad/Patrick_Studies/ML-Jet-BG-Subtraction/BKGD_ROOT_FILES");

  }
  else if( GRID ){

   getcwd(filepathstr, sizeof(filepathstr));

  }

  bkgd2->PassInSettings(filepathstr,  0 , 0 , 0 , seed6 );
  bkgd2->Load();
  bkgd2->PrintStatus();
  TClonesArray *BKGD = (TClonesArray *)bkgd2->GetEventBG(); //(outside particle loop)

  
  //_____________________________________________Starting Event Loop_______________________________________________________________//
  for (Int_t event = 0; event < nEvents; event++) {
    //if( (bkgd2->GetActiveFile()) == (bkgd2->GetNumFiles()) ) //GetActiveFile() starts at 0 

    if( (bkgd2->GetActiveFile()) == Num_BKGD_Files ){ //GetActiveFile() starts at 0  
      break; //break out of the loop if you have used up all the background events
    }
    if (event % 10 == 0) cout << "========================Event # " << event << endl;  



      if( event != 0){
        free(time);
        TDatime *time = new TDatime(); //creating new time stamp
        cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
        time_array[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        delete BKGD;
        TClonesArray *BKGD = (TClonesArray *)bkgd2->GetEventBG(); //(outside particle loop)
        cout<<"you got the new background event\n\n"<<endl;
      } 

      Double_t mass;
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

        if (K_F == 211 || K_F == -211) mass = 0.13957; //pi+ & pi- (GeV/c^2) (inside particle loop)
	else if (K_F  == 321 || K_F == -321) mass = 0.493677; //k+ && k- (GeV/c^2)
	else if (K_F == 2212 || K_F == -2212) mass = 0.938272; //p & pbar (GeV/c^2)
	else if (K_F == 111) mass = 0.134977; //pi0 (GeV/c^2)

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

  bkgd2->CloseFiles(); //close all the damn files, and print out QA histos to file if desired

  TFile *ff = new TFile(expression1,"RECREATE"); //make output files and directories

  TDirectory *particle_level_histos = ff->mkdir("particle_level_histos");
  TDirectory *all_charged_particles = ff->mkdir("all_charged_particles");
  TDirectory *all_particles = ff->mkdir("all_particles");


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
} //end Background_Test_Macro.C

