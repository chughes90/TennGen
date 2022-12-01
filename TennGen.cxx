//// Author: Charles Philip Hughes <chughe26@vols.utk.edu>
// Update: 2020-05-04 18:55:27-0500
// Copyright (c) 2017; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// University of Tennessee ALICE collaboration
      
#include "TMath.h"   
#include "TRandom.h"   
#include "TennGen.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TColor.h"
#include "Rtypes.h"
#include <fstream>
#include <iostream>
#include "TLatex.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TMCParticle.h"
#include "TRandom3.h"


ClassImp(TennGen)
//____________________________________________________________________________________________________________________//
//___________________________Gets the Integral in the BlastWave_____________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::MyIntegrandBG(Double_t *x, Double_t *p){
  // integrand for boltzman-gibbs blast wave
     // x[0] -> r (radius)
     // p[0] -> mass
     // p[1] -> pT (transverse momentum)
     // p[2] -> beta_max (surface velocity)
     // p[3] -> T (freezout temperature)
     // p[4] -> n (velocity profile)


  Double_t x0 = x[0]; 
  
  Double_t mass     = p[0];
  Double_t pT       = p[1];
  Double_t beta_max = p[2];
  Double_t temp     = p[3];
  Double_t n      = p[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  Double_t mT      = TMath::Sqrt(mass*mass+pT*pT);

  Double_t rho0   = TMath::ATanH(beta);  
  Double_t arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  Double_t arg01 = mT*TMath::CosH(rho0)/temp;
  Double_t f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);

  return f0;
}

//____________________________________________________________________________________________________________________//
//_________________________________________Get MyStaticBGdNdPtTimesPt________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::MyStaticBGdNdPtTimesPt(Double_t *x, Double_t *p) {
  // implementation of BGBW (1/pt dNdpt)

  Double_t pT = x[0];;
  

  Double_t mass    = p[0];
  Double_t beta    = p[1];
  Double_t temp    = p[2];
  Double_t n       = p[3];
  Double_t norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", MyIntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  Double_t result = fIntBG->Integral(0,1);
  //  printf ("[%4.4f], Int :%f\n", pT, result);
  return result*norm*pT;//*1e30;;
}

//____________________________________________________________________________________________________________________//
//_________________________________________Get BlastWavedNdptTimesPt________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::BlastWavedNdptTimesPt(Double_t *x, Double_t *p) {
  Double_t pT = x[0];
  
  Double_t mass    = p[0];
  Double_t A       = p[5];

  return MyStaticBGdNdPtTimesPt(x,p);
}
//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::fitf(Double_t *x, Double_t *par) 
{
   Double_t arg = x[0] /par[1];
   Double_t fitval = par[0]*TMath::Exp(-arg);
   //Double_t fitval = 10 * (par[1]/par[1]) * (par[0]/par[0]) * (x[0]/x[0]);
   return fitval;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi0(Double_t *phi, Double_t *par  ) //v1 - v5
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi1(Double_t *phi, Double_t *par  ) //v2 - v5
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi2(Double_t *phi, Double_t *par  )  //v3 - v5
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi3(Double_t *phi, Double_t *par  )  //v1 - v4 (exclude v5)
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi4(Double_t *phi, Double_t *par  ) //v1 - v3 
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi6(Double_t *phi, Double_t *par  )  //v1 -v2, v4 - v5 (exclude v3)
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 0.0;
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi7(Double_t *phi, Double_t *par  )  //v1 - v3, v5 (exclude v4)
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//


Double_t TennGen::dNdphi8(Double_t *phi, Double_t *par  )  //v1, v3 - v5 (exclude v2)
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi9(Double_t *phi, Double_t *par  )  //v1 only
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 2.0*par[0]*TMath::Cos( par[5] - phi[0] );
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 0.0;
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi10(Double_t *phi, Double_t *par  )  //v2 only
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 2.0*par[1]*TMath::Cos( 2.0*( 0.0 - phi[0] ) );
   Double_t Third_Term = 0.0;
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi11(Double_t *phi, Double_t *par  )  //v3 only
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 2.0*par[2]*TMath::Cos( 3.0*( par[6] - phi[0]) );
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi12(Double_t *phi, Double_t *par  )  //v4 only
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 0.0;
   Double_t Fourth_Term = 2.0*par[3]*TMath::Cos( 4.0*( 0.0 - phi[0] ) );
   Double_t Fifth_Term = 0.0;
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

Double_t TennGen::dNdphi13(Double_t *phi, Double_t *par  )  //v5 only
{

   Double_t Zero_Term = 1.0; 
   Double_t First_Term = 0.0;
   Double_t Second_Term = 0.0;
   Double_t Third_Term = 0.0;
   Double_t Fourth_Term = 0.0;
   Double_t Fifth_Term = 2.0*par[4]*TMath::Cos( 5.0*( par[7] - phi[0]) );
   
   Double_t Phi_Dist_Harmonics = Zero_Term + First_Term + Second_Term + Third_Term + Fourth_Term + Fifth_Term;
   
   return Phi_Dist_Harmonics;
}

//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//

void TennGen::Set_Phi_Dist(Double_t pT , Double_t KF , Double_t Psi_1_event , Double_t Psi_3_event, Double_t Psi_5_event) {

  if (Harmonics_Phi_Dist!=0) delete Harmonics_Phi_Dist;

  if( HarmonicsFlag == 0){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi0,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }
  
  else if( HarmonicsFlag == 1){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi1,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }
  
  else if( HarmonicsFlag == 2){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi2,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }
  
  else if( HarmonicsFlag == 3){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi3,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }
  
  else if( HarmonicsFlag == 4){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi4,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 6){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi6,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 7){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi7,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 8){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi8,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 9){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi9,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 10){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi10,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 11){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi11,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 12){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi12,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

  else if( HarmonicsFlag == 13){
    Harmonics_Phi_Dist = new TF1("Harmonics_Phi_Dist",dNdphi13,0,2.0*TMath::Pi(),8); 
    Harmonics_Phi_Dist->SetNpx(1000);
  }

/*_________________________________________________________________________________*/
   
    if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v2_pi_eval = 0;
    }

    else if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v2_pi_eval = (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v3_pi_eval = 0;
    }

    else if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v3_pi_eval = (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v4_pi_eval = 0;
    }

    else if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v4_pi_eval = (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v5_pi_eval = 0;
    }

    else if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v5_pi_eval = (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

/*_________________________________________________________________________________*/

    if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v2_K_eval = 0;
    }

    else if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v2_K_eval = (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v3_K_eval = 0;
    }

    else if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v3_K_eval = (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v4_K_eval = 0;
    }

    else if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v4_K_eval = (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v5_K_eval = 0;
    }

    else if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v5_K_eval = (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

/*_________________________________________________________________________________*/

    if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v2_P_eval = 0;
    }

    else if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v2_P_eval = (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v3_P_eval = 0;
    }

    else if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v3_P_eval = (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v4_P_eval = 0;
    }

    else if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v4_P_eval = (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

    if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
      v5_P_eval = 0;
    }

    else if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
      v5_P_eval = (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
    }

/*_________________________________________________________________________________*/ 
  
    if( KF == 111 || KF == 211 || KF == -211){
  
    Harmonics_Phi_Dist->SetParameter(0, v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
    Harmonics_Phi_Dist->SetParameter(1, v2_pi_eval                           );
    Harmonics_Phi_Dist->SetParameter(2, v3_pi_eval                           );
    Harmonics_Phi_Dist->SetParameter(3, v4_pi_eval                           );
    Harmonics_Phi_Dist->SetParameter(4, v5_pi_eval                           );
    Harmonics_Phi_Dist->SetParameter(5, Psi_1_event                          );
    Harmonics_Phi_Dist->SetParameter(6, Psi_3_event                          );
    Harmonics_Phi_Dist->SetParameter(7, Psi_5_event                          );

    }
  
    else if( KF == 321 || KF == -321){
  
    Harmonics_Phi_Dist->SetParameter(0, v1_K->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
    Harmonics_Phi_Dist->SetParameter(1, v2_K_eval                           );
    Harmonics_Phi_Dist->SetParameter(2, v3_K_eval                           );
    Harmonics_Phi_Dist->SetParameter(3, v4_K_eval                           );
    Harmonics_Phi_Dist->SetParameter(4, v5_K_eval                           );
    Harmonics_Phi_Dist->SetParameter(5, Psi_1_event                         );
    Harmonics_Phi_Dist->SetParameter(6, Psi_3_event                         );
    Harmonics_Phi_Dist->SetParameter(7, Psi_5_event                         );
  
    }
  
    else if( KF == 2212 || KF == -2212){
  
    Harmonics_Phi_Dist->SetParameter(0, v1_P->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
    Harmonics_Phi_Dist->SetParameter(1, v2_P_eval                           );
    Harmonics_Phi_Dist->SetParameter(2, v3_P_eval                           );
    Harmonics_Phi_Dist->SetParameter(3, v4_P_eval                           );
    Harmonics_Phi_Dist->SetParameter(4, v5_P_eval                           );
    Harmonics_Phi_Dist->SetParameter(5, Psi_1_event                         );
    Harmonics_Phi_Dist->SetParameter(6, Psi_3_event                         );
    Harmonics_Phi_Dist->SetParameter(7, Psi_5_event                         );
  
    }
    
  } 

//______________________________________________________________________________
TennGen::TennGen()//constructor
{

  //CentralityBin = 0;
  BkgdArray = new TClonesArray("TMCParticle" , 1000);
  cout<<"about to make TF1's"<<endl;
  pTdist_piPlus = new TF1("pTdist_piPlus",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_piPlus->SetNpx (1000);
  pTdist_piMinus = new TF1("pTdist_piMinus",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_piMinus->SetNpx (1000);  
  pTdist_kPlus = new TF1("pTdist_kPlus",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_kPlus->SetNpx (1000);  
  pTdist_kMinus = new TF1("pTdist_kMinus",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_kMinus->SetNpx (1000);  
  pTdist_p = new TF1("pTdist_p",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_p->SetNpx (1000);  
  pTdist_pbar = new TF1("pTdist_pbar",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_pbar->SetNpx (1000);  
  pTdist_piZero = new TF1("pTdist_piZero",BlastWavedNdptTimesPt ,0,100,5); 
  pTdist_piZero->SetNpx (1000);  
  
  //The next lines are for the vn functions
  //vn pions (piZero , piPlus, piMinus )
  v1_pi = new TF1("v1_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
  //v1_pi = new TF1("v1_pi", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + 0.0", 0, 100 );
  v1_pi->SetNpx(1000);
  v2_pi = new TF1("v2_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) ) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v2_pi = new TF1("v2_pi", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) ) + 0.0", 0, 100 );
  v2_pi->SetNpx(1000);
  v3_pi = new TF1("v3_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v3_pi = new TF1("v3_pi", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v3_pi->SetNpx(1000);
  v4_pi = new TF1("v4_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v4_pi = new TF1("v4_pi", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v4_pi->SetNpx(1000);
  v5_pi = new TF1("v5_pi", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v5_pi = new TF1("v5_pi", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.5", 0, 100 );
  v5_pi->SetNpx(1000);
  //vn kaons (K+ , K-)
  v1_K = new TF1("v1_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
  //v1_K = new TF1("v1_K", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + 0.0", 0, 100 );
  v1_K->SetNpx(1000);
  v2_K = new TF1("v2_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v2_K = new TF1("v2_K", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) ) + 0.0", 0, 100 );
  v2_K->SetNpx(1000);
  v3_K = new TF1("v3_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v3_K = new TF1("v3_K", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v3_K->SetNpx(1000);
  v4_K = new TF1("v4_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v4_K = new TF1("v4_K", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v4_K->SetNpx(1000);
  v5_K = new TF1("v5_K", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v5_K = new TF1("v5_K", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.5", 0, 100 );
  v5_K->SetNpx(1000);
  //vn protons (P , Pbar)
  v1_P = new TF1("v1_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4) - 0.02)", 0, 100 );
  //v1_P = new TF1("v1_P", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) - 0.02) + 0.0", 0, 100 );
  v1_P->SetNpx(1000);
  v2_P = new TF1("v2_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v2_P = new TF1("v2_P", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4) ) + 0.0", 0, 100 );
  v2_P->SetNpx(1000);
  v3_P = new TF1("v3_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v3_P = new TF1("v3_P", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v3_P->SetNpx(1000);
  v4_P = new TF1("v4_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v4_P = new TF1("v4_P", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.0", 0, 100 );
  v4_P->SetNpx(1000);
  v5_P = new TF1("v5_P", "(x<[5])*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + (x>[5])*([0] + ([1]*[5]) + ([2]*[5]^2) + ([3]*[5]^3) + ([4]*[5]^4))", 0, 100 );
  //v5_P = new TF1("v5_P", "(0)*([0] + ([1]*x) + ([2]*x^2) + ([3]*x^3) + ([4]*x^4)) + 0.5", 0, 100 );
  v5_P->SetNpx(1000);
  

  //Harmonics_Phi_Dist_Rand = 0;

  Harmonics_Phi_Dist = 0;

  New_Rand = 0;

  //TDirectory *vn_QA_plots = vn_QA_file->mkdir("vn_QA_plots");
  //TDirectory *Event_Plane_QA_plots = vn_QA_file->mkdir("Event_Plane_QA_plots");

  v1_pi_h->SetXTitle("p_{T} (GeV/c)");
  v1_pi_h->SetYTitle("v_{1}^{#pi}");
  v2_pi_h->SetXTitle("p_{T} (GeV/c)");
  v2_pi_h->SetYTitle("v_{2}^{#pi}");
  v3_pi_h->SetXTitle("p_{T} (GeV/c)");
  v3_pi_h->SetYTitle("v_{3}^{#pi}");
  v4_pi_h->SetXTitle("p_{T} (GeV/c)");
  v4_pi_h->SetYTitle("v_{4}^{#pi}");
  v5_pi_h->SetXTitle("p_{T} (GeV/c)");
  v5_pi_h->SetYTitle("v_{5}^{#pi}");
  v1_K_h->SetXTitle("p_{T} (GeV/c)");
  v1_K_h->SetYTitle("v_{1}^{K}");
  v2_K_h->SetXTitle("p_{T} (GeV/c)");
  v2_K_h->SetYTitle("v_{2}^{K}");
  v3_K_h->SetXTitle("p_{T} (GeV/c)");
  v3_K_h->SetYTitle("v_{3}^{K}");
  v4_K_h->SetXTitle("p_{T} (GeV/c)");
  v4_K_h->SetYTitle("v_{4}^{K}");
  v5_K_h->SetXTitle("p_{T} (GeV/c)");
  v5_K_h->SetYTitle("v_{5}^{K}");
  v1_P_h->SetXTitle("p_{T} (GeV/c)");
  v1_P_h->SetYTitle("v_{1}^{P}");
  v2_P_h->SetXTitle("p_{T} (GeV/c)");
  v2_P_h->SetYTitle("v_{2}^{P}");
  v3_P_h->SetXTitle("p_{T} (GeV/c)");
  v3_P_h->SetYTitle("v_{3}^{P}");
  v4_P_h->SetXTitle("p_{T} (GeV/c)");
  v4_P_h->SetYTitle("v_{4}^{P}");
  v5_P_h->SetXTitle("p_{T} (GeV/c)");
  v5_P_h->SetYTitle("v_{5}^{P}");

  v1_sq_h->SetXTitle("v_{1}^{2}");
  v1_sq_h->SetYTitle("counts");
  v1_sq_h->Sumw2();

  v2_sq_h->SetXTitle("v_{2}^{2}");
  v2_sq_h->SetYTitle("counts");
  v2_sq_h->Sumw2();

  v3_sq_h->SetXTitle("v_{3}^{2}");
  v3_sq_h->SetYTitle("counts");
  v3_sq_h->Sumw2();

  v4_sq_h->SetXTitle("v_{4}^{2}");
  v4_sq_h->SetYTitle("counts");
  v4_sq_h->Sumw2();

  v5_sq_h->SetXTitle("v_{5}^{2}");
  v5_sq_h->SetYTitle("counts");
  v5_sq_h->Sumw2();


  Psi_1_h->SetXTitle("#Psi_{1} (radians)");
  Psi_1_h->SetYTitle("dN/d#Psi_{1}");
  Psi_1_h->Sumw2();
  Psi_2_h->SetXTitle("#Psi_{2} (radians)");
  Psi_2_h->SetYTitle("dN/d#Psi_{2}");
  Psi_2_h->Sumw2();
  Psi_3_h->SetXTitle("#Psi_{3} (radians)");
  Psi_3_h->SetYTitle("dN/d#Psi_{3}");
  Psi_3_h->Sumw2();
  Psi_4_h->SetXTitle("#Psi_{4} (radians)");
  Psi_4_h->SetYTitle("dN/d#Psi_{4}");
  Psi_4_h->Sumw2();
  Psi_5_h->SetXTitle("#Psi_{5} (radians)");
  Psi_5_h->SetYTitle("dN/d#Psi_{5}");
  Psi_5_h->Sumw2();

  
  //Test kMinus
  //
  // mass = 0.139570
  //
  // beta = 0.937493
  //
  // temp = 0.0838151
  //
  // n = 0.151064
  //
  // norm = 37751.6
  //


}



//______________________________________________________________________________
TennGen::~TennGen()//destructor
{
  cout<<"You stepped into the destructor !!!"<<endl;
  delete BkgdArray;
  delete pTdist_piPlus;
  delete pTdist_piMinus;
  delete pTdist_kPlus;
  delete pTdist_kMinus;
  delete pTdist_p;
  delete pTdist_pbar;
  delete pTdist_piZero;
  delete etadist;
  //delete phidist;
  delete Psi_1;
  delete Psi_3;
  delete Psi_5;
  delete New_Rand;
  delete v1_pi;
  delete v2_pi;
  delete v3_pi;
  delete v4_pi;
  delete v5_pi;
  delete v1_K;
  delete v2_K;
  delete v3_K;
  delete v4_K;
  delete v5_K;
  delete v1_P;
  delete v2_P;
  delete v3_P;
  delete v4_P;
  delete v5_P;
  delete Harmonics_Phi_Dist;
  delete Harmonics_Phi_Dist_Rand;

/*
  //delete vn_QA_file;
  //delete vn_QA_plots;
  //delete Event_Plane_QA_plots;

  delete v1_pi_h;
  delete v2_pi_h;
  delete v3_pi_h;
  delete v4_pi_h;
  delete v5_pi_h;
  delete v1_K_h;
  delete v2_K_h;
  delete v3_K_h;
  delete v4_K_h;
  delete v5_K_h;
  delete v1_P_h;
  delete v2_P_h;
  delete v3_P_h;
  delete v4_P_h;
  delete v5_P_h;
  delete Psi_1_h;
  delete Psi_2_h;
  delete Psi_3_h;
  delete Psi_4_h;
  delete Psi_5_h;
*/
}

void TennGen::InitializeBackground(){
  cout<<"Now I am initializing background for centrality bin "<<CentralityBin<<endl;
  MyPrivateFunction(); //function for debugging purposes
  SetpiPlusParams();
  SetpiMinusParams();
  SetkPlusParams();
  SetkMinusParams();
  SetpParams();
  SetpbarParams();
  SetpiZeroParams();
  Setv1_Params();
  Setv2_Params();
  Setv3_Params();
  Setv4_Params();
  Setv5_Params();
  //etadist = new TRandom3();
  //etadist = 0;
  UInt_t seed_eta = Seed * 7;
  etadist = new TRandom3(seed_eta);
  //phidist = new TRandom3(974);
  
  //TRandom3 for 1st, 3rd, and 5th order event planes
  UInt_t seed_psi_1 = (Seed*2)+3;
  Psi_1 = new TRandom3(seed_psi_1);
  //Psi_1->Clear();
  //Psi_1->SetSeed((Seed*2.)+3);
  cout << "Random Psi_1 Seed from .cxx = "<< (Psi_1->GetSeed()) <<" for Event " << endl;
  //Psi_1 = 0;
  UInt_t seed_psi_3 = (Seed*3)+2;
  Psi_3 = new TRandom3(seed_psi_3);
  //Psi_3->Clear();
  //Psi_3->SetSeed((Seed*3.)+2);
  cout << "Random Psi_3 Seed from .cxx = "<< (Psi_3->GetSeed()) <<" for Event " << endl;
  //Psi_3 = 0;
  UInt_t seed_psi_5 = (Seed*5)*9;
  Psi_5 = new TRandom3(seed_psi_5);
  //Psi_5->Clear();
  //Psi_5->SetSeed((Seed*5.)*9);
  cout << "Random Psi_5 Seed from .cxx = "<< (Psi_5->GetSeed()) <<" for Event " << endl;
  //Psi_5 = 0;

  UInt_t seed_uniform_random_phi = (Seed*3)*4;
  Harmonics_Phi_Dist_Rand = new TRandom3(seed_uniform_random_phi);
}
void TennGen::CloserFunction(){
  if(PRINTQA == kFALSE){
    vn_QA_file->Close();
    //delete vn_QA_file;
  }
  else if(PRINTQA == kTRUE){
    //vn_QA_plots->cd();
    v1_pi_h->Write();
    v2_pi_h->Write();
    v3_pi_h->Write();
    v4_pi_h->Write();
    v5_pi_h->Write();
    v1_K_h->Write();
    v2_K_h->Write();
    v3_K_h->Write();
    v4_K_h->Write();
    v5_K_h->Write();
    v1_P_h->Write();
    v2_P_h->Write();
    v3_P_h->Write();
    v4_P_h->Write();
    v5_P_h->Write();
    //Event_Plane_QA_plots->cd();
    v1_sq_h->Write();
    v2_sq_h->Write();
    v3_sq_h->Write();
    v4_sq_h->Write();
    v5_sq_h->Write();
    Psi_1_h->Write();
    Psi_2_h->Write();
    Psi_3_h->Write();
    Psi_4_h->Write();
    Psi_5_h->Write();
    vn_QA_file->Write();
    vn_QA_file->Close();
  }
}
void TennGen::MyPrivateFunction(){
  cout<<"MY PRIVATE FUNCTION"<<endl;
}

void TennGen::GetRandomSeed(){
  cout<<"\n\n\n\n\nCurrent Seed = "<<Seed<<"\n\n\n\n\n"<<endl;
}

/*
   UInt_t seed = (jobID+1)*17;
  if( (seed>=0) && (seed<=900000000) ) {
    pythia->SetMRPY(1, seed);                   // set seed
    pythia->SetMRPY(2, 0);                      // use new seed
    cout<<"Random Seed : "<<seed<<endl;
  } else {cout << "error: time " << seed << " is not valid" << endl; exit(2);}
   //____________________________________________________________________
  */
//TF1 *pXdist;

	
//____________FITTING___PARAMETER ARRAYS_____x_dim = centrality bin________y_dim = parameter________________________________//

///*

void TennGen::Setv1_Params(){

  const Double_t v1_pi_params[6][6] = { {-6.29928e-03 , 6.29213e-02 , -1.69603e-02 , 1.24408e-03 , 0.0 , 5.6864 } , { -8.15353e-03 ,  9.46349e-02 , -2.27846e-02 , 1.44930e-03 , 0.0 , 5.68688} , { -1.27599e-02 , 1.37523e-01 , -3.44225e-02 , 2.40203e-03 , 0.0 , 5.6862 } , { -1.66651e-02 , 1.85760e-01 , -4.85808e-02 , 3.52947e-03 , 0.0 ,  5.68643 } , {  -2.05442e-02 , 2.19814e-01 , -6.16465e-02 , 4.91569e-03 , 0.0 ,  5.68753 }  , { -2.08654e-02 , 2.34196e-01 , -6.95590e-02 , 5.88573e-03 , 0.0 , 5.68962 } };
  
  const Double_t v1_K_params[6][6] = {{ -2.28987e-02 , 6.91034e-02 , -1.48611e-02  ,  5.67949e-04 , 0.0 , 3.79629 } , { -3.28047e-02 ,  1.02630e-01  , -2.03472e-02 , 7.34138e-04 , 0.0 , 3.79652 } , { -4.44470e-02 , 1.43257e-01 , -2.79470e-02  , 8.88236e-04 , 0.0 , 3.79639 } , { -5.50614e-02 , 1.89560e-01 , -3.80298e-02 , 1.13837e-03 , 0.0 , 3.79629 } , { -6.46262e-02 ,  2.32903e-01 , -5.62942e-02 , 3.39183e-03 , 0.0 ,  3.79624 }  , { -6.68228e-02 , 2.58028e-01 , -7.21831e-02 ,  5.84160e-03 , 0.0 , 3.79576 } };
  
  const Double_t v1_P_params[6][6] = {{1.34354e-03 ,   -2.69847e-02 , 4.68928e-02 , -1.30276e-02 , 1.03297e-03 , 5.74363 } , { 1.68188e-02 , -5.83743e-02 , 7.93535e-02 , -2.04909e-02 , 1.53971e-03  , 5.74345 } , { 2.33248e-02 , -7.58826e-02 , 1.08507e-01 , -2.81839e-02 , 2.11781e-03 , 5.74522 } , { 1.17865e-02 , -5.86004e-02 , 1.24604e-01 , -3.45561e-02 , 2.68779e-03 ,  5.74372 } , { -1.40410e-03 , -3.01182e-02 , 1.26896e-01 , -3.80857e-02 , 3.12556e-03 ,  5.74289 }  , { -2.61424e-02 , 3.52815e-02 , 9.56710e-02  , -3.30046e-02 , 2.85671e-03 ,  5.74307 } };
  
  v1_pi->SetParameter(0,v1_pi_params[CentralityBin][0]);
  v1_pi->SetParameter(1,v1_pi_params[CentralityBin][1]);
  v1_pi->SetParameter(2,v1_pi_params[CentralityBin][2]);
  v1_pi->SetParameter(3,v1_pi_params[CentralityBin][3]);
  v1_pi->SetParameter(4,v1_pi_params[CentralityBin][4]);
  v1_pi->SetParameter(5,v1_pi_params[CentralityBin][5]);
  
  v1_K->SetParameter(0,v1_K_params[CentralityBin][0]);
  v1_K->SetParameter(1,v1_K_params[CentralityBin][1]);
  v1_K->SetParameter(2,v1_K_params[CentralityBin][2]);
  v1_K->SetParameter(3,v1_K_params[CentralityBin][3]);
  v1_K->SetParameter(4,v1_K_params[CentralityBin][4]);
  v1_K->SetParameter(5,v1_K_params[CentralityBin][5]);
  
  v1_P->SetParameter(0,v1_P_params[CentralityBin][0]);
  v1_P->SetParameter(1,v1_P_params[CentralityBin][1]);
  v1_P->SetParameter(2,v1_P_params[CentralityBin][2]);
  v1_P->SetParameter(3,v1_P_params[CentralityBin][3]);
  v1_P->SetParameter(4,v1_P_params[CentralityBin][4]);
  v1_P->SetParameter(5,v1_P_params[CentralityBin][5]);
}

void TennGen::Setv2_Params(){

  const Double_t v2_pi_params[6][6] = { {-6.29928e-03 , 6.29213e-02 , -1.69603e-02 , 1.24408e-03 , 0.0 , 5.6864 } , { -8.15353e-03 ,  9.46349e-02 , -2.27846e-02 , 1.44930e-03 , 0.0 , 5.68688} , { -1.27599e-02 , 1.37523e-01 , -3.44225e-02 , 2.40203e-03 , 0.0 , 5.6862 } , { -1.66651e-02 , 1.85760e-01 , -4.85808e-02 , 3.52947e-03 , 0.0 ,  5.68643 } , {  -2.05442e-02 , 2.19814e-01 , -6.16465e-02 , 4.91569e-03 , 0.0 ,  5.68753 }  , { -2.08654e-02 , 2.34196e-01 , -6.95590e-02 , 5.88573e-03 , 0.0 , 5.68962 } };
  
  const Double_t v2_K_params[6][6] = { { -2.28987e-02 , 6.91034e-02 , -1.48611e-02  ,  5.67949e-04 , 0.0 , 3.79629 } , { -3.28047e-02 ,  1.02630e-01  , -2.03472e-02 , 7.34138e-04 , 0.0 , 3.79652 } , { -4.44470e-02 , 1.43257e-01 , -2.79470e-02  , 8.88236e-04 , 0.0 , 3.79639 } , { -5.50614e-02 , 1.89560e-01 , -3.80298e-02 , 1.13837e-03 , 0.0 , 3.79629 } , { -6.46262e-02 ,  2.32903e-01 , -5.62942e-02 , 3.39183e-03 , 0.0 ,  3.79624 }  , { -6.68228e-02 , 2.58028e-01 , -7.21831e-02 ,  5.84160e-03 , 0.0 , 3.79576 } };
  
  const Double_t v2_P_params[6][6] = { {1.34354e-03 ,   -2.69847e-02 , 4.68928e-02 , -1.30276e-02 , 1.03297e-03 , 5.74363 } , { 1.68188e-02 , -5.83743e-02 , 7.93535e-02 , -2.04909e-02 , 1.53971e-03  , 5.74345 } , { 2.33248e-02 , -7.58826e-02 , 1.08507e-01 , -2.81839e-02 , 2.11781e-03 , 5.74522 } , { 1.17865e-02 , -5.86004e-02 , 1.24604e-01 , -3.45561e-02 , 2.68779e-03 ,  5.74372 } , { -1.40410e-03 , -3.01182e-02 , 1.26896e-01 , -3.80857e-02 , 3.12556e-03 ,  5.74289 }  , { -2.61424e-02 , 3.52815e-02 , 9.56710e-02  , -3.30046e-02 , 2.85671e-03 ,  5.74307 } };
  
  v2_pi->SetParameter(0,v2_pi_params[CentralityBin][0]);
  v2_pi->SetParameter(1,v2_pi_params[CentralityBin][1]);
  v2_pi->SetParameter(2,v2_pi_params[CentralityBin][2]);
  v2_pi->SetParameter(3,v2_pi_params[CentralityBin][3]);
  v2_pi->SetParameter(4,v2_pi_params[CentralityBin][4]);
  v2_pi->SetParameter(5,v2_pi_params[CentralityBin][5]);
  
  v2_K->SetParameter(0,v2_K_params[CentralityBin][0]);
  v2_K->SetParameter(1,v2_K_params[CentralityBin][1]);
  v2_K->SetParameter(2,v2_K_params[CentralityBin][2]);
  v2_K->SetParameter(3,v2_K_params[CentralityBin][3]);
  v2_K->SetParameter(4,v2_K_params[CentralityBin][4]);
  v2_K->SetParameter(5,v2_K_params[CentralityBin][5]);
  
  v2_P->SetParameter(0,v2_P_params[CentralityBin][0]);
  v2_P->SetParameter(1,v2_P_params[CentralityBin][1]);
  v2_P->SetParameter(2,v2_P_params[CentralityBin][2]);
  v2_P->SetParameter(3,v2_P_params[CentralityBin][3]);
  v2_P->SetParameter(4,v2_P_params[CentralityBin][4]);
  v2_P->SetParameter(5,v2_P_params[CentralityBin][5]);

}

void TennGen::Setv3_Params(){

  const Double_t v3_pi_params[6][6] = { { -2.53410e-02 , 7.65854e-02  , -1.64087e-02 , 9.27325e-04 , 0.0 , 5.6864  } , { -2.47587e-02  ,  8.05937e-02 , -1.60442e-02 , 7.45643e-04 , 0.0 , 5.68688 } , {  -2.40356e-02 , 8.39185e-02 ,  -1.60066e-02 , 6.23750e-04 , 0.0, 5.6862 } , { -2.75520e-02  ,  9.74475e-02 , -2.11485e-02 , 1.13159e-03 , 0.0 , 5.68643 } , { -2.89726e-02  , 1.05326e-01 , -2.52450e-02 ,  1.65368e-03 , 0.0 , 5.68753 }  , { -2.29852e-02  , 9.41400e-02  , -1.86962e-02 , -5.14131e-04  , 2.43303e-04 , 5.68962 } };
 
  const Double_t v3_K_params[6][6] = { { -2.00467e-02 , 4.04763e-02 , 4.87640e-03 , -2.25452e-03 , 0.0 , 3.79629 } , {  -2.12149e-02 , 4.44332e-02 , 5.71953e-03 , -2.55308e-03 , 0.0 , 3.79652  } , { -2.35466e-02 , 5.22392e-02 , 4.10174e-03  ,  -2.57112e-03 , 0.0 , 3.79639 } , { -2.76820e-02 , 6.67519e-02 , -1.63168e-03  , -2.03702e-03 , 0.0 , 3.79629  } , { -3.02956e-02 , 7.58570e-02  , -5.81264e-03 , -1.57544e-03  , 0.0 , 3.79624 }  , { -2.46987e-02 , 6.96050e-02 , -4.33787e-03 , -1.83356e-03 , 0.0 ,  3.79576 } };
  
  const Double_t v3_P_params[6][6] = { { -4.69818e-03 , -1.24501e-02 , 2.67730e-02 ,  -4.02833e-03 , 0.0 , 5.74363 } , { -9.67248e-03  , -4.64947e-03 , 2.65129e-02 , -4.13887e-03 , 0.0 , 5.74345 } , { -8.40764e-03 , -9.29454e-03  , 3.21732e-02 ,  -5.16580e-03 , 0.0 , 5.74522 } , { -2.37472e-02 , 1.73042e-02 , 2.41911e-02 , -4.40537e-03 , 0.0 , 5.74372 } , { -2.26647e-02 , 2.18590e-02 , 2.34411e-02 , -4.53662e-03  , 0.0 , 5.74289 }  , { -4.06477e-02 , 5.58970e-02 ,  9.23320e-03 , -2.90313e-03 , 0.0 ,  5.74307 } };
  
  v3_pi->SetParameter(0,v3_pi_params[CentralityBin][0]);
  v3_pi->SetParameter(1,v3_pi_params[CentralityBin][1]);
  v3_pi->SetParameter(2,v3_pi_params[CentralityBin][2]);
  v3_pi->SetParameter(3,v3_pi_params[CentralityBin][3]);
  v3_pi->SetParameter(4,v3_pi_params[CentralityBin][4]);
  v3_pi->SetParameter(5,v3_pi_params[CentralityBin][5]);
  
  v3_K->SetParameter(0,v3_K_params[CentralityBin][0]);
  v3_K->SetParameter(1,v3_K_params[CentralityBin][1]);
  v3_K->SetParameter(2,v3_K_params[CentralityBin][2]);
  v3_K->SetParameter(3,v3_K_params[CentralityBin][3]);
  v3_K->SetParameter(4,v3_K_params[CentralityBin][4]);
  v3_K->SetParameter(5,v3_K_params[CentralityBin][5]);
  
  v3_P->SetParameter(0,v3_P_params[CentralityBin][0]);
  v3_P->SetParameter(1,v3_P_params[CentralityBin][1]);
  v3_P->SetParameter(2,v3_P_params[CentralityBin][2]);
  v3_P->SetParameter(3,v3_P_params[CentralityBin][3]);
  v3_P->SetParameter(4,v3_P_params[CentralityBin][4]);
  v3_P->SetParameter(5,v3_P_params[CentralityBin][5]);
}

void TennGen::Setv4_Params(){

  const Double_t v4_pi_params[6][6] = { { -2.61775e-02 , 5.41027e-02 , -8.11629e-03 , 1.07636e-04 , 0.0 , 5.56759 } , { -2.84327e-02  , 6.15741e-02 , -1.08541e-02 , 3.74468e-04 , 0.0 , 5.56826 } , { -2.50069e-02 ,  5.64826e-02  , -7.82525e-03 , -5.60331e-05 , 0.0 , 5.56784 } , { -3.66736e-02 , 8.52490e-02 , -2.18799e-02 , 1.76828e-03 , 0.0 ,  5.56903 } , { -2.67396e-02 , 6.76440e-02 , -1.25867e-02 , 3.54348e-04 , 0.0 , 5.57063 }  , { -2.83923e-02 ,  7.37708e-02 , -1.73843e-02  , 1.07088e-03 , 0.0 , 5.57323 } };
  
  const Double_t v4_K_params[6][6] = { { -1.96656e-02 , 2.61057e-02 , 7.38473e-03  , -2.14579e-03 , 0.0 , 3.88195 } , { -2.09796e-02  , 2.93941e-02 ,  6.44392e-03 , -2.01671e-03 , 0.0 , 3.88229 } , { -2.32744e-02 , 3.62767e-02  ,  3.93014e-03  , -1.89347e-03 , 0.0 , 3.88197 } , { -2.02513e-02 , 2.93811e-02 ,  1.03741e-02 , -3.26730e-03 , 0.0 , 3.8817 } , { -1.45235e-02 , 1.84089e-02 , 1.63240e-02 , -4.16862e-03 , 0.0 , 3.88098 }  , { -1.96178e-02 , 3.49596e-02 , 5.84558e-03 , -2.54327e-03 , 0.0 , 3.8814 } };
  
  const Double_t v4_P_params[6][6] = { {  1.68565e-02 , -4.84043e-02 , 3.57252e-02 , -4.69793e-03 , 0.0 , 5.62931 } , { 1.84244e-02  , -5.34851e-02 , 4.00556e-02 , -5.30060e-03 , 0.0 , 5.62921 } , { 1.53098e-02 , -4.93400e-02 , 3.99301e-02 , -5.45571e-03 , 0.0 , 5.63246 } , { 1.31460e-02 , -4.65074e-02 , 4.21324e-02 , -5.94997e-03 , 0.0 ,  5.6322 } , { 9.26040e-03 , -4.29033e-02 , 4.25759e-02 , -6.31299e-03 , 0.0 , 5.6311 }  , { -2.62813e-02 , 2.86358e-02  , 8.47114e-03 , -1.90998e-03 , 0.0 , 5.63034 } };
  
  v4_pi->SetParameter(0,v4_pi_params[CentralityBin][0]);
  v4_pi->SetParameter(1,v4_pi_params[CentralityBin][1]);
  v4_pi->SetParameter(2,v4_pi_params[CentralityBin][2]);
  v4_pi->SetParameter(3,v4_pi_params[CentralityBin][3]);
  v4_pi->SetParameter(4,v4_pi_params[CentralityBin][4]);
  v4_pi->SetParameter(5,v4_pi_params[CentralityBin][5]);
  
  v4_K->SetParameter(0,v4_K_params[CentralityBin][0]);
  v4_K->SetParameter(1,v4_K_params[CentralityBin][1]);
  v4_K->SetParameter(2,v4_K_params[CentralityBin][2]);
  v4_K->SetParameter(3,v4_K_params[CentralityBin][3]);
  v4_K->SetParameter(4,v4_K_params[CentralityBin][4]);
  v4_K->SetParameter(5,v4_K_params[CentralityBin][5]);
  
  v4_P->SetParameter(0,v4_P_params[CentralityBin][0]);
  v4_P->SetParameter(1,v4_P_params[CentralityBin][1]);
  v4_P->SetParameter(2,v4_P_params[CentralityBin][2]);
  v4_P->SetParameter(3,v4_P_params[CentralityBin][3]);
  v4_P->SetParameter(4,v4_P_params[CentralityBin][4]);
  v4_P->SetParameter(5,v4_P_params[CentralityBin][5]);

}

void TennGen::Setv5_Params(){

  const Double_t v5_pi_params[6][6] = { {  -2.27633e-02 , 3.62812e-02 , -5.09327e-03 , 3.04512e-06 , 0.0 , 5.3046 } , { -2.11594e-02 , 3.40405e-02 , -4.51201e-03 ,  1.65056e-04 , 0.0 , 5.30516 } , { -1.79688e-02 , 2.73746e-02 , 1.11420e-03  , -1.03606e-03 , 0.0 , 5.30769 } , { -1.64322e-02 , 2.88779e-02 , 3.15650e-05 , -6.48281e-04  , 0.0 , 5.31082 } , { -2.00310e-02 , 3.98351e-02 , -7.27710e-03  , 5.18693e-04, 0.0 ,  5.31308 }  , {  -2.14945e-02 , 4.30007e-02  , -8.88504e-03 , 4.45375e-04 , 0.0 , 5.31678 } };
  
  const Double_t v5_K_params[6][6] = { { -2.93359e-02 , 4.03216e-02 , -6.72350e-03 , 2.49656e-04 , 0.0 , 3.68094 } , { -2.45592e-02 ,  3.08562e-02 , -1.20602e-03 , -5.88864e-04 , 0.0 , 3.68144 } , { -2.23697e-02 , 2.31996e-02 , 5.23514e-03, -1.79606e-03 , 0.0 , 3.68252 } , { -2.38504e-02 , 3.07280e-02 , 3.01792e-04 , -1.04716e-03 , 0.0 , 3.68158 } , { -3.33063e-02 , 5.52523e-02 , -1.62386e-02 , 1.98994e-03 , 0.0 , 3.68023 }  , { -3.34578e-02 , 5.08559e-02 , -1.28618e-02 , 1.32781e-03 , 0.0 , 3.68031 } };
  
  const Double_t v5_P_params[6][6] = { { 3.20157e-02 ,  -8.18906e-02 , 5.85698e-02 ,  -1.34328e-02  ,  1.01128e-03  , 5.37351 } , { 1.08684e-02 , -3.07064e-02 , 2.07656e-02 , -2.71231e-03  , 4.65723e-05 , 5.37475 } , { 3.93613e-02 , -1.01226e-01 , 7.84549e-02 , -1.98730e-02 , 1.67558e-03  , 5.37907 } , { 5.87097e-02 ,  -1.42471e-01 , 1.07133e-01 , -2.72022e-02 ,  2.30510e-03 , 5.3852 } , {  5.54421e-02 , -1.37351e-01 , 1.02622e-01 , -2.52028e-02  , 1.99010e-03 ,  5.3828 }  , { 7.55687e-02 , -2.06518e-01 , 1.57785e-01 , -4.01594e-02 , 3.28859e-03 ,  5.38206 } };
  
  v5_pi->SetParameter(0,v5_pi_params[CentralityBin][0]);
  v5_pi->SetParameter(1,v5_pi_params[CentralityBin][1]);
  v5_pi->SetParameter(2,v5_pi_params[CentralityBin][2]);
  v5_pi->SetParameter(3,v5_pi_params[CentralityBin][3]);
  v5_pi->SetParameter(4,v5_pi_params[CentralityBin][4]);
  v5_pi->SetParameter(5,v5_pi_params[CentralityBin][5]);
  
  v5_K->SetParameter(0,v5_K_params[CentralityBin][0]);
  v5_K->SetParameter(1,v5_K_params[CentralityBin][1]);
  v5_K->SetParameter(2,v5_K_params[CentralityBin][2]);
  v5_K->SetParameter(3,v5_K_params[CentralityBin][3]);
  v5_K->SetParameter(4,v5_K_params[CentralityBin][4]);
  v5_K->SetParameter(5,v5_K_params[CentralityBin][5]);
  
  v5_P->SetParameter(0,v5_P_params[CentralityBin][0]);
  v5_P->SetParameter(1,v5_P_params[CentralityBin][1]);
  v5_P->SetParameter(2,v5_P_params[CentralityBin][2]);
  v5_P->SetParameter(3,v5_P_params[CentralityBin][3]);
  v5_P->SetParameter(4,v5_P_params[CentralityBin][4]);
  v5_P->SetParameter(5,v5_P_params[CentralityBin][5]);

}

//*/

void TennGen::SetpiPlusParams(){

  const Double_t piPlus_params[8][5] = { { 0.139570 , 0.947908 , 0.0737422 , 0.968789 , 2207060 } , { 0.139570 , 0.943360 , 0.078973 , 1.02539 , 1500940 } , { 0.139570 , 0.95381 , 0.0708164 , 0.987509 , 1603690 } , { 0.139570 , 0.95409 , 0.0673556 , 1.00806 , 1319880 } , { 0.139570 , 0.959811 , 0.0676647 , 1.08715 , 861341 } , { 0.139570 , 0.957266 , 0.08703763 , 1.21600, 468819 } , { 0.139570 , 0.949426 , 0.078257 , 1.46891 , 190742 } , { 0.139570 , 0.954145 , 0.0744700 , 1.58732, 119812 }};

  pTdist_piPlus->SetParameter(0,piPlus_params[CentralityBin][0]);
  pTdist_piPlus->SetParameter(1,piPlus_params[CentralityBin][1]);
  pTdist_piPlus->SetParameter(2,piPlus_params[CentralityBin][2]);
  pTdist_piPlus->SetParameter(3,piPlus_params[CentralityBin][3]);
  pTdist_piPlus->SetParameter(4,piPlus_params[CentralityBin][4]);

}

void TennGen::SetpiMinusParams(){

  const Double_t piMinus_params[8][5] = { { 0.139570 , 0.942100 , 0.072976 , 0.967010 , 1856920 } , { 0.139570 , 0.944275 , 0.0766907 , 0.986230 , 1583560 } , { 0.139570 , 0.945046 , 0.072708 , 1.01734 , 1160050 } , { 0.139570 , 0.952089 , 0.0729955 , 1.03405 , 970911 } , { 0.139570 , 0.952325 , 0.0737461 , 1.12383 , 620421 } , { 0.139570 , 0.940441 , 0.0844949 , 1.39199, 240963 } , { 0.139570 , 0.933892 , 0.0902015 , 1.63995 , 112371 } , { 0.139570 , 0.927508 , 0.0966857 , 2.00, 47185.9 }};

  pTdist_piMinus->SetParameter(0,piMinus_params[CentralityBin][0]);
  pTdist_piMinus->SetParameter(1,piMinus_params[CentralityBin][1]);
  pTdist_piMinus->SetParameter(2,piMinus_params[CentralityBin][2]);
  pTdist_piMinus->SetParameter(3,piMinus_params[CentralityBin][3]);
  pTdist_piMinus->SetParameter(4,piMinus_params[CentralityBin][4]);

}

void TennGen::SetkPlusParams(){

  const Double_t kPlus_params[8][5] = { { 0.493677 , 0.826760 , 0.129893 , 0.703341 , 230444 } , { 0.493677 , 0.820456 , 0.135416 , 0.754479 , 154330 } , { 0.493677 , 0.822368 , 0.136030 , 0.790558 , 113729 } , { 0.493677 , 0.812882 , 0.14257 , 0.868716 , 60544.6 } , { 0.493677 , 0.808265 , 0.148363 , 1.01543 , 32570.8 } , { 0.493677 , 0.741293 , 0.185352 , 1.36225, 7164.93 } , { 0.493677 , 0.771795 , 0.180722, 1.77221 , 4536.61 } , { 0.493677 , 0.823896 , 0.156262 , 1.94719, 4555.18 }};

  pTdist_kPlus->SetParameter(0,kPlus_params[CentralityBin][0]);
  pTdist_kPlus->SetParameter(1,kPlus_params[CentralityBin][1]);
  pTdist_kPlus->SetParameter(2,kPlus_params[CentralityBin][2]);
  pTdist_kPlus->SetParameter(3,kPlus_params[CentralityBin][3]);
  pTdist_kPlus->SetParameter(4,kPlus_params[CentralityBin][4]);

}

void TennGen::SetkMinusParams(){

  const Double_t kMinus_params[8][5] = { { 0.493677 , 0.821301 , 0.133049 , 0.748879 , 203970 } , { 0.493677 , 0.821564 , 0.134495 , 0.781053 , 160695 } , { 0.493677 , 0.816101 , 0.138573 , 0.816978 , 103440 } , { 0.493677 , 0.807483 , 0.146798 , 0.936718 , 52652.1 } , { 0.493677 , 0.802314 , 0.152261 , 1.08377 , 28720.5 } , { 0.493677 , 0.777818 , 0.171125 , 1.41012, 10213.3 } , { 0.493677 , 0.805967 , 0.162440 , 1.76754 , 7498.70 } , { 0.493677 , 0.830498 , 0.151675 , 1.97103, 5279.26 }};
  
  pTdist_kMinus->SetParameter(0,kMinus_params[CentralityBin][0]);
  pTdist_kMinus->SetParameter(1,kMinus_params[CentralityBin][1]);
  pTdist_kMinus->SetParameter(2,kMinus_params[CentralityBin][2]);
  pTdist_kMinus->SetParameter(3,kMinus_params[CentralityBin][3]);
  pTdist_kMinus->SetParameter(4,kMinus_params[CentralityBin][4]);

}

void TennGen::SetpParams(){

  const Double_t p_params[8][5] = { { 0.938272 , 0.861255 , 0.113053 , 0.646603 , 3540660 } , { 0.938272 , 0.829143 , 0.134133 , 0.583048 , 608721 } , { 0.93827 , 0.784519 , 0.158416 , 0.489724 , 121635 } , { 0.938272 , 0.755180 , 0.170715 , 0.472479 , 49103.6 } , { 0.938272 , 0.738171 , 0.178862 , 0.549364 , 24070.7 } , { 0.938272 , 0.718834 , 0.188403 , 0.698602 , 11050.3 } , { 0.938272 , 0.685109 , 0.204819 , 0.914109 , 3943.15 } , { 0.938272 , 0.638029 , 0.244535 , 1.78914 , 768.023 }};
  
        
  pTdist_p->SetParameter(0,p_params[CentralityBin][0]);
  pTdist_p->SetParameter(1,p_params[CentralityBin][1]);
  pTdist_p->SetParameter(2,p_params[CentralityBin][2]);
  pTdist_p->SetParameter(3,p_params[CentralityBin][3]);
  pTdist_p->SetParameter(4,p_params[CentralityBin][4]);

}

void TennGen::SetpbarParams(){

  const Double_t pbar_params[8][5] = { { 0.938272 , 0.863601 , 0.11476 , 0.658438 , 4073300 } , { 0.938272 , 0.855439 , 0.116446 , 0.647756 , 2212850 } , { 0.938272 , 0.832520 , 0.131665 , 0.622898 , 552092 } , { 0.938272 , 0.778317 , 0.160757 , 0.539450 , 77034.0 } , { 0.938272 , 0.739188 , 0.174989 , 0.532042 , 28702.3 } , { 0.938272 , 0.718917 , 0.186137 , 0.680732 , 12351.2 } , { 0.938272 , 0.688116 , 0.202278 , 0.939983 , 4397.05 } , { 0.938272 , 0.637255 , 0.244434 , 1.83999 , 799.620 }};
  
  pTdist_pbar->SetParameter(0,pbar_params[CentralityBin][0]);
  pTdist_pbar->SetParameter(1,pbar_params[CentralityBin][1]);
  pTdist_pbar->SetParameter(2,pbar_params[CentralityBin][2]);
  pTdist_pbar->SetParameter(3,pbar_params[CentralityBin][3]);
  pTdist_pbar->SetParameter(4,pbar_params[CentralityBin][4]);

}

void TennGen::SetpiZeroParams(){

  const Double_t piZero_params[8][5] = { { 0.139977 , 0.942100 , 0.072976 , 0.967010 , 1856920 } , { 0.139977 , 0.944275 , 0.0766907 , 0.986230 , 1583560 } , { 0.139977, 0.945046 , 0.072708 , 1.01734 , 1160050 } , { 0.139977 , 0.952089 , 0.0729955 , 1.03405 , 970911 } , { 0.139977 , 0.952325 , 0.0737461 , 1.12383 , 620421 } , { 0.139977 , 0.940441 , 0.0844949 , 1.39199, 240963 } , { 0.139977 , 0.933892 , 0.0902015 , 1.63995 , 112371 } , { 0.139977 , 0.927508 , 0.0966857 , 2.00, 47185.9 }};

  pTdist_piZero->SetParameter(0,piZero_params[CentralityBin][0]);
  pTdist_piZero->SetParameter(1,piZero_params[CentralityBin][1]);
  pTdist_piZero->SetParameter(2,piZero_params[CentralityBin][2]);
  pTdist_piZero->SetParameter(3,piZero_params[CentralityBin][3]);
  pTdist_piZero->SetParameter(4,piZero_params[CentralityBin][4]);

}


TClonesArray *TennGen::GetBackground(){
  cout<<"Got Background"<<endl;

//TF1 *TennGen::Random1(){

//___________________________________________________________________________________________________________//
//____________YIELD__ARRAYS_____x_dim = centrality bin________y_dim = PARTICLE________________________________//
//____________y_diM -- INDEX___ 0 = pi+____1 = pi-____2 = K-____3 = K+____4 = p____5 = pbar_____6 = pi0_________//
//const Double_t yield_arr[8][7] = { { (TMath::Ceil((EtaRange * 733)/0.5)) , (TMath::Ceil((EtaRange * 732)/0.5)) , (TMath::Ceil((EtaRange * 109)/0.5)) , (TMath::Ceil((EtaRange * 109)/0.5)) , (TMath::Ceil((EtaRange * 34)/0.5)) , (TMath::Ceil((EtaRange * 33)/0.5)) , (TMath::Ceil((EtaRange * 732)/0.5)) } , { (TMath::Ceil((EtaRange * 606)/0.5)) , (TMath::Ceil((EtaRange * 604)/0.5)) , (TMath::Ceil((EtaRange * 91)/0.5)) , (TMath::Ceil((EtaRange * 90)/0.5)) , (TMath::Ceil((EtaRange * 28)/0.5)) , (TMath::Ceil((EtaRange * 28)/0.5)) , (TMath::Ceil((EtaRange * 604)/0.5)) } , { (TMath::Ceil((EtaRange * 455)/0.5)) , (TMath::Ceil((EtaRange * 453)/0.5)) , (TMath::Ceil((EtaRange * 68)/0.5)) , (TMath::Ceil((EtaRange * 68)/0.5)) , (TMath::Ceil((EtaRange * 21.0)/0.5)) , (TMath::Ceil((EtaRange * 21)/0.5)) , (TMath::Ceil((EtaRange * 453)/0.5)) } , { (TMath::Ceil((EtaRange * 307)/0.5)) , (TMath::Ceil((EtaRange * 306 )/0.5)), (TMath::Ceil((EtaRange * 46)/0.5)) , (TMath::Ceil((EtaRange * 46)/0.5)) , (TMath::Ceil((EtaRange * 14)/0.5)) , (TMath::Ceil((EtaRange * 14)/0.5)) , (TMath::Ceil((EtaRange * 306)/0.5)) } , { (TMath::Ceil((EtaRange * 201)/0.5)) , (TMath::Ceil((EtaRange * 200)/0.5)) , (TMath::Ceil((EtaRange * 30)/0.5)) , (TMath::Ceil((EtaRange * 30)/0.5)) , (TMath::Ceil((EtaRange * 10)/0.5)) , (TMath::Ceil((EtaRange * 10)/0.5)) , (TMath::Ceil((EtaRange * 200)/0.5)) } , { (TMath::Ceil((EtaRange * 124)/0.5)) , (TMath::Ceil((EtaRange * 123)/0.5)) , (TMath::Ceil((EtaRange * 18)/0.5)) , (TMath::Ceil(EtaRange/0.5)) * 18 , (TMath::Ceil(EtaRange/0.5)) * 6 , (TMath::Ceil(EtaRange/0.5)) * 6 , (TMath::Ceil((EtaRange * 123)/0.5)) } , { (TMath::Ceil((EtaRange * 71)/0.5)) , (TMath::Ceil((EtaRange * 71)/0.5)) , (TMath::Ceil((EtaRange * 10)/0.5)) , (TMath::Ceil((EtaRange * 10)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 71)/0.5)) } , { (TMath::Ceil((EtaRange * 37)/0.5)) , (TMath::Ceil((EtaRange * 37)/0.5)) , ((TMath::Ceil(EtaRange * 5)/0.5)) , (TMath::Ceil((EtaRange * 5)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 37)/0.5)) } };

const Double_t yield_arr[8][7] = { { (TMath::Ceil((EtaRange * 654)/0.5)) , (TMath::Ceil((EtaRange * 654)/0.5)) , (TMath::Ceil((EtaRange * 97)/0.5)) , (TMath::Ceil((EtaRange * 97)/0.5)) , (TMath::Ceil((EtaRange * 30)/0.5)) , (TMath::Ceil((EtaRange * 29)/0.5)) , (TMath::Ceil((EtaRange * 654)/0.5)) } , { (TMath::Ceil((EtaRange * 541)/0.5)) , (TMath::Ceil((EtaRange * 539)/0.5)) , (TMath::Ceil((EtaRange * 81)/0.5)) , (TMath::Ceil((EtaRange * 80)/0.5)) , (TMath::Ceil((EtaRange * 25)/0.5)) , (TMath::Ceil((EtaRange * 25)/0.5)) , (TMath::Ceil((EtaRange * 539)/0.5)) } , { (TMath::Ceil((EtaRange * 406)/0.5)) , (TMath::Ceil((EtaRange * 404)/0.5)) , (TMath::Ceil((EtaRange * 61)/0.5)) , (TMath::Ceil((EtaRange * 61)/0.5)) , (TMath::Ceil((EtaRange * 19)/0.5)) , (TMath::Ceil((EtaRange * 19)/0.5)) , (TMath::Ceil((EtaRange * 404)/0.5)) } , { (TMath::Ceil((EtaRange * 274)/0.5)) , (TMath::Ceil((EtaRange * 273 )/0.5)), (TMath::Ceil((EtaRange * 41)/0.5)) , (TMath::Ceil((EtaRange * 41)/0.5)) , (TMath::Ceil((EtaRange * 13)/0.5)) , (TMath::Ceil((EtaRange * 13)/0.5)) , (TMath::Ceil((EtaRange * 273)/0.5)) } , { (TMath::Ceil((EtaRange * 179)/0.5)) , (TMath::Ceil((EtaRange * 179)/0.5)) , (TMath::Ceil((EtaRange * 27)/0.5)) , (TMath::Ceil((EtaRange * 27)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 179)/0.5)) } , { (TMath::Ceil((EtaRange * 111)/0.5)) , (TMath::Ceil((EtaRange * 110)/0.5)) , (TMath::Ceil((EtaRange * 16)/0.5)) , (TMath::Ceil((EtaRange * 16)/0.5)) , (TMath::Ceil((EtaRange * 5)/0.5)) , (TMath::Ceil((EtaRange * 5)/0.5)) , 	(TMath::Ceil((EtaRange * 110)/0.5)) } , { (TMath::Ceil((EtaRange * 63)/0.5)) , (TMath::Ceil((EtaRange * 63)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 63)/0.5)) } , { (TMath::Ceil((EtaRange * 33)/0.5)) , (TMath::Ceil((EtaRange * 33)/0.5)) , ((TMath::Ceil(EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 33)/0.5)) } }; //setting multiplicities for particles



//const Double_t yield_arr[8][7] = { { 5000 , 5000 , 5000 , 5000 , 5000 , 5000 , 5000 } , { 606 , 604 , 91 , 90 , 28 , 28 , 604 } , { 455 , 453 , 68 , 68 , 21.0 , 21 , 453 } , { 307 , 306 , 46 , 46 , 14 , 14 , 306 } , { 201 , 200 , 30 , 30 , 10 , 10 , 200 } , { 124 , 123 , 18 , 18 , 6 , 6 , 123 } , { 71 , 71 , 10 , 10 , 4 , 4 , 71 } , { 37 , 37 , 5 , 5 , 2 , 2 , 37 } }; 

  Double_t v2_pi_eval1;
  Double_t v3_pi_eval1;
  Double_t v4_pi_eval1;
  Double_t v5_pi_eval1;

  Double_t v2_K_eval1;
  Double_t v3_K_eval1;
  Double_t v4_K_eval1;
  Double_t v5_K_eval1;

  Double_t v2_P_eval1;
  Double_t v3_P_eval1;
  Double_t v4_P_eval1;
  Double_t v5_P_eval1;


/*
  if (etadist!=0) delete etadist;
  if (Psi_1!=0) delete Psi_1;
  if (Psi_3!=0) delete Psi_3;
  if (Psi_5!=0) delete Psi_5;
  if (Harmonics_Phi_Dist_Rand!=0) delete Harmonics_Phi_Dist_Rand;
  if (New_Rand!=0) delete New_Rand;

  New_Rand = new TRandom3(Seed);

  etadist = new TF1("etadist","[0]",-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
  etadist->SetParameter(0,1.);
  //cout << "\n\n\nRandom etadist Seed from .cxx = "<< (etadist->GetSeed()) <<" for Event " << endl;
  
  Psi_1 = new TRandom3(New_Rand->Uniform(0 , 1000));
  cout << "Random Psi_1 Seed from .cxx = "<< (Psi_1->GetSeed()) <<" for Event " << endl;

  Psi_3 = new TRandom3(New_Rand->Uniform(0 , 1000));
  cout << "Random Psi_3 Seed from .cxx = "<< (Psi_3->GetSeed()) <<" for Event " << endl;

  Psi_5 = new TRandom3(New_Rand->Uniform(0 , 1000));
  cout << "Random Psi_5 Seed from .cxx = "<< (Psi_5->GetSeed()) <<" for Event " << endl;

  Harmonics_Phi_Dist_Rand = new TRandom3(New_Rand->Uniform(0 , 1000));
  cout << "Random Harmonics_Phi_Dist_Rand Seed from .cxx = "<< (Harmonics_Phi_Dist_Rand->GetSeed()) <<" for Event \n\n\n" << endl;
  */
  //Psi_1_event.clear();
  //Psi_3_event.clear();
  //Psi_5_event.clear();
  Double_t Psi_1_event = Psi_1->Uniform( 0 , 2.0*TMath::Pi());
  cout << "\n\n\nPsi_1 from .cxx = "<< Psi_1_event <<" for Event " << endl;
  cout << "Random Psi_1 Seed from .cxx = "<< (Psi_1->GetSeed()) <<" for Event " << endl;
  Psi_1_h->Fill(Psi_1_event);
  Psi_2_h->Fill(0.);
  //Psi_1_event = (TMath::Pi()/2);
  Double_t Psi_3_event = Psi_3->Uniform( 0 , 2.0*TMath::Pi());
  cout << "Psi_3 from .cxx = "<< Psi_3_event <<" for Event " << endl;
  cout << "Random Psi_3 Seed from .cxx = "<< (Psi_3->GetSeed()) <<" for Event " << endl;
  Psi_3_h->Fill(Psi_3_event);
  Psi_4_h->Fill(0.);
  //Psi_3_event = (3*TMath::Pi()/2);
  Double_t Psi_5_event = Psi_5->Uniform( 0 , 2.0*TMath::Pi());
  cout << "Psi_5 from .cxx = "<< Psi_5_event <<" for Event " << endl;
  cout << "Random Psi_5 Seed from .cxx = "<< (Psi_5->GetSeed()) <<" for Event \n\n\n" << endl;
  Psi_5_h->Fill(Psi_5_event);
  //Psi_5_event = (5*TMath::Pi()/2);

  Int_t Npart = 0;
  for (Int_t n = 0 ; n < 7 ; n++){
    Npart = Npart + yield_arr[CentralityBin][n];
  }
  //cout<<"Number of Total Particles for Centrality Bin "<<CentralityBin<<" = "<<Npart<<endl;
  cout<<"\n\n\n\n\nSeed From the .cxx = "<<Seed<<"\n\n\n\n\n"<<endl;
  for( Int_t N = 0 ; N < yield_arr[CentralityBin][0] ; N++){ //piPlus , kF = 211
     
      if( N==0 ){
        //cout<<"piPlus getting created now"<<endl;
      }  
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
  
      TMCParticle *piPlus = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,piPlus_params[CentralityBin][0]);
      pTdist->SetParameter(1,piPlus_params[CentralityBin][1]);
      pTdist->SetParameter(2,piPlus_params[CentralityBin][2]);
      pTdist->SetParameter(3,piPlus_params[CentralityBin][3]);
      pTdist->SetParameter(4,piPlus_params[CentralityBin][4]);
      */
      //cout<<"Params for iteration "<<N<<" set now"<<endl;
      Double_t pT = pTdist_piPlus->GetRandom();

/*_________________________________________________________________________________*/

      if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_pi_eval1 = 0;
      }

      else if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_pi_eval1 = (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_pi_eval1 = 0;
      }

      else if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_pi_eval1 = (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_pi_eval1 = 0;
      }

      else if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_pi_eval1 = (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_pi_eval1 = 0;
      }

      else if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_pi_eval1 = (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/
      //cout<<"pT for iteration "<<N<<" set now"<<endl;
      //phi.clear();
      if( HarmonicsFlag != 5){
        Set_Phi_Dist(pT , 211 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for Pi+ now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_pi_h->Fill( pT , v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_pi_h->Fill( pT , v2_pi_eval1 );
        v3_pi_h->Fill( pT , v3_pi_eval1 );
        v4_pi_h->Fill( pT , v4_pi_eval1 );
        v5_pi_h->Fill( pT , v5_pi_eval1 );

        v1_sq_h->Fill( v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_pi_eval1*v2_pi_eval1 );
        v3_sq_h->Fill( v3_pi_eval1*v3_pi_eval1 );
        v4_sq_h->Fill( v4_pi_eval1*v4_pi_eval1 );
        v5_sq_h->Fill( v5_pi_eval1*v5_pi_eval1 );
        //cout<<"you got your phi from the fourier"<<endl;
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );  

        v1_pi_h->Fill( pT ,  0. );
        v2_pi_h->Fill( pT ,  0. );
        v3_pi_h->Fill( pT ,  0. );
        v4_pi_h->Fill( pT ,  0. );
        v5_pi_h->Fill( pT ,  0. ); 
        //cout<<"phi is random noise"<<endl;
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      piPlus -> SetPx( pT * TMath::Cos( phi ) );
      piPlus -> SetPy( pT * TMath::Sin( phi ) );
      piPlus -> SetPz( pT*  TMath::SinH( eta ) );
      piPlus -> SetMass(0.139570);
      piPlus -> SetKF(211);
      piPlus -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.139570*0.139570  ) );
      
    }
    
    for( Int_t N = yield_arr[CentralityBin][0] ; N < (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1]) ; N++  ){  //piMinus , kF = -211
    
      if( N== yield_arr[CentralityBin][0] ){
        //cout<<"piMinus getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;

      TMCParticle *piMinus = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,piMinus_params[CentralityBin][0]);
      pTdist->SetParameter(1,piMinus_params[CentralityBin][1]);
      pTdist->SetParameter(2,piMinus_params[CentralityBin][2]);
      pTdist->SetParameter(3,piMinus_params[CentralityBin][3]);
      pTdist->SetParameter(4,piMinus_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_piMinus->GetRandom();

/*_________________________________________________________________________________*/

      if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_pi_eval1 = 0;
      }

      else if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_pi_eval1 = (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_pi_eval1 = 0;
      }

      else if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_pi_eval1 = (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_pi_eval1 = 0;
      }

      else if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_pi_eval1 = (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_pi_eval1 = 0;
      }

      else if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_pi_eval1 = (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/
      //phi.clear();
      if (HarmonicsFlag != 5){
        Set_Phi_Dist(pT , -211 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for Pi- now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_pi_h->Fill( pT , v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_pi_h->Fill( pT , v2_pi_eval1 );
        v3_pi_h->Fill( pT , v3_pi_eval1 );
        v4_pi_h->Fill( pT , v4_pi_eval1 );
        v5_pi_h->Fill( pT , v5_pi_eval1 );

        v1_sq_h->Fill( v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_pi_eval1*v2_pi_eval1 );
        v3_sq_h->Fill( v3_pi_eval1*v3_pi_eval1 );
        v4_sq_h->Fill( v4_pi_eval1*v4_pi_eval1 );
        v5_sq_h->Fill( v5_pi_eval1*v5_pi_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_pi_h->Fill( pT ,  0. );
        v2_pi_h->Fill( pT ,  0. );
        v3_pi_h->Fill( pT ,  0. );
        v4_pi_h->Fill( pT ,  0. );
        v5_pi_h->Fill( pT ,  0. ); 
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      piMinus -> SetPx( pT * TMath::Cos( phi ) );
      piMinus -> SetPy( pT * TMath::Sin( phi ) );
      piMinus -> SetPz( pT*  TMath::SinH( eta ) );    
      piMinus -> SetMass(0.139570);
      piMinus -> SetKF(-211);
      piMinus -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.139570*0.139570  ) );

    }
    
    for( Int_t N = (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1]) ; N < (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2]) ; N++){ //kPlus , kF = 321
    
      if( N == (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1]) ){
        //cout<<"kPlus getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
    
      TMCParticle *k_Plus = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,kPlus_params[CentralityBin][0]);
      pTdist->SetParameter(1,kPlus_params[CentralityBin][1]);
      pTdist->SetParameter(2,kPlus_params[CentralityBin][2]);
      pTdist->SetParameter(3,kPlus_params[CentralityBin][3]);
      pTdist->SetParameter(4,kPlus_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_kPlus->GetRandom();
/*_________________________________________________________________________________*/

      if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_K_eval1 = 0;
      }

      else if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_K_eval1 = (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_K_eval1 = 0;
      }

      else if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_K_eval1 = (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_K_eval1 = 0;
      }

      else if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_K_eval1 = (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_K_eval1 = 0;
      }

      else if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_K_eval1 = (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/

      //phi.clear();
      if (HarmonicsFlag != 5){
        //Set_Phi_Dist(pT , 321 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for K+ now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_K_h->Fill( pT , v1_K->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_K_h->Fill( pT , v2_K_eval1 );
        v3_K_h->Fill( pT , v3_K_eval1 );
        v4_K_h->Fill( pT , v4_K_eval1 );
        v5_K_h->Fill( pT , v5_K_eval1 );

        v1_sq_h->Fill( v1_K->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_K->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_K_eval1*v2_K_eval1 );
        v3_sq_h->Fill( v3_K_eval1*v3_K_eval1 );
        v4_sq_h->Fill( v4_K_eval1*v4_K_eval1 );
        v5_sq_h->Fill( v5_K_eval1*v5_K_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_K_h->Fill( pT ,  0. );
        v2_K_h->Fill( pT ,  0. );
        v3_K_h->Fill( pT ,  0. );
        v4_K_h->Fill( pT ,  0. );
        v5_K_h->Fill( pT ,  0. ); 
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      k_Plus -> SetPx( pT * TMath::Cos( phi ) );
      k_Plus -> SetPy( pT * TMath::Sin( phi ) );
      k_Plus -> SetPz( pT*  TMath::SinH( eta ) );
      k_Plus -> SetMass(0.493677);
      k_Plus -> SetKF(321);
      k_Plus -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.493677*0.493677  ) );
    }
    
    for( Int_t N = (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2]) ; N < ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3]) ; N++ ){ //kMinus , kF = -321
    
      if( N == (yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2]) ){
        //cout<<"kMinus getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
    
      TMCParticle *kMinus = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,kMinus_params[CentralityBin][0]);
      pTdist->SetParameter(1,kMinus_params[CentralityBin][1]);
      pTdist->SetParameter(2,kMinus_params[CentralityBin][2]);
      pTdist->SetParameter(3,kMinus_params[CentralityBin][3]);
      pTdist->SetParameter(4,kMinus_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_kMinus->GetRandom();
/*_________________________________________________________________________________*/

      if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_K_eval1 = 0;
      }

      else if( (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_K_eval1 = (v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_K_eval1 = 0;
      }

      else if( (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_K_eval1 = (v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_K_eval1 = 0;
      }

      else if( (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_K_eval1 = (v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_K_eval1 = 0;
      }

      else if( (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_K_eval1 = (v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/

      //phi.clear();
      if (HarmonicsFlag != 5){
        Set_Phi_Dist(pT , -321 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for K- now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_K_h->Fill( pT , v1_K->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_K_h->Fill( pT , v2_K_eval1 );
        v3_K_h->Fill( pT , v3_K_eval1 );
        v4_K_h->Fill( pT , v4_K_eval1 );
        v5_K_h->Fill( pT , v5_K_eval1 );

        v1_sq_h->Fill( v1_K->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_K->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_K_eval1*v2_K_eval1 );
        v3_sq_h->Fill( v3_K_eval1*v3_K_eval1 );
        v4_sq_h->Fill( v4_K_eval1*v4_K_eval1 );
        v5_sq_h->Fill( v5_K_eval1*v5_K_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_K_h->Fill( pT ,  0. );
        v2_K_h->Fill( pT ,  0. );
        v3_K_h->Fill( pT ,  0. );
        v4_K_h->Fill( pT ,  0. );
        v5_K_h->Fill( pT ,  0. ); 
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      kMinus -> SetPx( pT * TMath::Cos( phi ) );
      kMinus -> SetPy( pT * TMath::Sin( phi ) );
      kMinus -> SetPz( pT*  TMath::SinH( eta ) );
      kMinus -> SetMass(0.493677);
      kMinus -> SetKF(-321);
      kMinus -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.493677*0.493677 ) );
    }
    
    for( Int_t N = ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3])  ; N < ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] ) ; N++ ){ //proton , kF = 2212
    
      if( N == ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3]) ){
        //cout<<"p getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
    
      TMCParticle *p = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,p_params[CentralityBin][0]);
      pTdist->SetParameter(1,p_params[CentralityBin][1]);
      pTdist->SetParameter(2,p_params[CentralityBin][2]);
      pTdist->SetParameter(3,p_params[CentralityBin][3]);
      pTdist->SetParameter(4,p_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_p->GetRandom();
/*_________________________________________________________________________________*/

      if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_P_eval1 = 0;
      }

      else if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_P_eval1 = (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_P_eval1 = 0;
      }

      else if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_P_eval1 = (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_P_eval1 = 0;
      }

      else if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_P_eval1 = (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_P_eval1 = 0;
      }

      else if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_P_eval1 = (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/ 
      //phi.clear();
      if(HarmonicsFlag != 5){
        Set_Phi_Dist(pT , 2212 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for P now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_P_h->Fill( pT , v1_P->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_P_h->Fill( pT , v2_P_eval1 );
        v3_P_h->Fill( pT , v3_P_eval1 );
        v4_P_h->Fill( pT , v4_P_eval1 );
        v5_P_h->Fill( pT , v5_P_eval1 );

        v1_sq_h->Fill( v1_P->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_P->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_P_eval1*v2_P_eval1 );
        v3_sq_h->Fill( v3_P_eval1*v3_P_eval1 );
        v4_sq_h->Fill( v4_P_eval1*v4_P_eval1 );
        v5_sq_h->Fill( v5_P_eval1*v5_P_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_P_h->Fill( pT ,  0. );
        v2_P_h->Fill( pT ,  0. );
        v3_P_h->Fill( pT ,  0. );
        v4_P_h->Fill( pT ,  0. );
        v5_P_h->Fill( pT ,  0. ); 
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      p -> SetPx( pT * TMath::Cos( phi ) );
      p -> SetPy( pT * TMath::Sin( phi ) );
      p -> SetPz( pT*  TMath::SinH( eta ) );
      p -> SetMass(0.938272);
      p -> SetKF(2212);
      p -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.938272*0.938272  ) );
    }
    
    for( Int_t N = ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] )  ; N < ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] + yield_arr[CentralityBin][5] ) ; N++ ){ //anti-proton , kF = -2212
    
      if( N == ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] ) ){
        //cout<<"pbar getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
 
      TMCParticle *pbar = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,pbar_params[CentralityBin][0]);
      pTdist->SetParameter(1,pbar_params[CentralityBin][1]);
      pTdist->SetParameter(2,pbar_params[CentralityBin][2]);
      pTdist->SetParameter(3,pbar_params[CentralityBin][3]);
      pTdist->SetParameter(4,pbar_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_pbar->GetRandom();
/*_________________________________________________________________________________*/

      if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_P_eval1 = 0;
      }

      else if( (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_P_eval1 = (v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_P_eval1 = 0;
      }

      else if( (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_P_eval1 = (v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_P_eval1 = 0;
      }

      else if( (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_P_eval1 = (v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_P_eval1 = 0;
      }

      else if( (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_P_eval1 = (v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

/*_________________________________________________________________________________*/ 
      //phi.clear();
      if(HarmonicsFlag != 5){
        Set_Phi_Dist(pT , -2212 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for Pbar now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_P_h->Fill( pT , v1_P->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_P_h->Fill( pT , v2_P_eval1 );
        v3_P_h->Fill( pT , v3_P_eval1 );
        v4_P_h->Fill( pT , v4_P_eval1 );
        v5_P_h->Fill( pT , v5_P_eval1 );

        v1_sq_h->Fill( v1_P->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_P->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_P_eval1*v2_P_eval1 );
        v3_sq_h->Fill( v3_P_eval1*v3_P_eval1 );
        v4_sq_h->Fill( v4_P_eval1*v4_P_eval1 );
        v5_sq_h->Fill( v5_P_eval1*v5_P_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_P_h->Fill( pT ,  0. );
        v2_P_h->Fill( pT ,  0. );
        v3_P_h->Fill( pT ,  0. );
        v4_P_h->Fill( pT ,  0. );
        v5_P_h->Fill( pT ,  0. );
      }
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t eta = etadist->GetRandom();
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      pbar -> SetPx( pT * TMath::Cos( phi ) );
      pbar -> SetPy( pT * TMath::Sin( phi ) );
      pbar -> SetPz( pT*  TMath::SinH( eta ) );
      pbar -> SetMass(0.938272);
      pbar -> SetKF(-2212);
      pbar -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.938272*0.938272 ) );
    
    }
    
    for( Int_t N = ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] + yield_arr[CentralityBin][5] )  ; N < Npart ; N++ ){ //piZer0 , kF = 111
    
          if( N == ( yield_arr[CentralityBin][0] + yield_arr[CentralityBin][1] + yield_arr[CentralityBin][2] + yield_arr[CentralityBin][3] + yield_arr[CentralityBin][4] + yield_arr[CentralityBin][5] ) ){
        //cout<<"piZero getting created now"<<endl;
      }
      
      //cout<<"Loop Iteration in .cxx  is at "<<N<<"th iteration"<<endl;
    
        
      TMCParticle *piZero = (TMCParticle*)BkgdArray->ConstructedAt(N);
      /*
      pTdist->SetParameter(0,piZero_params[CentralityBin][0]);
      pTdist->SetParameter(1,piZero_params[CentralityBin][1]);
      pTdist->SetParameter(2,piZero_params[CentralityBin][2]);
      pTdist->SetParameter(3,piZero_params[CentralityBin][3]);
      pTdist->SetParameter(4,piZero_params[CentralityBin][4]);
      */
      Double_t pT = pTdist_piZero->GetRandom();

/*_________________________________________________________________________________*/

      if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v2_pi_eval1 = 0;
      }

      else if( (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v2_pi_eval1 = (v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v3_pi_eval1 = 0;
      }

      else if( (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v3_pi_eval1 = (v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v4_pi_eval1 = 0;
      }

      else if( (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v4_pi_eval1 = (v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }

      if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) < 0 ){
        v5_pi_eval1 = 0;
      }

      else if( (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 )) >= 0 ){
        v5_pi_eval1 = (v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ));
      }


/*_________________________________________________________________________________*/
      //phi.clear();
      if(HarmonicsFlag != 5){
        Set_Phi_Dist(pT , 111 , Psi_1_event , Psi_3_event , Psi_5_event);
        //cout<<"\nGetting Phi for PiZero now\n"<<endl;
        phi = Harmonics_Phi_Dist->GetRandom();

        v1_pi_h->Fill( pT , v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0 ) );
        v2_pi_h->Fill( pT , v2_pi_eval1 );
        v3_pi_h->Fill( pT , v3_pi_eval1 );
        v4_pi_h->Fill( pT , v4_pi_eval1 );
        v5_pi_h->Fill( pT , v5_pi_eval1 );

        v1_sq_h->Fill( v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0)*v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0) );
        v2_sq_h->Fill( v2_pi_eval1*v2_pi_eval1 );
        v3_sq_h->Fill( v3_pi_eval1*v3_pi_eval1 );
        v4_sq_h->Fill( v4_pi_eval1*v4_pi_eval1 );
        v5_sq_h->Fill( v5_pi_eval1*v5_pi_eval1 );
      }
      else if(HarmonicsFlag == 5){
        phi = Harmonics_Phi_Dist_Rand->Uniform(0 , 2.0*TMath::Pi() );

        v1_pi_h->Fill( pT ,  0. );
        v2_pi_h->Fill( pT ,  0. );
        v3_pi_h->Fill( pT ,  0. );
        v4_pi_h->Fill( pT ,  0. );
        v5_pi_h->Fill( pT ,  0. ); 
      }
      //Double_t eta = etadist->GetRandom();
      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      //Double_t phi = phidist ->Uniform( 0. , 2.0*TMath::Pi() );
      piZero -> SetPx( pT * TMath::Cos( phi ) );
      piZero -> SetPy( pT * TMath::Sin( phi ) );
      piZero -> SetPz( pT*  TMath::SinH( eta ) );    
      piZero -> SetMass(0.139977);
      piZero -> SetKF(111);
      piZero -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + 0.139977*0.139977  ) );
    
    }    
  //bob->SetPx(gRandom->Uniform( 1.));
  //bob->SetPy(gRandom->Gaus( 0.5 , 0.2 ) );
  //bob->SetPz(gRandom->BreitWigner( 0.5 , 0.3 ) );
  //part -> SetPx( pXdist->GetRandom() );
  //Add part to TClonesArray


  if(PRINTQA == kFALSE ){ //empty out the histograms if you don't want to print them out, we don't want to take up too much memory
    
    v1_pi_h->Reset();
    v2_pi_h->Reset();
    v3_pi_h->Reset();
    v4_pi_h->Reset();
    v5_pi_h->Reset();

    v1_K_h->Reset();
    v2_K_h->Reset();
    v3_K_h->Reset();
    v4_K_h->Reset();
    v5_K_h->Reset();

    v1_P_h->Reset();
    v2_P_h->Reset();
    v3_P_h->Reset();
    v4_P_h->Reset();
    v5_P_h->Reset();

    v1_sq_h->Reset();
    v2_sq_h->Reset();
    v3_sq_h->Reset();
    v4_sq_h->Reset();
    v5_sq_h->Reset();

    Psi_1_h->Reset();
    Psi_2_h->Reset();
    Psi_3_h->Reset();
    Psi_4_h->Reset();
    Psi_5_h->Reset();

  }

  else if(PRINTQA == kTRUE ){ //void out the harmonics whic are 0
    //cout<<"did you step inside the reset IF ??"<<endl;
     //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)  // HF ---> HarmonicsFlag

    if( HarmonicsFlag == 1){ // v2-v5
      v1_pi_h->Reset();
      v1_K_h->Reset();
      v1_P_h->Reset();
      v1_sq_h->Reset();
    }
    else if(HarmonicsFlag == 2){ //v3-v5
      v1_pi_h->Reset();
      v2_pi_h->Reset();
      v1_K_h->Reset();
      v2_K_h->Reset();
      v1_P_h->Reset();
      v2_P_h->Reset();
      v1_sq_h->Reset();
      v2_sq_h->Reset();
    }
    else if(HarmonicsFlag == 3){ //v1-v4
      v5_pi_h->Reset();
      v5_K_h->Reset();
      v5_P_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 4){ //v1-v3
      v4_pi_h->Reset();
      v5_pi_h->Reset();
      v4_K_h->Reset();
      v5_K_h->Reset();
      v4_P_h->Reset();
      v5_P_h->Reset();
      v4_sq_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 6){ //v1-v2, v4-v5
      v3_pi_h->Reset();
      v3_K_h->Reset();
      v3_P_h->Reset();
      v3_sq_h->Reset();
    }
    else if(HarmonicsFlag == 7){ //v1-v3, v5
      v4_pi_h->Reset();
      v4_K_h->Reset();
      v4_P_h->Reset();
      v4_sq_h->Reset();
    }
    else if(HarmonicsFlag == 8){ //v1, v3-v5
      v2_pi_h->Reset();
      v2_K_h->Reset();
      v2_P_h->Reset();
      v2_sq_h->Reset();
    }
    else if(HarmonicsFlag == 9){ //v1 only
      v2_pi_h->Reset();
      v3_pi_h->Reset();
      v4_pi_h->Reset();
      v5_pi_h->Reset();
      v2_K_h->Reset();
      v3_K_h->Reset();
      v4_K_h->Reset();
      v5_K_h->Reset();
      v2_P_h->Reset();
      v3_P_h->Reset();
      v4_P_h->Reset();
      v5_P_h->Reset();  
      v2_sq_h->Reset();
      v3_sq_h->Reset();
      v4_sq_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 10){ //v2 only
      v1_pi_h->Reset();
      v3_pi_h->Reset();
      v4_pi_h->Reset();
      v5_pi_h->Reset();
      v1_K_h->Reset();
      v3_K_h->Reset();
      v4_K_h->Reset();
      v5_K_h->Reset();
      v1_P_h->Reset();
      v3_P_h->Reset();
      v4_P_h->Reset();
      v5_P_h->Reset();
      v1_sq_h->Reset();
      v3_sq_h->Reset();
      v4_sq_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 11){ //v3 only
      v1_pi_h->Reset();
      v2_pi_h->Reset();
      v4_pi_h->Reset();
      v5_pi_h->Reset();
      v1_K_h->Reset();
      v2_K_h->Reset();
      v4_K_h->Reset();
      v5_K_h->Reset();
      v1_P_h->Reset();
      v2_P_h->Reset();
      v4_P_h->Reset();
      v5_P_h->Reset();
      v1_sq_h->Reset();
      v2_sq_h->Reset();
      v4_sq_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 12){ //v4 only
      v1_pi_h->Reset();
      v2_pi_h->Reset();
      v3_pi_h->Reset();
      v5_pi_h->Reset();
      v1_K_h->Reset();
      v2_K_h->Reset();
      v3_K_h->Reset();
      v5_K_h->Reset();
      v1_P_h->Reset();
      v2_P_h->Reset();
      v3_P_h->Reset();
      v5_P_h->Reset();
      v1_sq_h->Reset();
      v2_sq_h->Reset();
      v3_sq_h->Reset();
      v5_sq_h->Reset();
    }
    else if(HarmonicsFlag == 13){ //v5 only
      v1_pi_h->Reset();
      v2_pi_h->Reset();
      v3_pi_h->Reset();
      v4_pi_h->Reset();
      v1_K_h->Reset();
      v2_K_h->Reset();
      v3_K_h->Reset();
      v4_K_h->Reset();
      v1_P_h->Reset();
      v2_P_h->Reset();
      v3_P_h->Reset();
      v4_P_h->Reset();
      v1_sq_h->Reset();
      v2_sq_h->Reset();
      v3_sq_h->Reset();
      v4_sq_h->Reset();
    }


  }

  cout<<"\n\nthe number of entries in the bkgd array per event = "<<BkgdArray->GetEntries()<<endl;
  //return BkgdArray;
  return BkgdArray;


}


Double_t TennGen::dNdPhi_VonNeumann(Double_t phiPart, Double_t pT, 
  Double_t Psi1 , Double_t Psi2, Double_t Psi3, Double_t Psi4, Double_t Psi5,
  Double_t v1,Double_t v2,Double_t v3,Double_t v4,Double_t v5){
    return (1.0+2.0*(v1*TMath::Cos(phiPart-Psi1) + v2*TMath::Cos(2.0*(phiPart-Psi2)) + v3*TMath::Cos(3.0*(phiPart-Psi3))+ v4*TMath::Cos(4.0*(phiPart-Psi4)) + v5*TMath::Cos(5.0*(phiPart-Psi5))));
}

TClonesArray *TennGen::GetBackground_TannerTest(){
  
  cout<<"Got Background"<<endl;

  const Double_t yield_arr[8][7] = { { (TMath::Ceil((EtaRange * 654)/0.5)) , (TMath::Ceil((EtaRange * 654)/0.5)) , (TMath::Ceil((EtaRange * 97)/0.5)) , (TMath::Ceil((EtaRange * 97)/0.5)) , (TMath::Ceil((EtaRange * 30)/0.5)) , (TMath::Ceil((EtaRange * 29)/0.5)) , (TMath::Ceil((EtaRange * 654)/0.5)) } , { (TMath::Ceil((EtaRange * 541)/0.5)) , (TMath::Ceil((EtaRange * 539)/0.5)) , (TMath::Ceil((EtaRange * 81)/0.5)) , (TMath::Ceil((EtaRange * 80)/0.5)) , (TMath::Ceil((EtaRange * 25)/0.5)) , (TMath::Ceil((EtaRange * 25)/0.5)) , (TMath::Ceil((EtaRange * 539)/0.5)) } , { (TMath::Ceil((EtaRange * 406)/0.5)) , (TMath::Ceil((EtaRange * 404)/0.5)) , (TMath::Ceil((EtaRange * 61)/0.5)) , (TMath::Ceil((EtaRange * 61)/0.5)) , (TMath::Ceil((EtaRange * 19)/0.5)) , (TMath::Ceil((EtaRange * 19)/0.5)) , (TMath::Ceil((EtaRange * 404)/0.5)) } , { (TMath::Ceil((EtaRange * 274)/0.5)) , (TMath::Ceil((EtaRange * 273 )/0.5)), (TMath::Ceil((EtaRange * 41)/0.5)) , (TMath::Ceil((EtaRange * 41)/0.5)) , (TMath::Ceil((EtaRange * 13)/0.5)) , (TMath::Ceil((EtaRange * 13)/0.5)) , (TMath::Ceil((EtaRange * 273)/0.5)) } , { (TMath::Ceil((EtaRange * 179)/0.5)) , (TMath::Ceil((EtaRange * 179)/0.5)) , (TMath::Ceil((EtaRange * 27)/0.5)) , (TMath::Ceil((EtaRange * 27)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 179)/0.5)) } , { (TMath::Ceil((EtaRange * 111)/0.5)) , (TMath::Ceil((EtaRange * 110)/0.5)) , (TMath::Ceil((EtaRange * 16)/0.5)) , (TMath::Ceil((EtaRange * 16)/0.5)) , (TMath::Ceil((EtaRange * 5)/0.5)) , (TMath::Ceil((EtaRange * 5)/0.5)) , 	(TMath::Ceil((EtaRange * 110)/0.5)) } , { (TMath::Ceil((EtaRange * 63)/0.5)) , (TMath::Ceil((EtaRange * 63)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 9)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 63)/0.5)) } , { (TMath::Ceil((EtaRange * 33)/0.5)) , (TMath::Ceil((EtaRange * 33)/0.5)) , ((TMath::Ceil(EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 4)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 2)/0.5)) , (TMath::Ceil((EtaRange * 33)/0.5)) } }; //setting multiplicities for particles

  Double_t Psi_1_event = Psi_1->Uniform( 0 , 2.0*TMath::Pi());
  cout << "\n\n\nPsi_1 from .cxx = "<< Psi_1_event <<" for Event " << endl;
  cout << "Random Psi_1 Seed from .cxx = "<< (Psi_1->GetSeed()) <<" for Event " << endl;
  Psi_1_h->Fill(Psi_1_event);
  Psi_2_h->Fill(0.);
  //Psi_1_event = (TMath::Pi()/2);
  Double_t Psi_3_event = Psi_3->Uniform( 0 , 2.0*TMath::Pi());
  cout << "Psi_3 from .cxx = "<< Psi_3_event <<" for Event " << endl;
  cout << "Random Psi_3 Seed from .cxx = "<< (Psi_3->GetSeed()) <<" for Event " << endl;
  Psi_3_h->Fill(Psi_3_event);
  Psi_4_h->Fill(0.);
  //Psi_3_event = (3*TMath::Pi()/2);
  Double_t Psi_5_event = Psi_5->Uniform( 0 , 2.0*TMath::Pi());
  cout << "Psi_5 from .cxx = "<< Psi_5_event <<" for Event " << endl;
  cout << "Random Psi_5 Seed from .cxx = "<< (Psi_5->GetSeed()) <<" for Event \n\n\n" << endl;
  Psi_5_h->Fill(Psi_5_event);
  //Psi_5_event = (5*TMath::Pi()/2);
  
  // Set even psi's 
  Double_t Psi_2_event = 0.0;
  Double_t Psi_4_event = 0.0;


  Int_t vN_bool_Array[5];
  if(HarmonicsFlag == 0){
		Int_t vN_bool_Array_temp[5] = {1,1,1,1,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 1){
		Int_t vN_bool_Array_temp[5] = {0,1,1,1,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 2){
		Int_t vN_bool_Array_temp[5] = {0,0,1,1,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 3){
		Int_t vN_bool_Array_temp[5] = {1,1,1,1,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 4){
		Int_t vN_bool_Array_temp[5] = {1,1,1,0,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 5){
		Int_t vN_bool_Array_temp[5] = {0,0,0,0,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 6){
		Int_t vN_bool_Array_temp[5] = {1,1,0,1,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 7){
		Int_t vN_bool_Array_temp[5] = {1,1,1,0,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 8){
		Int_t vN_bool_Array_temp[5] = {1,0,1,1,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 9){
		Int_t vN_bool_Array_temp[5] = {1,0,0,0,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 10){
		Int_t vN_bool_Array_temp[5] = {0,1,0,0,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 11){
		Int_t vN_bool_Array_temp[5] = {0,0,1,0,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 12){
		Int_t vN_bool_Array_temp[5] = {0,0,0,1,0};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}
  else if(HarmonicsFlag == 13){
		Int_t vN_bool_Array_temp[5] = {0,0,0,0,1};
		for(int j=0;j<4;j++){  vN_bool_Array[j] = vN_bool_Array_temp[j];}}


  Int_t Npart = 0;
  for (Int_t n = 0 ; n < 7 ; n++){
    Npart = Npart + yield_arr[CentralityBin][n];
  }
  //cout<<"Number of Total Particles for Centrality Bin "<<CentralityBin<<" = "<<Npart<<endl;
  cout<<"\n\n\n\n\nSeed From the .cxx = "<<Seed<<"\n\n\n\n\n"<<endl;
  Int_t Nparts_total=0;
  for(Int_t n=0; n<7; n++){
    cout << "Starting n=" <<n<<endl;
    for( Int_t N = 0 ; N < yield_arr[CentralityBin][n] ; N++){

      TMCParticle *newParticle = (TMCParticle*)BkgdArray->ConstructedAt(Nparts_total);
      Double_t pT;
      if(n == 0) pT = pTdist_piPlus->GetRandom();
      else if(n==1)  pT = pTdist_piMinus->GetRandom();
      else if(n==2)  pT = pTdist_kPlus->GetRandom();
      else if(n==3) pT = pTdist_kMinus->GetRandom();
      else if(n==4) pT = pTdist_p->GetRandom();
      else if(n==5) pT = pTdist_pbar->GetRandom();
      else if(n==6) pT = pTdist_piZero->GetRandom();

      Double_t v1_eval, v2_eval, v3_eval, v4_eval, v5_eval;

      

      if(n==0||n==1||n==6){ // get pion vN
        v1_eval = v1_pi->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v2_eval = v2_pi->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v3_eval = v3_pi->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v4_eval = v4_pi->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v5_eval = v5_pi->Eval( pT , 0.0 , 0.0 ,  0.0 );
      }
      else if(n==2||n==3){ // get kaon vN
        v1_eval = v1_K->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v2_eval = v2_K->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v3_eval = v3_K->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v4_eval = v4_K->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v5_eval = v5_K->Eval( pT , 0.0 , 0.0 ,  0.0 );
      }
      else if(n==4||n==5){ // get proton vN
        v1_eval = v1_P->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v2_eval = v2_P->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v3_eval = v3_P->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v4_eval = v4_P->Eval( pT , 0.0 , 0.0 ,  0.0 );
        v5_eval = v5_P->Eval( pT , 0.0 , 0.0 ,  0.0 );
      }

      if(v2_eval<0.0) v2_eval =0;
      if(v3_eval<0.0) v3_eval =0;
      if(v4_eval<0.0) v4_eval =0;
      if(v5_eval<0.0) v5_eval =0;


      // turn off harmonics which are not desired 
      Double_t v1 = vN_bool_Array[0]*v1_eval;
      Double_t v2 = vN_bool_Array[1]*v2_eval;
      Double_t v3 = vN_bool_Array[2]*v3_eval;
      Double_t v4 = vN_bool_Array[3]*v4_eval;
      Double_t v5 = vN_bool_Array[4]*v5_eval;

      Double_t Max_Phi = (1.0+2.0*(TMath::Abs(v1)+v2+v3+v4+v5));
      Int_t CHECK = 0.0;
     
      while(CHECK!=1){
         Double_t dndphi = gRandom->Uniform(0.0,Max_Phi);
         Double_t test_phi = gRandom->Uniform(0,2.0*TMath::Pi());
         Double_t VN_Value = dNdPhi_VonNeumann(test_phi,pT,Psi_1_event,Psi_2_event,Psi_3_event,Psi_4_event,Psi_5_event,v1,v2,v3,v4,v5);
         if(dndphi < VN_Value) {
          CHECK = 1;
          phi = test_phi;
         }
      }

      if(n==0||n==1||n==6){ // get pion vN
        v1_pi_h->Fill( pT , v1 );
        v2_pi_h->Fill( pT , v2 );
        v3_pi_h->Fill( pT , v3);
        v4_pi_h->Fill( pT , v4);
        v5_pi_h->Fill( pT , v5);
      }
      else if(n==2||n==3){ // get kaon vN
        v1_K_h->Fill( pT , v1 );
        v2_K_h->Fill( pT , v2 );
        v3_K_h->Fill( pT , v3);
        v4_K_h->Fill( pT , v4);
        v5_K_h->Fill( pT , v5);
      }
      else if(n==4||n==5){ // get proton vN
        v1_P_h->Fill( pT , v1 );
        v2_P_h->Fill( pT , v2 );
        v3_P_h->Fill( pT , v3);
        v4_P_h->Fill( pT , v4);
        v5_P_h->Fill( pT , v5);
      }

      v1_sq_h->Fill( v1*v1 );
      v2_sq_h->Fill( v2*v2 );
      v3_sq_h->Fill( v3*v3 );
      v4_sq_h->Fill( v4*v4 );
      v5_sq_h->Fill( v5*v5 );

      Double_t eta = etadist->Uniform(-0.5*(EtaRange/0.5) , 0.5*(EtaRange/0.5));
      Double_t particle_mass;
      Int_t particle_KF;   

      if(n == 0) particle_mass = 0.139570;
      else if(n==1) particle_mass = 0.139570;
      else if(n==2) particle_mass = 0.493677;
      else if(n==3) particle_mass = 0.493677;
      else if(n==4) particle_mass = 0.938272;
      else if(n==5) particle_mass = 0.938272;
      else if(n==6) particle_mass = 0.139977;

      if(n == 0) particle_KF = 211;
      else if(n==1) particle_KF = -211;
      else if(n==2) particle_KF = 321;
      else if(n==3) particle_KF = -321;
      else if(n==4) particle_KF = 2212;
      else if(n==5) particle_KF = -2212;
      else if(n==6) particle_KF = 111;

      newParticle -> SetPx( pT * TMath::Cos( phi ) );
      newParticle -> SetPy( pT * TMath::Sin( phi ) );
      newParticle -> SetPz( pT*  TMath::SinH( eta ) );
      newParticle -> SetMass(particle_mass);
      newParticle -> SetKF(particle_KF);
      newParticle -> SetEnergy( TMath::Sqrt( (pT * TMath::Cos( phi )*pT * TMath::Cos( phi )) + (pT * TMath::Sin( phi )*pT * TMath::Sin( phi )) + (pT*  TMath::SinH( eta )*pT*  TMath::SinH( eta )) + particle_mass*particle_mass  ) );
    

      Nparts_total++;
    }

  }

  if(PRINTQA == kFALSE ){ //empty out the histograms if you don't want to print them out, we don't want to take up too much memory
    
    v1_pi_h->Reset();
    v2_pi_h->Reset();
    v3_pi_h->Reset();
    v4_pi_h->Reset();
    v5_pi_h->Reset();

    v1_K_h->Reset();
    v2_K_h->Reset();
    v3_K_h->Reset();
    v4_K_h->Reset();
    v5_K_h->Reset();

    v1_P_h->Reset();
    v2_P_h->Reset();
    v3_P_h->Reset();
    v4_P_h->Reset();
    v5_P_h->Reset();

    v1_sq_h->Reset();
    v2_sq_h->Reset();
    v3_sq_h->Reset();
    v4_sq_h->Reset();
    v5_sq_h->Reset();

    Psi_1_h->Reset();
    Psi_2_h->Reset();
    Psi_3_h->Reset();
    Psi_4_h->Reset();
    Psi_5_h->Reset();

  }

  cout<<"\n\nthe number of entries in the bkgd array per event = "<<BkgdArray->GetEntries()<<endl;
  //return BkgdArray;
  return BkgdArray;


}