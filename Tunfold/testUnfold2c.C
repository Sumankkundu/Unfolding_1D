// Author: Stefan Schmitt
// DESY, July 2016

//  Version 17.9, example of using the SURE method

#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include "TUnfoldDensity.h"
#include "TUnfoldIterativeEM.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
// 
// Test program for the classes TUnfoldDensity and TUnfoldBinning
//
// The example comprizes several macros
//   testUnfold2a.C   create root files with TTree objects for
//                      signal, background and data
//            -> write files  testUnfold2a_MC.root
//                            testUnfold2a_data.root
//
//   testUnfold2b.C   loop over trees and fill histograms based on the
//                      TUnfoldBinning objects
//            -> read  testUnfold2a_MC.root
//                     testUnfold2a_data.root
//            -> write testUnfold2b_histograms.root
//            -> produce plots
//
//   testUnfold2c.C   run the unfolding
//            -> read  testUnfold2b_histograms.root
//            -> write testUnfold2c_unfolded.root
//            -> produce plots
// 
///////////////////////////////////////////////////////////////////////

static TH1 *generatePoissonToy(TH1 *base,int ntoy);
static void analyzeToy(TH1 const *hist_toy,
                       TH1 const *hist_truth,TProfile *&prof_pull,
                       TH1 *&coverage,int nToyTotal);

void testUnfold2c()
{
   //==================================================================  
   // switch on histogram errors
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
  
  //Input Data and MC histogram 
   TFile *inputData=new TFile("Data2017_GT102_3sigmaBin_16April.root");
   TFile *inputMC=new TFile("PY8_C5_1523_14PU_GT102X_3sig_16April.root");

  //Unfolder Data and Covarince matrix
   TFile *outputFile=new TFile("testunfold2c_unfolded.root","recreate");
   outputFile->cd();
   TDirectoryFile *inputDir=new TDirectoryFile("input","input");
   inputDir->cd();



   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetPadRightMargin(0.06);
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadLeftMargin(0.19);
   gStyle->SetPadTopMargin(0.07);
   int font=42;
   gStyle->SetLabelFont(font,"XYZ");
   gStyle->SetLegendFont(font);
   gStyle->SetStatFont(font);
   gStyle->SetTextFont(font);
   gStyle->SetTitleFont(font,"XYZ");
   gStyle->SetTitleOffset(1.5,"y");
   gStyle->SetTitleOffset(1.2,"x");
   gStyle->SetTitleSize(0.8,"P");
   gStyle->SetTitleSize(0.06,"xy");
   gStyle->SetLabelSize(0.05,"xy");
   gStyle->SetLabelOffset(0.012,"xy");


  int const type = 2;         // Jet & Charage particles
  int const itype[type]={0,1};   //{0}--->Jet ; {1}---> Charged Particles
  const  char* itypeN[type]={"Jets","Charged Particles"};
  char histname1[100] ,histname2[100], histname3[100], histname4[100], histname5[100], histname6[100], histname7[100],  histname8[100], histname9[100], histname10[100];
  char title1[100], title2[100], title3[100], title4[100];
  const int nHLTmx=8; //HT2 Range 
  const int njetetamn=1;  //eta value used 2.4
  double etarange[njetetamn] ={2.4}; //
  static const int nvar=32;  // Total number of eventshape variables
  static const int nusedvar = 5;   //Event Shape variables used
  Int_t var[nusedvar]={3,9,15,18,24};   // Names 3 Thrust , 9 Jet mass , 15 Y23 , 18 Jet Boardening , 24 Total Jet mass
  static const int nhist=10; //We need 8 But define 10
  const int njetptmn = nHLTmx;
  double leadingPtThreshold[njetptmn+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger
  
  const char* vartitle[nvar]={"Anti-Y_{23,C} ", "Anti-Y_{23,E} ", "Anti-Y_{23,R} ",     
                              "#tau_{_{#perp}} ", "#tau_{_{#perp} _{   ,E}} ", "#tau_{_{#perp} _{   ,R}} ",
                              "T_{ m,C} ", "T_{ m,E} ", "T_{ m,R} ",
                              "#rho_{Tot} ", "#rho_{Tot,E} ", "#rho_{Tot,R} ",
                              "#rho_{H,C} ", "#rho_{H,E} ", "#rho_{H,R} ",
                              "Y_{23} ", "Y_{23,E} ", "Y_{23,R} ",
                              "B_{ T,C} ", "B_{ T,E} ", "B_{ T,R} ",
                              "B_{ W,C} ", "B_{ W,E} ", "B_{ W,R} ",
                              "#rho^{T}_{Tot} ", "#rho^{T}_{Tot,E} ", "#rho^{T}_{Tot,R} ",
                              "#rho^{T}_{H,C} ", "#rho^{T}_{H,E} ", "#rho^{T}_{H,R} ",
                              "S_{_{#perp} _{   ,C}}", "C-parameter_{C}"};
  TH1F *PY8_Reco[type][nusedvar][njetptmn];  //Reconstracted MC
 //TH1F *PY8_Gen[type][nusedvar][njetptmn];   //Generator level MC
  TH1F *Data_Reco[type][nusedvar][njetptmn];    //Reconstracted Data
  TH2F *h2dGenDetMC[type][nusedvar][njetptmn];   // MC generator Vs Reco

   
//Read the MC and data
//Reco fine MC
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
        sprintf(histname1, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
        //      sprintf(histname2, "analyzeBasicPat/reco_typ_1_pt%i_eta0_%i", ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
        PY8_Reco[ity][ivar][ipt] = (TH1F*) inputMC->Get(histname1);
        //      PY8_Reco_Char[ivar][ipt] = (TH1F*) file1->Get(histname2);
        cout << histname1 << endl;
      }
    }
  }

/*  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
       for(int ipt = 0 ; ipt < njetptmn ; ipt++){
          sprintf(histname2, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
        PY8_Gen[ity][ivar][ipt] = (TH1F*) inputMC->Get(histname2);
      }
    }
  }
*/

  //copy the data
  //Reco_fine_Data
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
         sprintf(histname3, "analyzeBasicPat/reco_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //reco_typ_1_pt4_eta0_24
        Data_Reco[ity][ivar][ipt] = (TH1F*) inputData->Get(histname3);
      }
    }
  }
  //copy the correlation matrix
  //Reco fine coarse MC
  for(int ity=0; ity <type; ity++){
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){
        sprintf(histname4, "analyzeBasicPat/corr_typ_%i_pt%i_eta0_%i", ity, ipt, var[ivar]); //corr_typ_0_pt2_eta0_3
        h2dGenDetMC[ity][ivar][ipt] = (TH2F*) inputMC->Get(histname4);   //Xgen(coarse) , Yreco(fine)
      }
    }
  }


//Define Unfolded Histogram for outputs
  TH1 *Data_Unfolded[type][nusedvar][njetptmn];
 // TH1F *Data_Unfolded[type][nusedvar][njetptmn];
  TH2 *rhoij[type][nusedvar][njetptmn];

  

   
   inputFile->GetObject("CoarseGen",binningCoarseGen);
   inputFile->GetObject("FineReco",binningFineReco);


 for(int ity=0; ity <type; ity++){   //
    for(int ivar=0; ivar < nusedvar ; ivar ++){
      for(int ipt = 0 ; ipt < njetptmn ; ipt++){



   //========================================================
   // outer loop: unfold data, then unfold toys from unfolded data

   // histograms to draw toys from
   TH1 *hist_toybase_noRegularisation=0;
   TH1 *hist_toybase_TikhonovSURE=0;
   TH1 *hist_toybase_TikhonovLCurve=0;
   TH1 *hist_toybase_IterativeSURE=0;
   TH1 *hist_toybase_IterativeFixed=0;

   // analysis of toys
   TH1 *hist_unfolded_noRegularisation=0;
   TH1 *hist_unfolded_TikhonovSURE=0;
   TH1 *hist_unfolded_TikhonovLCurve=0;
   TH1 *hist_unfolded_IterativeSURE=0;
   TH1 *hist_unfolded_IterativeFixed=0;
   TProfile *prof_pull_noRegularisation=0;
   TProfile *prof_pull_TikhonovSURE=0;
   TProfile *prof_pull_TikhonovLCurve=0;
   TProfile *prof_pull_IterativeSURE=0;
   TProfile *prof_pull_IterativeFixed=0;
   TH1 *hist_coverage_noRegularisation=0;
   TH1 *hist_coverage_TikhonovSURE=0;
   TH1 *hist_coverage_TikhonovLCurve=0;
   TH1 *hist_coverage_IterativeSURE=0;
   TH1 *hist_coverage_IterativeFixed=0;

   // auxillary output
   // Tikhonov SURE scan
   TGraph *graph_logTauSURE_TikhonovSURE,*graph_dfChi2A_TikhonovSURE;
   int iBest_TikhonovSURE=-1;
   // Tikhonov, minimum L-Curve curvature
   TGraph *graph_LCurve_TikhonovLCurve;
   TSpline *spline_Curvature_TikhonovLCurve;
   int iBest_TikhonovLCurve=-1;
   double tauBest_TikhonovLCurve=-1.,DF_TikhonovLCurve=-1.;
   // Iterative SURE scan
   TGraph *graph_SURE_IterativeSURE,*graph_DFdeviance_IterativeSURE;
   int iBest_IterativeSURE=-1;

   double biasScale=1.0;

   int MAXTOY=300;

   int NPOINT_TikhonovSURE=50;
   int NPOINT_TikhonovLCurve=50;
   int NUM_FIXED_ITERATION=4;
   int NITER_Iterative=100;

   for(int itoy=0;itoy<=MAXTOY;itoy++) {
      cout<<"================== itoy="<<itoy<<" ==========================\n";
      TH1 *input;
      int iBest;
      //============================= (1) ============================
      // chi**2 fit without regularisation
      {
      if(itoy==0) {
         input=hist_unfoldingRecoFine_data;
      } else {
         input=generatePoissonToy
            (hist_toybase_noRegularisation,itoy);
      }

      TUnfoldDensity tunfoldNoRegularisation
         (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputHoriz,
          TUnfoldDensity::kRegModeSize,TUnfoldDensity::kEConstraintNone,
          TUnfoldDensity::kDensityModeNone,binningCoarseGen,
          binningFineReco);
      tunfoldNoRegularisation.SubtractBackground(hist_unfoldingRecoFine_bgr,"bgr");
      tunfoldNoRegularisation.SetInput(input,biasScale);
      tunfoldNoRegularisation.DoUnfold(0.0);
      //
      if(itoy==0) {
         // all unfolded bins (here including fakes normalisation)
         hist_unfolded_noRegularisation=
            tunfoldNoRegularisation.GetOutput
            ("hist_unfolded_noRegularisation",";bin",0,0,false);
         hist_unfolded_noRegularisation->Write();
         // unfolding result folded back -> for generating toys
         hist_toybase_noRegularisation=tunfoldNoRegularisation.GetFoldedOutput
            ("hist_toybase_noRegularisation",";bin",0,0,false,true);
         hist_toybase_noRegularisation->Write();
         // unfolding result, signal only
         TH1 *hist_PTunfolded_noRegularisation=
            tunfoldNoRegularisation.GetOutput
            ("hist_PTunfolded_noRegularisation",
             "P_{T,unfolded} [GeV]","signal");
         hist_PTunfolded_noRegularisation->Write();
      } else {
         // analyze toys
         TH1 *hist_unfolded_noRegularisation_toy=
            tunfoldNoRegularisation.GetOutput
            (TString::Format("hist_unfolded_noRegularisation_toy%d",itoy),
             ";bin",0,0,false);
         analyzeToy(hist_unfolded_noRegularisation_toy,
                    hist_unfolded_noRegularisation,
                    prof_pull_noRegularisation,
                    hist_coverage_noRegularisation,MAXTOY);
      }
      }
      //=========================================================
      // Tikhonov, minimum of the SURE variable
      {
      if(itoy==0) {
         input=hist_unfoldingRecoFine_data;
      } else {
         input=generatePoissonToy(hist_toybase_TikhonovSURE,itoy);
      }

      TUnfoldDensity tunfoldTikhonovSURE
         (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputHoriz,
          TUnfoldDensity::kRegModeSize,TUnfoldDensity::kEConstraintNone,
          TUnfoldDensity::kDensityModeNone,binningCoarseGen,
          binningFineReco);
      tunfoldTikhonovSURE.SubtractBackground(hist_unfoldingRecoFine_bgr,"bgr");
      tunfoldTikhonovSURE.SetInput(input,biasScale);
      iBest=tunfoldTikhonovSURE.ScanSURE
         (NPOINT_TikhonovSURE,0.,0.,
          (itoy==0) ? &graph_logTauSURE_TikhonovSURE : 0,
          (itoy==0) ? &graph_dfChi2A_TikhonovSURE : 0,
          0);
      if(itoy==0) {
         iBest_TikhonovSURE=iBest;
         // all unfolded bins (here including fakes normalisation)
         hist_unfolded_TikhonovSURE=
            tunfoldTikhonovSURE.GetOutput("hist_unfolded_TikhonovSURE",
                                          ";bin",0,0,false);
         hist_unfolded_TikhonovSURE->Write();
         // unfolding result folded back -> for generating toys
         hist_toybase_TikhonovSURE=tunfoldTikhonovSURE.GetFoldedOutput
            ("hist_toybase_TikhonovSURE",";bin",0,0,false,true);
         hist_toybase_TikhonovSURE->Write();
         // unfolding result, signal only
         TH1 *hist_PTunfolded_TikhonovSURE=
            tunfoldTikhonovSURE.GetOutput("hist_PTunfolded_TikhonovSURE",
                                          "P_{T,unfolded} [GeV]","signal");
         hist_PTunfolded_TikhonovSURE->Write();
         // save auxillary plots
         graph_logTauSURE_TikhonovSURE->Write("graph_logTauSURE_TikhonovSURE");
         graph_dfChi2A_TikhonovSURE->Write("graph_dfChi2A_TikhonovSURE");
      } else {
         // analyze toys
         TH1 *hist_unfolded_TikhonovSURE_toy=
            tunfoldTikhonovSURE.GetOutput
            (TString::Format("hist_PTunfolded_TikhonovSURE_toy%d",itoy),
             ";bin",0,0,false);
         analyzeToy(hist_unfolded_TikhonovSURE_toy,
                    hist_unfolded_TikhonovSURE,
                    prof_pull_TikhonovSURE,
                    hist_coverage_TikhonovSURE,MAXTOY);
      }
      }
      //=========================================================
      // Tikhonov, minimum L-Curve curvature
      {
      if(!hist_toybase_TikhonovLCurve) {
         input=hist_unfoldingRecoFine_data;
      } else {
         input=generatePoissonToy(hist_toybase_TikhonovLCurve,itoy);
      }
      TUnfoldDensity tunfoldTikhonovLCurve
         (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputHoriz,
          TUnfoldDensity::kRegModeSize,TUnfoldDensity::kEConstraintNone,
          TUnfoldDensity::kDensityModeNone,binningCoarseGen,
          binningFineReco);
      tunfoldTikhonovLCurve.SubtractBackground(hist_unfoldingRecoFine_bgr,"bgr");
      tunfoldTikhonovLCurve.SetInput(input,biasScale);
      iBest=tunfoldTikhonovLCurve.ScanLcurve
         (NPOINT_TikhonovLCurve,0.,0.,
          (itoy==0) ? &graph_LCurve_TikhonovLCurve : 0,
          0,0,
          (itoy==0) ? &spline_Curvature_TikhonovLCurve : 0);
      if(itoy==0) {
         iBest_TikhonovLCurve=iBest;
         tauBest_TikhonovLCurve=tunfoldTikhonovLCurve.GetTau();
         DF_TikhonovLCurve=tunfoldTikhonovLCurve.GetDF();
         // all unfolded bins (here including fakes normalisation)
         hist_unfolded_TikhonovLCurve=
            tunfoldTikhonovLCurve.GetOutput("hist_unfolded_TikhonovLCurve",
             ";bin",0,0,false);
         hist_unfolded_TikhonovLCurve->Write();
         // unfolding result folded back -> for generating toys
         hist_toybase_TikhonovLCurve=tunfoldTikhonovLCurve.GetFoldedOutput
            ("hist_toybase_TikhonovLCurve",";bin",0,0,false,true);
         // unfolding result, signal only
         TH1 *hist_PTunfolded_TikhonovLCurve=
            tunfoldTikhonovLCurve.GetOutput("hist_PTunfolded_TikhonovLCurve",
                                     "P_{T,unfolded} [GeV]","signal");
         hist_PTunfolded_TikhonovLCurve->Write();
         // save auxillary plots
         graph_LCurve_TikhonovLCurve->Write("graph_LCurve_TikhonovLCurve");
         spline_Curvature_TikhonovLCurve->Write
            ("spline_Curvature_TikhonovLCurve");
      } else {
         TH1 *hist_unfolded_TikhonovLCurve_toy=
            tunfoldTikhonovLCurve.GetOutput
            (TString::Format("hist_PTunfolded_TikhonovLCurve_toy%d",itoy),
             ";bin",0,0,false);
         analyzeToy(hist_unfolded_TikhonovLCurve_toy,
                    hist_unfolded_TikhonovLCurve,
                    prof_pull_TikhonovLCurve,
                    hist_coverage_TikhonovLCurve,MAXTOY);

      }
      }
      //============================================================
      // Iterative method, SURE scan
      {
      if(itoy==0) {
         input=hist_unfoldingRecoFine_data;
      } else {
         input=generatePoissonToy(hist_toybase_IterativeSURE,itoy);
      }
      TUnfoldIterativeEM tunfoldIterative
         (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputHoriz,
          binningCoarseGen,binningFineReco);
      tunfoldIterative.SubtractBackground(hist_unfoldingRecoFine_bgr,"bgr");
      tunfoldIterative.SetInput(input,biasScale);
      iBest=tunfoldIterative.ScanSURE
         (NITER_Iterative,
          (itoy==0) ? &graph_SURE_IterativeSURE : 0,
          (itoy==0) ? &graph_DFdeviance_IterativeSURE :0);
      if(itoy==0) {
         iBest_IterativeSURE=iBest;
         // all unfolded bins (here including fakes normalisation)
         hist_unfolded_IterativeSURE=
            tunfoldIterative.GetOutput("hist_unfolded_IterativeSURE",
                                       ";bin",0,0,false);
         hist_unfolded_IterativeSURE->Write();
         // unfolding result folded back -> for generating toys
         hist_toybase_IterativeSURE=tunfoldIterative.GetFoldedOutput
            ("hist_toybase_IterativeSURE",";bin",0,0,false,true);
         // unfolding result, signal only
         TH1 *hist_PTunfolded_IterativeSURE=
            tunfoldIterative.GetOutput("hist_PTunfolded_IterativeSURE",
                                       "P_{T,unfolded} [GeV]","signal");
         hist_PTunfolded_IterativeSURE->Write();
         // save auxillary plots
         graph_SURE_IterativeSURE->Write("graph_SURE_IterativeSURE");
         graph_DFdeviance_IterativeSURE->Write
            ("graph_DFdeviance_IterativeSURE");
      } else {
         TH1 *hist_unfolded_IterativeSURE_toy=
            tunfoldIterative.GetOutput
            (TString::Format("hist_unfolded_IterativeSURE_toy%d",itoy),
             ";bin",0,0,false);
         analyzeToy(hist_unfolded_IterativeSURE_toy,
                    hist_unfolded_IterativeSURE,
                    prof_pull_IterativeSURE,
                    hist_coverage_IterativeSURE,MAXTOY);
      }
      }
      //============================================================
      // Iterative method, four iterations
      {
      if(itoy==0) {
         input=hist_unfoldingRecoFine_data;
      } else {
         input=generatePoissonToy(hist_toybase_IterativeFixed,itoy);
      }
      TUnfoldIterativeEM tunfoldIterative
         (hist_migrationCoarseFine_MC,TUnfoldDensity::kHistMapOutputHoriz,
          binningCoarseGen,binningFineReco);
      tunfoldIterative.SubtractBackground(hist_unfoldingRecoFine_bgr,"bgr");
      tunfoldIterative.SetInput(input,biasScale);
      tunfoldIterative.DoUnfold(NUM_FIXED_ITERATION);
      if(itoy==0) {
         // all unfolded bins (here including fakes normalisation)
         hist_unfolded_IterativeFixed=
            tunfoldIterative.GetOutput("hist_unfolded_IterativeFixed",
                                       ";bin",0,0,false);
         hist_unfolded_IterativeFixed->Write();
         // unfolding result folded back -> for generating toys
         hist_toybase_IterativeFixed=tunfoldIterative.GetFoldedOutput
            ("hist_toybase_IterativeFixed",";bin",0,0,false,true);
         // unfolding result, signal only
         TH1 *hist_PTunfolded_IterativeFixed=
            tunfoldIterative.GetOutput("hist_PTunfolded_IterativeFixed",
                                       "P_{T,unfolded} [GeV]","signal");
         hist_PTunfolded_IterativeFixed->Write();
      } else {
         TH1 *hist_unfolded_IterativeFixed_toy=
            tunfoldIterative.GetOutput
            (TString::Format("hist_unfolded_IterativeFixed_toy%d",itoy),
             ";bin",0,0,false);
         analyzeToy(hist_unfolded_IterativeFixed_toy,
                    hist_unfolded_IterativeFixed,
                    prof_pull_IterativeFixed,
                    hist_coverage_IterativeFixed,MAXTOY);
      }
      }//end of iterative method
      
   
   }

   prof_pull_noRegularisation->Write();
   prof_pull_TikhonovSURE->Write();
   prof_pull_TikhonovLCurve->Write();
   prof_pull_IterativeSURE->Write();
   prof_pull_IterativeFixed->Write();

   hist_coverage_noRegularisation->Write();
   hist_coverage_TikhonovSURE->Write();
   hist_coverage_TikhonovLCurve->Write();
   hist_coverage_IterativeSURE->Write();
   hist_coverage_IterativeFixed->Write();

   // produce comparison plots
   // compare: (1) Tikhonov L-curve scan [L-curve]
   //          (2) Tikhonov SURE scan [SURE,DF,chi**2]
   //          (3) Iterative SURE scan [SURE,DF,deviance]

   
   delete outputFile;

}

