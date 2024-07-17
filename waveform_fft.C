using std::vector;


#include<fstream>
#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <functional>
#include <sys/stat.h>
#include <dirent.h>
#include <TROOT.h>
#include <TMath.h>
#include <TTree.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TThread.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>



void waveform_fft(){

  TFile *myFile = TFile::Open("Run3_EXFLU1_Compensated_W13_pad1.3_2.5E15neq_-30C/stats_Sr_Run3_250V_trig2750V.root");
  TTree *tree = dynamic_cast<TTree*>(myFile->Get("Analysis"));
  TTreeReader myReader("Analysis", myFile);
  TTreeReaderValue<std::vector<std::vector<double>>> myw(myReader, "w");
  TTreeReaderValue<std::vector<std::vector<double>>> myt(myReader, "t");



  // prompt event number (between 0 and 249999 for these)

  
  int event;
  std::cout << "event number?" << std::endl;
  scanf("%i",&event);
  if(event > myReader.GetEntries() || event < 0){
    std::cout << "please enter an integer between 0 and " << myReader.GetEntries() - 1 << "." << std::endl;
    std::cout << "\n event number?" << std::endl;
    scanf("%i",&event);
  }
  
  
  std::cout << "Analyzing event number: " << event << std::endl;

  myReader.SetEntry(event);
  std::cout << myReader.GetCurrentEntry() << std::endl;
   
  int n = myt->at(0).size();
  std::cout << "number of bins: " << n << std::endl;
  std::cout << "t[0]: " << myt->at(0).at(1) << std::endl;
  std::cout << "t[1]: " << myt->at(0).at(0) << std::endl;
  Double_t width = myt->at(0).at(1) - myt->at(0).at(0) ;
  Double_t fullwidth = myt->at(0).at(n-1) - myt->at(0).at(0) ;
  std::cout << "bin width: " << width << std::endl;
  std::cout << "full time window: " << fullwidth << std::endl;

  Double_t start = myt->at(0).at(0) - width/2;
  Double_t end = myt->at(0).at(n-1) + width/2;

  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;

  

  TCanvas *c1 = new TCanvas("c1","Waveform");
  TH1D *hwave = new TH1D("hwave", "hwave", n, start, end);
  hwave->SetTitle("Waveforms;time [s];amplitude [mV]");
  hwave->SetStats(0);
  hwave->SetFillColor(kBlue);
  hwave->SetLineColor(kBlue);
  hwave->SetFillStyle(3004);
    
  std::cout << "initialized hwave" << std::endl;


  for (Int_t i = 0; i < n; i++){
    
    hwave->Fill((myt->at(0).at(i) - width/2), (myw->at(0).at(i)/2));
        
			}

  // Analysis->Draw("w[0]:t[0]>>lwave","event[0]==7000","l*");

  // hwave->SetAxisRange(-0.01,0.025,"Y");
  hwave->Draw();

  TCanvas *c2 = new TCanvas("c2","FFT_Mag");
  TH1 *hfftm = nullptr;
  TVirtualFFT::SetTransform(nullptr);

  hfftm = hwave->FFT(hfftm, "MAG");


  TH1D* hfftm_rebin = new TH1D("Hist FFT","Histogram of FFT Magnitude", n-1,
			 hfftm->GetXaxis()->GetXmin() / fullwidth,
			 hfftm->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftm_rebin->SetBinContent(i, hfftm->GetBinContent(i));
  }

  hfftm_rebin->SetTitle("Magnitude of Fast Fourier Transform;frequency [Hz]");
  hfftm_rebin->SetStats(0);
  hfftm_rebin->SetFillColor(kBlue);
  hfftm_rebin->SetLineColor(kBlue);
  hfftm_rebin->SetFillStyle(3004);
  c2->cd();
  hfftm_rebin->Draw();
 
  

  
  TCanvas *c3 = new TCanvas("c3","FFT_Phase");
  TH1 *hfftp = nullptr;

  hfftp = hwave->FFT(hfftp, "PH");

  TH1D* hfftp_rebin = new TH1D("Hist FFT","Histogram of FFT Phase", n-1,
			 hfftp->GetXaxis()->GetXmin() / fullwidth,
			 hfftp->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftp_rebin->SetBinContent(i, hfftp->GetBinContent(i));
  }

  hfftp_rebin->SetTitle("Phase of Fast Fourier Transform;frequency [Hz]");
  hfftp_rebin->SetStats(0);
  hfftp_rebin->SetFillColor(kBlue);
  hfftp_rebin->SetLineColor(kBlue);
  hfftp_rebin->SetFillStyle(3004);
  c3->cd();
  hfftp_rebin->Draw();


  
}
