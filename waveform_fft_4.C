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



void waveform_fft_4(){


  TFile *myFile = TFile::Open("Run3_EXFLU1_Compensated_W13_pad1.3_2.5E15neq_-30C/stats_Sr_Run3_250V_trig2750V.root");
  TTree *tree = dynamic_cast<TTree*>(myFile->Get("Analysis"));
  TTreeReader myReader("Analysis", myFile);
  TTreeReaderValue<std::vector<std::vector<double>>> myw(myReader, "w");
  TTreeReaderValue<std::vector<std::vector<double>>> myt(myReader, "t");



  // SIGNAL WAVEFORMS (events 7111 and 7459)



  // SIGNAL WAVEFORM 1 (7111)
  
  int event = 7111;
  std::cout << "event number: " << event << std::endl;

  myReader.SetEntry(event);
  std::cout << myReader.GetCurrentEntry() << std::endl;
   
  int n = myt->at(0).size();
  std::cout << "number of bins: " << n << std::endl;
  std::cout << "t[0]: " << myt->at(0).at(1) << std::endl;
  std::cout << "t[1]: " << myt->at(0).at(0) << std::endl;
  Double_t width = myt->at(0).at(1) - myt->at(0).at(0) ;
  Double_t fullwidth = myt->at(0).at(n-1) - myt->at(0).at(0) ;
  std::cout << "bin width: " << width << std::endl;

  Double_t start = myt->at(0).at(0) - width/2;
  Double_t end = myt->at(0).at(n-1) + width/2;

  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;

  

  TCanvas *c1 = new TCanvas("c1","Waveform");
  TH1D *hwave1 = new TH1D("hwave1", "hwave1", n, start, end);
  hwave1->SetTitle("W13 Fluence 2.5E15, Bias 250V: Waveforms;time [s];amplitude [mV]");
  hwave1->SetStats(0);
  hwave1->SetFillColor(kBlue);
  hwave1->SetLineColor(kBlue);
  hwave1->SetFillStyle(3004);
    
  std::cout << "initialized hwave1" << std::endl;


  for (Int_t i = 0; i < n; i++){
    
    hwave1->Fill((myt->at(0).at(i) - width/2), (myw->at(0).at(i)/2));
        
			}

  // Analysis->Draw("w[0]:t[0]>>lwave","event[0]==7000","l*");

  hwave1->SetAxisRange(-0.01,0.025,"Y");
  hwave1->Draw();

  TCanvas *c2 = new TCanvas("c2","FFT_Mag");
  TH1 *hfftm1 = nullptr;
  TVirtualFFT::SetTransform(nullptr);

  hfftm1 = hwave1->FFT(hfftm1, "MAG");

  TH1D* hfftm_rebin1 = new TH1D("Hist FFT","Histogram of FFT Magnitude", n-1,
			 hfftm1->GetXaxis()->GetXmin() / fullwidth,
			 hfftm1->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftm_rebin1->SetBinContent(i, hfftm1->GetBinContent(i)/sqrt(n));
  }

  hfftm_rebin1->SetTitle("Magnitude of Fast Fourier Transform;frequency [Hz]");
  hfftm_rebin1->SetStats(0);
  hfftm_rebin1->SetFillColor(kBlue);
  hfftm_rebin1->SetLineColor(kBlue);
  hfftm_rebin1->SetFillStyle(3004);
  c2->cd();
  hfftm_rebin1->Draw("SAME");

  
  TCanvas *c3 = new TCanvas("c3","FFT_Phase");
  TH1 *hfftp1 = nullptr;

  hfftp1 = hwave1->FFT(hfftp1, "PH");
  TH1D* hfftp_rebin1 = new TH1D("Hist FFT","Histogram of FFT Phase", n-1,
			 hfftp1->GetXaxis()->GetXmin() / fullwidth,
			 hfftp1->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftp_rebin1->SetBinContent(i, hfftp1->GetBinContent(i));
  }

  hfftp_rebin1->SetTitle("Phase of Fast Fourier Transform;frequency [Hz]");
  hfftp_rebin1->SetStats(0);
  hfftp_rebin1->SetFillColor(kBlue);
  hfftp_rebin1->SetLineColor(kBlue);
  hfftp_rebin1->SetFillStyle(3004);
  c3->cd();
  hfftp_rebin1->Draw("SAME");





  // EVENT WAVEFORM 2


  event = 7459;
  std::cout << "event number: " << event << std::endl;

  myReader.SetEntry(event);
  std::cout << myReader.GetCurrentEntry() << std::endl;
   
  n = myt->at(0).size();
  std::cout << "number of bins: " << n << std::endl;
  std::cout << "t[0]: " << myt->at(0).at(1) << std::endl;
  std::cout << "t[1]: " << myt->at(0).at(0) << std::endl;
  width = myt->at(0).at(1) - myt->at(0).at(0) ;
  std::cout << "bin width: " << width << std::endl;

  start = myt->at(0).at(0) - width/2;
  end = myt->at(0).at(n-1) + width/2;

  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;

 
  TH1D *hwave2 = new TH1D("hwave2", "hwave2", n, start, end);
  hwave2->SetTitle("Waveforms;time [s];amplitude [mV]");
  hwave2->SetStats(0);
  hwave2->SetFillColor(kBlue);
  hwave2->SetLineColor(kBlue);
  hwave2->SetFillStyle(3005);
    
  std::cout << "initialized hwave2" << std::endl;


  for (Int_t i = 0; i < n; i++){
    
    hwave2->Fill((myt->at(0).at(i) - width/2), (myw->at(0).at(i)/2));
        
			}

  // Analysis->Draw("w[0]:t[0]>>lwave","event[0]==7000","l*");

  c1->cd();
  hwave2->Draw("SAME");

  TH1 *hfftm2 = nullptr;

  hfftm2 = hwave2->FFT(hfftm2, "MAG");


  TH1D* hfftm_rebin2 = new TH1D("Hist FFT","Histogram of FFT Magnitude", n-1,
			 hfftm2->GetXaxis()->GetXmin() / fullwidth,
			 hfftm2->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftm_rebin2->SetBinContent(i, hfftm2->GetBinContent(i)/sqrt(n));
    
  }

  hfftm_rebin2->SetTitle("Magnitude of Fast Fourier Transform;frequency [Hz]");
  hfftm_rebin2->SetStats(0);
  hfftm_rebin2->SetFillColor(kBlue);
  hfftm_rebin2->SetLineColor(kBlue);
  hfftm_rebin2->SetFillStyle(3005);
  c2->cd();
  hfftm_rebin2->Draw("SAME");
 
  
 
  TH1 *hfftp2 = nullptr;

  hfftp2 = hwave2->FFT(hfftp2, "PH");

  TH1D* hfftp_rebin2 = new TH1D("Hist FFT","Histogram of FFT Phase", n-1,
			 hfftp2->GetXaxis()->GetXmin() / fullwidth,
			 hfftp2->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftp_rebin2->SetBinContent(i, hfftp2->GetBinContent(i));
  }

  hfftp_rebin2->SetTitle("Phase of Fast Fourier Transform;frequency [Hz]");
  hfftp_rebin2->SetStats(0);
  hfftp_rebin2->SetFillColor(kBlue);
  hfftp_rebin2->SetLineColor(kBlue);
  hfftp_rebin2->SetFillStyle(3005);
  c3->cd();
  hfftp_rebin2->Draw("SAME");





  // NOISE WAVEFORMS (picked randomly:events 7000 and 20000)

  // NOISE WAVEFORM 1 (7000)



  event = 7000;
  std::cout << "event number: " << event << std::endl;

  myReader.SetEntry(event);
  std::cout << myReader.GetCurrentEntry() << std::endl;
   
  n = myt->at(0).size();
  std::cout << "number of bins: " << n << std::endl;
  std::cout << "t[0]: " << myt->at(0).at(1) << std::endl;
  std::cout << "t[1]: " << myt->at(0).at(0) << std::endl;
  width = myt->at(0).at(1) - myt->at(0).at(0) ;
  std::cout << "bin width: " << width << std::endl;

  start = myt->at(0).at(0) - width/2;
  end = myt->at(0).at(n-1) + width/2;

  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;

 
  TH1D *hwave3 = new TH1D("hwave3", "hwave3", n, start, end);
  hwave3->SetTitle("Waveforms;time [s];amplitude [mV]");
  hwave3->SetStats(0);
  hwave3->SetFillColor(kRed);
  hwave3->SetLineColor(kRed);
  hwave3->SetFillStyle(3006);
    
  std::cout << "initialized hwave3" << std::endl;


  for (Int_t i = 0; i < n; i++){
    
    hwave3->Fill((myt->at(0).at(i) - width/2), (myw->at(0).at(i)/2));
        
			}

  // Analysis->Draw("w[0]:t[0]>>lwave","event[0]==7000","l*");

  c1->cd();
  hwave3->Draw("SAME");

  TH1 *hfftm3 = nullptr;
  hfftm3 = hwave3->FFT(hfftm3, "MAG");


  TH1D* hfftm_rebin3 = new TH1D("Hist FFT","Histogram of FFT Magnitude", n-1,
			 hfftm3->GetXaxis()->GetXmin() / fullwidth,
			 hfftm3->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftm_rebin3->SetBinContent(i, hfftm3->GetBinContent(i)/sqrt(n));
  }

  hfftm_rebin3->SetTitle("Magnitude of Fast Fourier Transform;frequency [Hz]");
  hfftm_rebin3->SetStats(0);
  hfftm_rebin3->SetFillColor(kRed);
  hfftm_rebin3->SetLineColor(kRed);
  hfftm_rebin3->SetFillStyle(3006);
  c2->cd();
  hfftm_rebin3->Draw("SAME");
 
 
  TH1 *hfftp3 = nullptr;

  hfftp3 = hwave3->FFT(hfftp3, "PH");

  TH1D* hfftp_rebin3 = new TH1D("Hist FFT","Histogram of FFT Phase", n-1,
			 hfftp3->GetXaxis()->GetXmin() / fullwidth,
			 hfftp3->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftp_rebin3->SetBinContent(i, hfftp3->GetBinContent(i));
  }

  hfftp_rebin3->SetTitle("Phase of Fast Fourier Transform;frequency [Hz]");
  hfftp_rebin3->SetStats(0);
  hfftp_rebin3->SetFillColor(kRed);
  hfftp_rebin3->SetLineColor(kRed);
  hfftp_rebin3->SetFillStyle(3006);
  c3->cd();
  hfftp_rebin3->Draw("SAME");


  

  // NOISE WAVEFORM 2


  

  event = 20000;
  std::cout << "event number: " << event << std::endl;

  myReader.SetEntry(event);
  std::cout << myReader.GetCurrentEntry() << std::endl;
   
  n = myt->at(0).size();
  std::cout << "number of bins: " << n << std::endl;
  std::cout << "t[0]: " << myt->at(0).at(1) << std::endl;
  std::cout << "t[1]: " << myt->at(0).at(0) << std::endl;
  width = myt->at(0).at(1) - myt->at(0).at(0) ;
  std::cout << "bin width: " << width << std::endl;

  start = myt->at(0).at(0) - width/2;
  end = myt->at(0).at(n-1) + width/2;

  std::cout << "start: " << start << std::endl;
  std::cout << "end: " << end << std::endl;

 
  TH1D *hwave4 = new TH1D("hwave4", "hwave4", n, start, end);
  hwave4->SetTitle("W13 Fluence 2.5E15, Bias 250V: Waveforms;time [s];amplitude [mV]");
  hwave4->SetStats(0);
  hwave4->SetFillColor(kRed);
  hwave4->SetLineColor(kRed);
  hwave4->SetFillStyle(3007);
    
  std::cout << "initialized hwave4" << std::endl;


  for (Int_t i = 0; i < n; i++){
    
    hwave4->Fill((myt->at(0).at(i) - width/2), (myw->at(0).at(i)/2));
        
			}

  // Analysis->Draw("w[0]:t[0]>>lwave","event[0]==7000","l*");

  c1->cd();
  hwave4->Draw("SAME");

  TH1 *hfftm4 = nullptr;
  hfftm4 = hwave4->FFT(hfftm4, "MAG");


  TH1D* hfftm_rebin4 = new TH1D("Hist FFT","Histogram of FFT Magnitude", n-1,
			 hfftm4->GetXaxis()->GetXmin() / fullwidth,
			 hfftm4->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftm_rebin4->SetBinContent(i, hfftm4->GetBinContent(i)/sqrt(n));
  }

  hfftm_rebin4->SetTitle("Magnitude of Fast Fourier Transform;frequency [Hz]");
  hfftm_rebin4->SetStats(0);
  hfftm_rebin4->SetFillColor(kRed);
  hfftm_rebin4->SetLineColor(kRed);
  hfftm_rebin4->SetFillStyle(3007);
  c2->cd();
  hfftm_rebin4->Draw("SAME");
 
  

  TH1 *hfftp4 = nullptr;

  hfftp4 = hwave4->FFT(hfftp4, "PH");

  TH1D* hfftp_rebin4 = new TH1D("Hist FFT","Histogram of FFT Phase", n-1,
			 hfftp4->GetXaxis()->GetXmin() / fullwidth,
			 hfftp4->GetXaxis()->GetXmax() / fullwidth );

  for(int i = 0; i < n; i++)
  {
    hfftp_rebin4->SetBinContent(i, hfftp4->GetBinContent(i));
  }

  hfftp_rebin4->SetTitle("Phase of Fast Fourier Transform;frequency [Hz]");
  hfftp_rebin4->SetStats(0);
  hfftp_rebin4->SetFillColor(kRed);
  hfftp_rebin4->SetLineColor(kRed);
  hfftp_rebin4->SetFillStyle(3007);
  c3->cd();
  hfftp_rebin4->Draw("SAME");

  
  auto* legend1 = new TLegend(.6,.6,.9,.9);
  legend1->AddEntry(hwave1,"Signal Waveform (event 7111)");
  legend1->AddEntry(hwave2,"Signal Waveform (event 7459)");
  legend1->AddEntry(hwave3,"Noise Waveform (event 7000)");
  legend1->AddEntry(hwave4,"Noise Waveform (event 20000)");

  c1->cd();
  legend1->Draw();
  c2->cd();
  legend1->Draw();
  c3->cd();
  legend1->Draw();




  
}
