#include "TMath.h"
#include "TF1.h"
#include <iostream>

/*
 * This code Defines a mathematical Function describing
 * the response of a PhotoMultiplier Tube.
 * It is most effective when the average number of photons
 * expected is small, and can be described a Poisson distribution.
 * It also accounts for background contributionsfrom the electronic 
 * pedestal and thermal noise. 
 *
 * While it is useable over a full range, the fit is most useful,
 * when the 0 bin pedestal is still visible, mu <=6, and multiple
 * PEs contribute to the response, mu >=0.15.
 * Authors claim to measure the gain to within 1% accuracy.
 */

const double NumberOfPulses=10000; //Number of Laser Pulses in Sample
const double mu=1.68; //Mean PhotoElectrons for Poisson Distribution
const double SinglePE=35.04;// Position of Single PhotoElectron
const double SinglePESpread=11.73; // Spread of Single PhotoElectron
const double Pedestal=23.26 ; // Position of Pedestal
const double PedestalSpread=0.192;//Spread of Pedestal
const double ThermalFraction=0.383; //Fraction of Thermal Noise
const double ThermalParameter=0.034; //Thermal Noise Parameter

void DrawPMTResponse(){
  TF1 *f1 = new TF1("f1",pmtResponseFunction,0,10,1);
  f1->SetParameter(0,0.5);
  f1->Draw();
  //string PMTResponse = DefinePMTResponseFunction();
  //cout << PMTResponse << endl;
  //test();
}


double pmtResponseFunction(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Double_t f = signalFunction(x,par);
  return f;
}

double signalFunction(Double_t *x, Double_t *par){
  for(int n =0; n <12; n++){
    f += signalFunction_N(Double_t *x, Double_t *par, n);
  }
  return f;
}

double signalFunction_N(Double_t *x, Double_t *par){
  Float xx = x[0];
  Double_t f = SignalFunction(x,par):

    pois_i = "TMath::PoissonI(" + i_str + ",[0])";
    Q_i = "[3]+" + i_str + "*[1]"; 
    if(i==0) Sigma_i = "[4]";
    if(i!=0) Sigma_i = "sqrt(" + i_str + ")*[2]";

    gaus_i = "TMath::Gaus(x," + Q_i + "," + Sigma_i + ")/(sqrt(TMath::TwoPi())*"+Sigma_i+")"; //Have Included Normalization
    
    signal_i = "(1-[5])*" + gaus_i;

}

/*
 * Defines the PMT response function.
 * Different descriptions of the thermal background,
 * are defined depending on the average number of 
 * photons, mu.
 */
string DefinePMTResponseFunction(){
  
  string function;
  string pois_i, signal_i, background_i, gaus_i, exp_i, sign_i, erf1_i, erf2_i,  Q_i, Sigma_i;
  for(int i =0; i <7; i++){
    pois_i.clear();
    signal_i.clear();
    background_i.clear();
    gaus_i.clear();
    sign_i.clear();
    erf1_i.clear();
    erf2_i.clear();
    Q_i.clear();
    Sigma_i.clear();
    stringstream i_stream; 
    i_stream << i;
    string i_str = i_stream.str();
    if(i!=0){ 
      function += "+";
    }
    
    pois_i = "TMath::PoissonI(" + i_str + ",[0])";
    Q_i = "[3]+" + i_str + "*[1]"; 
    if(i==0) Sigma_i = "[4]";
    if(i!=0) Sigma_i = "sqrt(" + i_str + ")*[2]";

    gaus_i = "TMath::Gaus(x," + Q_i + "," + Sigma_i + ")/(sqrt(TMath::TwoPi())*"+Sigma_i+")"; //Have Included Normalization
    
    signal_i = "(1-[5])*" + gaus_i;
    
    exp_i = "TMath::Exp(-1*[6]*(x -" + Q_i + "-[6]*TMath::Power(" + Sigma_i +",2)))";
    erf1_i = "TMath::Erf(TMath::Abs([3]-" + Q_i + "-[6]*TMath::Power(" + Sigma_i+",2))/("+Sigma_i+"*sqrt(2))" ")";
    sign_i = "TMath::Sign(1, x-" + Q_i + "-[6]*TMath::Power(" + Sigma_i +",2))";
    erf2_i = "TMath::Erf(TMath::Abs(x-" + Q_i + "-[6]*TMath::Power(" + Sigma_i+",2))/("+Sigma_i+"*sqrt(2))" ")";
    
    background_i = "[5]*[6]/2.0 *" + exp_i + "*(" + erf1_i + "+" + sign_i + "*" +erf2_i +")";
    // cout << "working" << endl;
    function += pois_i + "*(" + signal_i + "+" + background_i + ")";
    
  }
  return function;
}


test(){
  //TMath math;
  cout << TMath::PoissonI(1, 0.1) << endl;
  cout << TMath::Gaus(1,2,1,kTRUE) << endl;

  //TF1 *f1 = new TF1("f1", "TMath::PoissonI(x/[1], [0])", -1, 10);
  
  string function, background;
  string pois_i, signal_i, background_i, gaus_i, exp_i, sign_i, erf1_i, erf2_i,  Q_i, Sigma_i;
  for(int i =0; i <7; i++){
    pois_i.clear();
    signal_i.clear();
    background_i.clear();
    gaus_i.clear();
    sign_i.clear();
    erf1_i.clear();
    erf2_i.clear();
    Q_i.clear();
    Sigma_i.clear();
    stringstream i_stream; 
    i_stream << i;
    string i_str = i_stream.str();
    if(i!=0){ 
      function += "+";
      background += "+";
    }
    
    pois_i = "TMath::PoissonI(" + i_str + ",[0])";
    Q_i = "[3]+" + i_str + "*[1]"; 
    if(i==0) Sigma_i = "[4]";
    if(i!=0) Sigma_i = "sqrt(" + i_str + ")*[2]";

    gaus_i = "TMath::Gaus(x," + Q_i + "," + Sigma_i + ")/(sqrt(TMath::TwoPi())*"+Sigma_i+")"; //Have Included Normalization

    signal_i = "(1-[5])*" + gaus_i;

    exp_i = "TMath::Exp(-1*[6]*(x -" + Q_i + "-[6]*TMath::Power(" + Sigma_i +",2)))";
    erf1_i = "TMath::Erf(TMath::Abs([3]-" + Q_i + "-[6]*TMath::Power(" + Sigma_i+",2))/("+Sigma_i+"*sqrt(2))" ")";
    sign_i = "TMath::Sign(1, x-" + Q_i + "-[6]*TMath::Power(" + Sigma_i +",2))";
    erf2_i = "TMath::Erf(TMath::Abs(x-" + Q_i + "-[6]*TMath::Power(" + Sigma_i+",2))/("+Sigma_i+"*sqrt(2))" ")";
    //cout << exp_i << endl;
    //cout << erf1_i << endl;
    //cout << erf2_i << endl;
    //cout << sign_i << endl;
    
    background_i = "[5]*[6]/2.0 *" + exp_i + "*(" + erf1_i + "+" + sign_i + "*" +erf2_i +")";
    // cout << "working" << endl;
    function += pois_i + "*(" + signal_i + "+" + background_i + ")";
    
    background += pois_i+ "*" + background_i;
  }

  cout << "Function: " << function << endl;
  TF1 *f1 = new TF1("f1",function.c_str(), 0 ,300); 
  f1->SetNpx(1000);
  f1->SetParameter(0,mu); //Mean of Poisson Distribution
  f1->SetParameter(1,SinglePE); // Position of Single PhotoElectron
  f1->SetParameter(2,SinglePESpread); // Spread of Single PhotoElectron
  f1->SetParameter(3,Pedestal); // Position of Pedestal
  f1->SetParameter(4,PedestalSpread); //Spread of Pedestal
  f1->SetParameter(5,ThermalFraction); //Fraction of Thermal Noise
  f1->SetParameter(6, ThermalParameter); //Thermal Noise Parameter
  f1->SetLineColor(kBlack);
  f1->Draw();
  
  //Background

  cout << "Background: " << background << endl;
  //For the Background
  TF2 *f2 = new TF1("f2",background.c_str(), 0 ,300); 
  f2->SetNumberFitPoints(600);
  f2->SetParameter(0,mu); //Mean of Poisson Distribution
  f2->SetParameter(1,35.04); // Position of Single PhotoElectron
  f2->SetParameter(2,11.73); // Spread of Single PhotoElectron
  f2->SetParameter(3,23.26); // Position of Pedestal
  f2->SetParameter(4,0.192); //Spread of Pedestal
  f2->SetParameter(5,0.383); //Fraction of Thermal Noise
  f2->SetParameter(6, 0.034); //Thermal Noise Parameter
  f2->SetLineColor(kBlue);
  f2->Draw("same");

  
  // Photoelectron Peaks and Pedestal

  for(int i =0; i <6; i++){
    pois_i.clear();
    signal_i.clear();
    gaus_i.clear();
    Q_i.clear();
    Sigma_i.clear();
    stringstream i_stream; 
    i_stream << i;
    string i_str = i_stream.str();
    if(i!=0){ 
      function += "+";
      background += "+";
    }
    
    pois_i = "TMath::PoissonI(" + i_str + ",[0])";
    Q_i = "[3]+" + i_str + "*[1]"; 
    if(i==0) Sigma_i = "[4]";
    if(i!=0) Sigma_i = "sqrt(" + i_str + ")*[2]";

    gaus_i = "TMath::Gaus(x," + Q_i + "," + Sigma_i + ")/(sqrt(TMath::TwoPi())*"+Sigma_i+")"; //Have Included Normalization

    signal_i = pois_i + "*" + "(1-[5])*" + gaus_i;
    
    TF3 *f3 = new TF1("f3",signal_i.c_str(), 0 ,600); 
    f3->SetNumberFitPoints(600);
    f3->SetParameter(0,mu); //Mean of Poisson Distribution
    f3->SetParameter(1,35.04); // Position of Single PhotoElectron
    f3->SetParameter(2,11.73); // Spread of Single PhotoElectron
    f3->SetParameter(3,23.26); // Position of Pedestal
    f3->SetParameter(4,0.192); //Spread of Pedestal
    f3->SetParameter(5,0.383); //Fraction of Thermal Noise
    f3->SetParameter(6, 0.034); //Thermal Noise Parameter
    f3->SetLineColor(kRed);
    f3->Draw("same"); 
		     
  }
}
