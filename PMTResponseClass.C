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


class PMTResponseObject {
  int n, numPEPeaks;
public:
  PMTResponseObject();
  //double Signal()(double *x, double *p){
  //return Signal(x,p);
  //}
  void setN(int i){
    n = i;
  }
  int getN(){
    return n;
  }

  double PhotoElectron_N(double *x, double *p);
  double Signal(double *x, double *p);  
  //double PhotoElectron_N()(double *x, double *p){
  //  return PhotoElectron_N(x,p);
  //}
    //double testFunction()(double *x, double*p){
  //  return n* testFunction(x,p);
  //}
};


/*
 * Constructor
 */
PMTResponseObject::PMTResponseObject(){
  n=1;
  numPEPeaks=12;
}

/*
 * Function For Signal, Sums Over PhotoElectron Series
 */
double PMTResponseObject::Signal(double *x, double *par){
  double signal = 0; // Initialize/Clear to 0;
  //cout << numPEPeaks;
  for (n=0; n < numPEPeaks; n++){    
    signal += PhotoElectron_N(x,par);
  }
  n = 0; //Set Back to Default
  return signal;
}

/*
 * Function for PhotoElectron peak N.
 */
double PMTResponseObject::PhotoElectron_N(double *x, double *par){
  Float_t xx = x[0];
  double Q_n, Sigma_n;
  Q_n = par[3]+n*par[1];
  n==0 ? Sigma_n =par[4] : Sigma_n = sqrt(n)*par[2]; 
  // cout << Sigma_n << endl;
  double f = TMath::PoissonI(n,par[0])*(1-par[5])*TMath::Gaus(xx,Q_n, Sigma_n, kTRUE);
  //double f = TMath::Gaus(xx, n*par[1], par[2], kTRUE);
  return f;
}

void DrawPMTResponse(){  
  
  // Draw the PMT Response Function
  PMTResponseObject* fObj = new PMTResponseObject();  
  //cout << fObj->getN();
  TF1 *f1 = new TF1("f1", fObj,&PMTResponseObject::Signal, 0, 300, 7, "PMTResponseObject", "Signal");
  setParameters(f1);
  f1->SetLineColor(kBlack);
  f1->Draw();    

  //Draw the Contribution from each of the Single PE 
  for(int i =0; i< 5; i++){
    PMTResponseObject* fObj = new PMTResponseObject();  //Class Object is strongly linked to TF1.
    fObj->setN(i);
    cout << fObj->getN() << endl;
    TF1 *f2 = new TF1(Form("f2_%i",i), fObj,&PMTResponseObject::PhotoElectron_N,
		      0, 300, 7, 
		      "PMTResponseObject", "PhotoElectron_N");
    setParameters(f2);
    f2->Draw("same");
  }
  //TF1 *f1 = new TF1("f2", PMTResponseObject::testFunction, 0, 10, 1);
}

void setParameters(TF1 *f){
  f->SetParameter(0,mu);
  f->SetParameter(1,SinglePE);
  f->SetParameter(2,SinglePESpread);
  f->SetParameter(3,Pedestal);
  f->SetParameter(4,PedestalSpread);
  f->SetParameter(5,ThermalFraction);
  f->SetParameter(6,ThermalParameter);
  f->SetNpx(1000); //Number of Points to Draw
}

