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
 *
 * The class is written, to be compatible with the a Root Fit function,
 * taking double *x and double *p.
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
  //}
  void setN(int i){
    n = i;
  }
  int getN(){
    return n;
  }

  
  double Response(double *x, double *p);
  double PhotoElectron_N(double *x, double *p);
  double Pedestal(double *x, double *p);
  double Background(double *x, double*p);
  //double Background_N(double *x, double *p);
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
  numPEPeaks=7;
}

/*
 * Function For Signal, Sums Over PhotoElectron Series
 */
double PMTResponseObject::Response(double *x, double *par){
  double signal = 0; // Initialize/Clear to 0;
  //cout << numPEPeaks;
  for (n=1; n < numPEPeaks; n++){    
    signal += PhotoElectron_N(x,par);
  }
  double pedestal = Pedestal(x, par);
  double background = Background(x, par);
  double response = signal + pedestal + background;
  return response;
}

/*
 * Function for PhotoElectron peak N.
 */
double PMTResponseObject::PhotoElectron_N(double *x, double *par){
  Float_t xx = x[0];
  double Q_n, Sigma_n;
  Q_n = par[3]+n*par[1] +par[5]/par[6];  //Last Part is Shift due to background
  Sigma_n = sqrt(n)*par[2]; 
  // cout << Sigma_n << endl;
  double f = TMath::PoissonI(n,par[0])*TMath::Gaus(xx,Q_n, Sigma_n, kTRUE);
  //double f = TMath::Gaus(xx, n*par[1], par[2], kTRUE);
  return f;
}
/*
 * Function for Pedestal
 */
double PMTResponseObject::Pedestal(double *x, double *par){
  Float_t xx = x[0];
  double pedestal = TMath::PoissonI(n,par[0])*(1-par[5])*TMath::Gaus(xx,par[3], par[4], kTRUE);
  return pedestal;
}
/*
 * Function for Thermal Background
 * Using Approximation Method
 */
double PMTResponseObject::Background(double *x, double *par){
  Float_t xx = x[0];
  double background = TMath::PoissonI(0,par[0])*par[5]*par[6]*TMath::Exp(-1*par[6]*(xx - par[3]));
  
  // Equation includes a heaviside function so:
  double heaviside = 0;
  if((xx - par[3]) < 0) heaviside = 0;
  if(xx - par[3] == 0) heaviside = 0.5;
  if(xx-par[3] > 0) heaviside =1;
  return heaviside*background;
}

void DrawPMTResponse(){  
  // Draw the PMT Response Function
  PMTResponseObject* fObj1 = new PMTResponseObject();  
  //cout << fObj->getN();
  TF1 *f1 = new TF1("f1", fObj1,&PMTResponseObject::Response, 0, 300, 7, "PMTResponseObject", "Response");
  setParameters(f1); //Helper Function Defined at end of the file
  f1->SetLineColor(kBlack);
  f1->Draw();    

  //Draw the Contribution from each of the Single PE 
  for(int i =1; i< 5; i++){
    PMTResponseObject* fObj2 = new PMTResponseObject();  //Class Object is strongly linked to TF1.
    fObj2->setN(i);
    cout << fObj2->getN() << endl;
    TF1 *f2 = new TF1(Form("f2_%i",i), fObj2,&PMTResponseObject::PhotoElectron_N,
		      0, 300, 7, 
		      "PMTResponseObject", "PhotoElectron_N");
    setParameters(f2);
    f2->Draw("same");
  }

  //Draw the Contribution from the Pedestal
  PMTResponseObject* fObj3 = new PMTResponseObject();  
  //cout << fObj->getN();
  TF1 *f3 = new TF1("f3", fObj3, &PMTResponseObject::Pedestal, 0, 300, 7, "PMTResponseObject", "Pedestal");
  setParameters(f3);
  f3->SetLineColor(kGreen);
  f3->Draw("same");


  //Draw the Contribution from Thermal Background
  PMTResponseObject* fObj4 = new PMTResponseObject();  
  //cout << fObj->getN();
  TF1 *f4 = new TF1("f4", fObj4, &PMTResponseObject::Background, 0, 300, 7, "PMTResponseObject", "Background");
  setParameters(f4);
  f4->SetLineColor(kBlue);
  f4->Draw("same");
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


/*
 * Code for Calculating the PMT Response, without currently using the
 * approximation method.  Currently commmented out.
 * NB Thermal Background is Convoluted with Signal.
 */

/*

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

double PMTResponseObject::Background(double *x, double *par){
  double background = 0;
  for (n=0; n < numPEPeaks; n++){    
    background += Background_N(x,par);
  }
  return background;
}

double PMTResponseObject::Background_N(double *x, double *par){
  Float_t xx = x[0];
  double background_N;
  double Q_n, Sigma_n, exp_n, erf1_n, sign_n, erf2_n;
  Q_n = par[3]+n*par[1];
  n==0 ? Sigma_n =par[4] : Sigma_n = sqrt(n)*par[2];
  //std::cout <<"Background_" << n << "Sigma_n: " << Sigma_n << endl;   
  std::cout <<"Background_" << n << endl; 
  
  exp_n = TMath::Exp(-1*par[6]*(xx-Q_n-par[6]*TMath::Power(Sigma_n,2)));
  erf1_n = TMath::Erf(TMath::Abs(par[3]-Q_n-par[6]*TMath::Power(Sigma_n,2))/(Sigma_n*sqrt(2)));
  sign_n = TMath::Sign(1.0, xx -Q_n - par[6]*TMath::Power(Sigma_n,2));
  erf2_n = TMath::Erf(TMath::Abs(xx-Q_n-par[6]*TMath::Power(Sigma_n,2))/(Sigma_n*sqrt(2)));

  background_N = TMath::PoissonI(n,par[0])*par[5]*par[6]/2.0 * exp_n*(erf1_n + sign_n*erf2_n);
  return background_N;
}
*/
