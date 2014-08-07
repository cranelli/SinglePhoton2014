l//double fact(int n);
//double Pois(int i,double  mu);
//Double_t myfunction(Double_t *x, Double_t *par);
#include "TMath.h"
#include <iostream>

/*
 *For the inputs:
 *percent_0 is the number of 0 photon events.  Important for the Poisson distribution.
 *Mean, is the mean position of the single photon, and sigma is it's resolution.
 *The Gaussian distributions of the other peaks. are scaled to the first.
 */

void plotMultiGauss(){
  cout << "hello" << endl;
}


  // Double_t mu = -1.*log(percent_0);
  
  /*
  std::stringstream oname;
  oname <<"MultiPeak_" <<percent_0 <<"_" << mean <<"_"<<sigma <<"_"<<numpeaks <<".png";
  TCanvas *c = new TCanvas("c","",1200,800);
  c->cd();
  */
  /*
  TF1* f1 = new TF1("my func", myfunction, 0,(numpeaks +1)*mean, 4);
  f1->SetParameters(mu, mean, sigma, numpeaks);
  f1->SetLineColor(kRed);
  f1->Draw();
  
  for(int i =1; i < (numpeaks+1); i++){
    TF1* f2 = new TF1("f2", "gaus(0)", 0,(numpeaks+1)*mean);
    Double_t amp = Pois(i, mu);
    cout << amp << endl;
    f2->SetParameters(amp,i*mean,sigma);
    f2->SetLineWidth(2);
    f2->SetLineColor(kBlue);
    f2->Draw("same");
  }

  if(do_print) c->Print((oname.str().c_str()),"png");
  */

/*
Double_t myfunction(Double_t *x, Double_t *par){
  Float_t xx = x[0];
  Double_t f = 0;
  for(int i = 1; i< (par[3]+1); i++){
    //cout << Pois(i, par[0]) << endl;
    Double_t amp = Pois(i, par[0]);
    f +=  amp*exp(-0.5*((xx-i*par[1])/par[2])**2);
  }
  return f;
}


double fact(int n){
  if( n==1||n==0) return 1.0;
  return n*fact(n-1);
}


double Pois(int i, double mu){
  return (1.0/fact(i)*mu**i*exp(-1*mu));
}

*/
