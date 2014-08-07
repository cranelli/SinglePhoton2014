#!/usr/bin/env python

#
# 1/24/2013 Basic code to plot a histogram for the single photon runs with
# the new Hamamatsu base.
#
# Can run from the command line by calling:
# 
# python PlotHist.py 


import string
import csv
from numpy import *
import matplotlib.pyplot as plt
import os
from math import *
from ROOT import TH1D,TF1, TSpectrum, TCanvas, TPaveText
from ROOT import TLegend
import sys

#Creates a histrogram from a file named by the user.

def hist(filename, dir ="", name="Single Photon Response", top="", ceiling=0.0):
    print filename

    #Depending on the Voltage of the Run, specify how to set up the histogram and fit.
    voltage = string.rsplit(filename,'_')[0]
    print voltage;
    if(voltage =="1500"):      
        
        mu =39.0
        sigma =11.0
        #Minimum Percent(of non-dark events), used to calculate num_peaks
        min_percent = 0.04   #Selected to have it only fit 1
        fit_range_min = 15.0
        floor = 0.1
        

    else: #Default is 835
        mu =0.4062
        sigma =0.1962
        min_percent = 0.01
        fit_range_min = 0.14
        floor = 0.08

    #Reading in our CSV file and place in array amps.
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        amps = [];
        for row in reader:
            for elem in row:
                amps.append(float(elem))

            
    nbins =150
    #nbins = Find_Nbin(amps)

    #h = TH1D("", "Histogram Title", 100, min(amps),max(amps))
    h = TH1D(top, name, nbins,0.5*((max(amps)-min(amps))/nbins), max(amps))
    h_ped = TH1D(top, "Pedestal", nbins,min(amps)+0.5*((max(amps)-min(amps))/nbins), max(amps))
    
    ceiling = Find_Ceiling(h,nbins)
    print "The ceiling is: " ,ceiling
    print "The maximum amplitude is: ",max(amps) 

    h.SetAxisRange(floor, ceiling)
    h_ped.SetAxisRange(min(amps), floor)
    
    zero_bin_count = 0.0
    total_events = len(amps)-1
    print "The total number of events is ", total_events
    bin_width =  h.GetBinWidth(1)
    
    for amp in amps:
        if amp < floor:
            zero_bin_count = zero_bin_count+1
            h_ped.Fill(amp)
        if amp > floor and amp < ceiling:
            h.Fill(amp)
    
    print "The percent of events in the zero bin is "
    zero_bin_per = zero_bin_count/total_events
    print  zero_bin_per

    
    print "The average number of photons, assuming a Poisson distribution" 
    mu_photons = -log(zero_bin_per)
    print mu_photons

    c2 = TCanvas()
    h_ped.Draw()

    c1=TCanvas()
    c1.cd()
    h.Draw()
        
    mean = h.GetMean()
    rms = h.GetRMS()
    pedestal = -1*abs(h_ped.GetMean()) #Some of the early runs were taking the absolute, instead of the negative.
    
    #Each additional peak requires ntimes this many events

    num_peaks = 1 #Doesn't includes the zero bin peak
    while(Pois(num_peaks+1,mu_photons)/(1-zero_bin_per)>min_percent):
        num_peaks = num_peaks+1

    #remember to remove this
    #num_peaks =7

    print num_peaks

    
    #num_peaks = 6 #Counts Zero Peak
    #These are the factors for each n-photon peak, assuming a poisson distribution.
    factor = [float(Pois(i,mu_photons)) for i in range (0,num_peaks+1)]

    norm = total_events*bin_width/sqrt(2*pi)
    print norm

    equation = str(norm)+"*("
    for i in range (1,num_peaks+1): #Our fit starts with the 1-photon peak   
        print "Factor", i ," = ", factor[i]
        if(i == 1):
            equation = equation + "(("+str(factor[i])+"*["+str(num_peaks+i)+"])/[1])*exp(-0.5*((x-([0]+"+str(pedestal)+"))/([1]))**2)"
        else :
            equation = equation + "+(("+str(factor[i])+"*["+str(num_peaks+i)+"])/([1]*[" + str(i)+"]))*exp(-0.5*((x-("+str(i)+"*[0])+"+str(pedestal)+")/([1]*["+str(i)+"]))**2)"

    equation = equation + ")"

    #print equation    
    #print "Total events ",total_events
    

    gfit_sum = TF1("Fit", equation,floor,ceiling);
    gfit_sum.SetRange(fit_range_min,(num_peaks+1)*mu)
    gfit_sum.SetParameter(0,mu);
    gfit_sum.SetParLimits(0,0.95*mu,1.05*mu)
    gfit_sum.SetParameter(1,sigma);
    gfit_sum.SetParLimits(1,0.95*sigma,1.05*sigma)
    #Set the ratio of sigma's
    for i in range (2, num_peaks+1):
        gfit_sum.SetParameter(i,sqrt(i))
        gfit_sum.SetParLimits(i,0.95,1.25*sqrt(i))#Should not be smaller than first sigma
    #Set the additional factor on the expected amplitude
    for i in range (num_peaks+1, 2*num_peaks+1):
        gfit_sum.SetParameter(i, 1.0)
        gfit_sum.SetParLimits(i,0.9,1.1)
        
    gfit_sum.SetLineColor(2)
    #gfit.SetMarkerColor(2)
    gfit_sum.SetLineWidth(2)
    
    h.Fit(gfit_sum, "LBNQR")
    
    gfit_sum.Draw("same")

    xmin = 0.69
    ymax = 0.75
    ymin = 0.25
    pave = TPaveText(xmin, ymin, 0.85,ymax, "NDC")
    
    muname = "#mu = %0.4f" %gfit_sum.GetParameter(0)+ " #pm %0.4f" %gfit_sum.GetParError(0)
    sigmaname = "#sigma_{1} = %0.4f" %gfit_sum.GetParameter(1)+ " #pm %0.4f" %gfit_sum.GetParError(1)
    
    
    pave.AddText(muname)
    pave.AddText(sigmaname)

    #Pave for Sigma Ratios
    for i in range(2, num_peaks+1):
        sigma_ratio_name = "#sigma_{"+str(i)+"}/#sigma_{1} = %0.2f" %gfit_sum.GetParameter(i) + " #pm %0.2f" %gfit_sum.GetParError(i)
        pave.AddText(sigma_ratio_name)
        
    
    #Pave for additional factor to the amplitude of each Gaussian
    #Aname = "A_{1} = %0.2f" %gfit_sum.GetParameter(num_peaks+1)+ " #pm %0.2f" %gfit_sum.GetParError(num_peaks+1) +" *A_{1 exp}"
    #pave.AddText(Aname)
    for i in range (1, num_peaks+1):
        Aname = "A_{"+str(i)+"} = %0.2f" %gfit_sum.GetParameter(num_peaks+i)+ " #pm %0.2f" %gfit_sum.GetParError(num_peaks+i) +" *A_{"+str(i)+" exp}"
        pave.AddText(Aname)

    #Pave for Number of Events within Each Gaussian
    for i in range (1, num_peaks+1):
        Nname = "N_{"+str(i)+"} = %i" %int(norm/bin_width*factor[i]*sqrt(2*pi)*gfit_sum.GetParameter(num_peaks+i))+ " #pm %i" %int(norm/bin_width*factor[i]*sqrt(2*pi)*gfit_sum.GetParError(num_peaks+i))
        #pave.AddText(Nname)

    
    pave.SetFillColor(0)
    pave.SetBorderSize(0)
    pave.SetTextFont(42)
    pave.SetTextSize(0.04)
    pave.Draw("same")


    for i in range (0,2*num_peaks+1):
        print gfit_sum.GetParameter(i)
        print gfit_sum.GetParError(i)

    
    #for i in range(1, numpeaks+1):
    #i =1;
    #gpeak = TF1("Peak", "gaus(0)", floor, ceiling)
    #gpeak.SetParameters(norm*factor[i]*gfit_sum.GetParameter(num_peaks+i)/gfit_sum.GetParameter(1), gfit_sum.GetParameter(0), gfit_sum.GetParameter(1))
    #gpeak.SetLineColor(3)
    #gpeak.SetLineWidth(2)
    #gpeak.Draw("same")

    gpeaks=[]
    for i in range(1, num_peaks+1):
        mean = (i*gfit_sum.GetParameter(0))+pedestal
        if(i==1):
            sigma_i =gfit_sum.GetParameter(1)
        else:
            sigma_i = gfit_sum.GetParameter(1)*gfit_sum.GetParameter(i)
        amp_gaus = norm*factor[i]*gfit_sum.GetParameter(num_peaks+i)/sigma_i
        gpeaks.append( TF1("Peak_%i" %i, "gaus(0)", floor, ceiling))
        gpeaks[-1].SetParameters(amp_gaus, mean, sigma_i)
    
        gpeaks[-1].SetLineColor(3)
        gpeaks[-1].SetLineWidth(2)
        gpeaks[-1].Draw("same")
    
    
    print "The total events in the fit are "
    print gfit_sum.Integral(floor, ceiling)/bin_width
    print "Expected total events are"
    print total_events*(1.0-zero_bin_per)
    
    oname = "../Plots_2013/" +voltage+"_%0.2f_fit_mu_sigma_bound_sum_mybase.gif" % zero_bin_per
    c1.SaveAs(oname)


    #Infinite Loop, so I can run from command line and keep Plots open.
    #When finished just kill the program
    while 1:
        pass

    return h
    
  #ceiling = Find_Ceiling(h,nbins)
    #ceiling = max(amps)
   # print "The ceiling is" ,ceiling

    #h.SetAxisRange(floor, ceiling)
    
    #zero_bin_count =0;
    
    #total_events = len(amps)-1
    #print "The total number of events is ", total_events
    #bin_width =  h.GetBinWidth(1)
    
    #for amp in amps:
        #if amp < floor:
            #zero_bin_count = zero_bin_count+1
        #if amp > floor and amp < ceiling:
        #    h.Fill(amp)
                
   # c1=TCanvas()
   # h.Draw()


def Pois(i , mu):
    return 1.0/factorial(i)*pow(mu,i)*exp(-1*mu)


def Most_Common(a):
    max_count = 0
    most_common_elem = 0
    for i in a:
        elem_count = 0
        for j in a:
            if i == j: elem_count= elem_count+1
            if(elem_count > max_count):
                max_count = elem_count
                most_common_elem = i
    return most_common_elem

    
#Uses the most common (non-zero) seperation between sorted elements to define
#the bin_width, and from there to get the number of bins.


def Find_Nbin(a):
    length = max(a)-min(a)
    a.sort()
    a = diff(a)
    a.sort()
    a = compress(a>1e-5,a) #Rounding Error (Don't Fully Understand)
    bin_width = Most_Common(a)
    nbins = int(length/bin_width)
    return nbins

def Find_Ceiling(h, nbins):
    for i in range(1, nbins):
        if(h.GetBinContent(nbins-i)!=0 and h.GetBinContent(nbins-i-1)!=0):
            return h.GetBinLowEdge(i)
        

    return h.GetBinLowEdge(nbins+1)

if __name__=="__main__":
    hist(sys.argv[1])

