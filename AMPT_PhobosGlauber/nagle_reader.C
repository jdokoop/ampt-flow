#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TLine.h"
#include "TProfile.h"

#include <iostream>
#include <fstream>

void nagle_reader(char *filename = "ana/ampt.dat") {

  gROOT->Reset();
  gROOT->SetStyle("Plain");

  ifstream myfile;
  myfile.open(filename);
  // add check to see if file was opened
  if (! myfile) {
    printf("No file with that name found\n");
    return;
  }

  TH1F *hphi = new TH1F("hphi","hphi",100,-TMath::Pi(),+TMath::Pi());

  TProfile *hv1mid = new TProfile("hv1mid","hv1mid",8,0.0,2.0);
  TProfile *hv1for = new TProfile("hv1for","hv1for",8,0.0,2.0);
  TProfile *hv1bac = new TProfile("hv1bac","hv1bac",8,0.0,2.0);

  int   nevent_perbin[10] = {0,0,0,0,0,0,0,0,0,0};
  TH1F *hdndeta[10];
  for (int i=0;i<10;i++) {
    char foobar[100];
    sprintf(foobar,"hdndeta%d",i);
    hdndeta[i] = new TH1F(foobar,foobar,48,-6.,6.);
  }
  int evtnumber = 0;
  while (myfile) {

    // read the event header
    int testnum;
    int nlist; // number of particles in out put list
    double impactpar;
    int npartproj;
    int nparttarg;
    int npartprojelas;
    int npartprojinelas;
    int nparttargelas;
    int nparttarginelas;
    myfile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas;

    cout << "Event number = " << evtnumber << " with b = " << impactpar << " [fm]" << endl;

    int impactbin = (int) (impactpar/2.0);
    if (impactbin<10 && impactbin >=0) nevent_perbin[impactbin]++;

    // store the particles....
    double temp_store[10000][2]; // px, py
    int temp_nstore = nlist;
    TLorentzVector *particle = new TLorentzVector();
    for (int i=0;i<nlist;i++) {
      int partid;
      double pvec[3];
      double mass;
      double spacetime[4];
      myfile >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> 
	spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];
      double energy =  sqrt (pow(pvec[0],2)+pow(pvec[1],2)+pow(pvec[2],2)+pow(mass,2));
      particle->SetPxPyPzE(pvec[0],pvec[1],pvec[2],energy);
      double pseudorapidity = 0.0;
      double phi = 0.0;
      if (pvec[0] == 0.0 && pvec[1] == 0.0) {
	pseudorapidity = -9999.0;
      } else {
	pseudorapidity = particle->Eta();
	phi = particle->Phi();
	if (impactbin < 10 && impactbin >=0) hdndeta[impactbin]->Fill(pseudorapidity);
	
	double pt = sqrt(pvec[0]*pvec[0]+pvec[1]*pvec[1]);
	if (pseudorapidity > 0.0 && pseudorapidity < 1.2)  hv1for->Fill(pt, pvec[0] / pt);
	if (pseudorapidity < 0.0 && pseudorapidity > -1.2) hv1bac->Fill(pt, pvec[0] / pt);
	if (pseudorapidity < 0.35 && pseudorapidity > -0.35) hv1mid->Fill(pt, pvec[0] / pt);

      }
      if (pseudorapidity > -1.0 && pseudorapidity < +1.0) {
	temp_store[i][0] = pvec[0];
	temp_store[i][1] = pvec[1];
	hphi->Fill(phi);
      }
    } // end loop over particles

    // calculations for the particles in this event...
    bool verbosity = true;

    float meanpartx = 0.0;
    float meanparty = 0.0;
    // calculate the mean values first
    for (int ii=0;ii<temp_nstore;ii++) {
      meanpartx += temp_store[ii][0];
      meanparty += temp_store[ii][1];
    }
    meanpartx = meanpartx / ((float) temp_nstore);
    meanparty = meanparty / ((float) temp_nstore);
    
    // then loop again and determine the moments
    // now that we have the mean, calculate the various eccentricity angles/moments
    // A = epsilon_N cos(N*psiN) = sum (r^2 cos(N*phi)) / sum(r^2)
    // B = epsilon_N sin(N*psiN) = sum (r^2 sin(N*phi)) / sum(r^2)
    // psiN = (atan2(B/A))/N
    // epsilon_N = A / cos(N*psiN)
    float A1 = 0.0;
    float B1 = 0.0;
    float A2 = 0.0;
    float B2 = 0.0;
    float A3 = 0.0;
    float B3 = 0.0;
    float A4 = 0.0;
    float B4 = 0.0;
    float A5 = 0.0;
    float B5 = 0.0;
    float sumr2 = 0.0;
    float impact = 0.0;
    for (int ii=0;ii<temp_nstore;ii++) {       
      float r2 = pow(temp_store[ii][0]-meanpartx,2) + pow(temp_store[ii][1]-meanparty,2);

      // atan2(y,x) returns angle -pi to +pi (with zero being straight up)
      /*      
      A1 += r2 * cos(1.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B1 += r2 * sin(1.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A2 += r2 * cos(2.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B2 += r2 * sin(2.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A3 += r2 * cos(3.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B3 += r2 * sin(3.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A4 += r2 * cos(4.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B4 += r2 * sin(4.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A5 += r2 * cos(5.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B5 += r2 * sin(5.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      */

      // should there be some momentum weighting (like a Q vector ???)

      A1 += cos(1.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B1 += sin(1.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A2 += cos(2.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B2 += sin(2.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A3 += cos(3.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B3 += sin(3.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A4 += cos(4.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B4 += sin(4.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      A5 += cos(5.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      B5 += sin(5.0*atan2(temp_store[ii][1]-meanparty,temp_store[ii][0]-meanpartx));
      
      //      sumr2 += r2;
      sumr2 += 1.0;
    }
    
    A1 = A1 / sumr2;
    B1 = B1 / sumr2;
    float psi1 = atan2(B1,A1)/2.0;
    float epsilon1 = A1 / cos(2.0*psi1);
    if (verbosity) cout << "Psi1 angle = " << psi1 << endl;
    if (verbosity) cout << "Epsilon1 = " << epsilon1 << endl;
    
    A2 = A2 / sumr2;
    B2 = B2 / sumr2;
    float psi2 = atan2(B2,A2)/2.0;
    float epsilon2 = A2 / cos(2.0*psi2);
    if (verbosity) cout << "Psi2 angle = " << psi2 << endl;
    if (verbosity) cout << "Epsilon2 = " << epsilon2 << endl;
    
    A3 = A3 / sumr2;
    B3 = B3 / sumr2;
    float psi3 = atan2(B3,A3)/3.0;
    float epsilon3 = A3 / cos(3.0*psi3);
    if (verbosity) cout << "Psi3 angle = " << psi3 << endl;
    if (verbosity) cout << "Epsilon3 = " << epsilon3 << endl;
    
    A4 = A4 / sumr2;
    B4 = B4 / sumr2;
    float psi4 = atan2(B4,A4)/4.0;
    float epsilon4 = A4 / cos(4.0*psi4);
    if (verbosity) cout << "Psi4 angle = " << psi4 << endl;
    if (verbosity) cout << "Epsilon4 = " << epsilon4 << endl;
    
    A5 = A5 / sumr2;
    B5 = B5 / sumr2;
    float psi5 = atan2(B5,A5)/5.0;
    float epsilon5 = A5 / cos(5.0*psi5);
    if (verbosity) cout << "Psi5 angle = " << psi5 << endl;
    if (verbosity) cout << "Epsilon5 = " << epsilon5 << endl;
    
    
  } // end loop over events

  hphi->SetLineColor(2);
  hphi->Draw();
  
  for (int i=0;i<10;i++) hdndeta[i]->Scale(1.0/ (double) nevent_perbin[i]);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  hdndeta[0]->SetXTitle("#eta");
  hdndeta[0]->SetYTitle("dN/d#eta");
  hdndeta[0]->SetMaximum(200.0);
  TLine *tl1 = new TLine(-2.2,0.0,-2.2,200.0);
  TLine *tl2 = new TLine(-1.2,0.0,-1.2,200.0);
  TLine *tl3 = new TLine( 0.0,0.0, 0.0,200.0);
  TLine *tl4 = new TLine( 1.2,0.0, 1.2,200.0);
  TLine *tl5 = new TLine( 2.2,0.0, 2.2,200.0);
  tl1->SetLineStyle(2);
  tl2->SetLineStyle(2);
  tl4->SetLineStyle(2);
  tl5->SetLineStyle(2);

  hdndeta[0]->DrawCopy("l");
  for (int i=0;i<6;i++) {
    hdndeta[i]->SetLineColor(i+1);
    hdndeta[i]->SetLineWidth(4);
    hdndeta[i]->Draw("l,same");
  }

  tl1->Draw("same");
  tl2->Draw("same");
  tl3->Draw("same");
  tl4->Draw("same");
  tl5->Draw("same");


  TLegend *tleg = new TLegend(0.6,0.7,0.9,0.9,"AMPT Cu+Au @ 200 GeV","brNDC");
  tleg->SetFillColor(0);
  for (int i=0;i<6;i++) {
    char foobar[100];
    sprintf(foobar,"b=%d-%d R(-(1.2-2.2)/+(1.2-2.2))=%3.2f",i*2,(i+1)*2,hdndeta[i]->Integral(16,19)/hdndeta[i]->Integral(30,33));
    tleg->AddEntry(hdndeta[i],foobar,"l");
  }
  tleg->Draw("same");

  TCanvas *c5 = new TCanvas("c5","c5",10,10,500,500);
  c5->cd();
  hv1mid->Scale(100.0);
  hv1bac->Scale(100.0);
  hv1for->Scale(100.0);
  hv1bac->SetMaximum(3.0);
  hv1bac->SetMinimum(-3.0);
  hv1bac->SetXTitle("p_{T} [GeV/c]");
  hv1bac->SetYTitle("v_{1} (%)");
  hv1bac->SetMarkerStyle(20);
  hv1bac->SetLineColor(1);
  hv1bac->DrawCopy("p,e,l");
  hv1mid->SetMarkerStyle(20);
  hv1mid->SetMarkerColor(2);
  hv1mid->SetLineColor(2);
  hv1mid->Draw("p,e,l,same");
  hv1for->SetMarkerStyle(20);
  hv1for->SetMarkerColor(4);
  hv1for->SetLineColor(4);
  hv1for->Draw("p,e,l,same");

  

}
