//----------------------------------------------------
// Code to parse the partonic history files from AMPT
// and trace the scattering evolution of quarks prior
// to hadronization.
//
// Author: J. Orjuela-Koop
// 07-29-2015
//----------------------------------------------------
#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

using namespace std;

//-----------------------------------
// Parton Data Structure
//-----------------------------------
struct parton
{
  //Three-momentum
  float px;
  float py;
  float pz;

  //Four-position
  float x;
  float y;
  float z;
  float t;

  //Mass
  float m;
};

//-----------------------------------
// Variables
//-----------------------------------

//Categories for the number of scattering events. NSCAT = 0, 1, 2, >=3
const int NSCATT = 4;

//Vector of vectors to contain parton evolution for a given event
vector<vector<parton> > eventParticles;

//Output histograms
TProfile *hNscatt_pT;
TH1F *hdN_dpT[NSCATT];
TH1F *hDeltaRapidity[NSCATT];
TH1F *hRapidity[NSCATT];

//Counter for the number of events processed
int evtnumber = 0;

//-----------------------------------
// Functions
//-----------------------------------

/*
 * Draw spectra and v_2 for different number of scattering events
 */
void draw()
{
  gStyle->SetOptStat(0);

  TCanvas *cDeltaRapidity = new TCanvas("cDeltaRapidity","cDeltaRapidity",500,500);
  hDeltaRapidity[0]->Rebin(2);
  hDeltaRapidity[1]->Rebin(2);
  hDeltaRapidity[2]->Rebin(2);
  hDeltaRapidity[3]->Rebin(2);

  hDeltaRapidity[0]->SetLineColor(kViolet-3);
  hDeltaRapidity[1]->SetLineColor(kAzure-3);
  hDeltaRapidity[2]->SetLineColor(kSpring-6);
  hDeltaRapidity[3]->SetLineColor(kOrange-3);

  hDeltaRapidity[0]->Scale(1/(hDeltaRapidity[0]->GetMaximum()));
  hDeltaRapidity[1]->Scale(1/(hDeltaRapidity[1]->GetMaximum()));
  hDeltaRapidity[2]->Scale(1/(hDeltaRapidity[2]->GetMaximum()));
  hDeltaRapidity[3]->Scale(1/(hDeltaRapidity[3]->GetMaximum()));

  hDeltaRapidity[0]->SetTitle("");
  hDeltaRapidity[0]->Draw();
  hDeltaRapidity[1]->Draw("same");
  hDeltaRapidity[2]->Draw("same");
  hDeltaRapidity[3]->Draw("same");

  TCanvas *cRapidity = new TCanvas("cRapidity","cRapidity",500,500);
  hRapidity[0]->Rebin(2);
  hRapidity[1]->Rebin(2);
  hRapidity[2]->Rebin(2);
  hRapidity[3]->Rebin(2);

  hRapidity[0]->SetLineColor(kViolet-3);
  hRapidity[1]->SetLineColor(kAzure-3);
  hRapidity[2]->SetLineColor(kSpring-6);
  hRapidity[3]->SetLineColor(kOrange-3);

  hRapidity[0]->SetTitle("");
  hRapidity[0]->Draw();
  hRapidity[1]->Draw("same");
  hRapidity[2]->Draw("same");
  hRapidity[3]->Draw("same");

  TCanvas *cdNdpT = new TCanvas("cdNdpT","cdNdpT",500,500);
  hdN_dpT[0]->Rebin(2);
  hdN_dpT[1]->Rebin(2);
  hdN_dpT[2]->Rebin(2);
  hdN_dpT[3]->Rebin(2);

  hdN_dpT[0]->SetLineColor(kViolet-3);
  hdN_dpT[1]->SetLineColor(kAzure-3);
  hdN_dpT[2]->SetLineColor(kSpring-6);
  hdN_dpT[3]->SetLineColor(kOrange-3);

  hdN_dpT[0]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[1]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[2]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[3]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));

  hdN_dpT[0]->Scale(1.0/(float)evtnumber);
  hdN_dpT[1]->Scale(1.0/(float)evtnumber);
  hdN_dpT[2]->Scale(1.0/(float)evtnumber);
  hdN_dpT[3]->Scale(1.0/(float)evtnumber);

  hdN_dpT[0]->SetTitle("");
  hdN_dpT[0]->Draw();
  hdN_dpT[1]->Draw("same");
  hdN_dpT[2]->Draw("same");
  hdN_dpT[3]->Draw("same");

  TCanvas *cNscatt_pT = new TCanvas("cNscatt_pT","cNscatt_pT",500,500);
  hNscatt_pT->Rebin(2);
  hNscatt_pT->Draw();
}

/*
 * Compute the rapidity of a parton given its mass and momentum
 */
float computeRapidity(parton p)
{
  float px = p.px;
  float py = p.py;
  float pz = p.pz;
  float mass = p.m;

  float energy =  sqrt (pow(px,2)+pow(py,2)+pow(pz,2)+pow(mass,2));
  TLorentzVector ev(px,py,pz,energy);

  return ev.Rapidity();
}

/*
 * Loop over partons in each event, filling desired histograms
 */
void processEvent()
{
  //Each entry in the 'eventParticles' vector corresponds to a single parton and its history
  for(int i=0; i<eventParticles.size(); i++)
    {   
      vector<parton> v = eventParticles[i];
      int numstages = v.size();
      int numscatterings = v.size()-1;

      if(numscatterings > 2) numscatterings = 3;

      //Fill histograms
      float y = computeRapidity(v[numstages-1]);
      float delta_y = y - computeRapidity(v[0]);
      float pT = sqrt(pow(v[numstages-1].py,2) + pow(v[numstages-1].pz,2));

      hRapidity[numscatterings]->Fill(y);
      hDeltaRapidity[numscatterings]->Fill(delta_y);
      hdN_dpT[numscatterings]->Fill(pT);
      hNscatt_pT->Fill(pT,numscatterings,1);
    }
}

void parseHistoryFiles(char *initialInfoFile = "/direct/phenix+u/jdok/work/ampt/Ampt-v1.26t4-v2.26t4/ana/parton-initial-afterPropagation.dat", char *evolInfoFile = "/direct/phenix+u/jdok/work/ampt/Ampt-v1.26t4-v2.26t4/ana/parton-collisionsHistory.dat", char *outputFile = "evolution_out.root")
{
  //Initialize histograms
  for(int i=0; i<NSCATT; i++)
    {
      hDeltaRapidity[i] = new TH1F(Form("hDeltaRapidity_%i",i),Form("hDeltaRapidity_%i;#Delta y",i),500,-10,10);
      hRapidity[i] = new TH1F(Form("hRapidity_%i",i),Form("hRapidity_%i;y",i),500,-5,5);
      hdN_dpT[i] = new TH1F(Form("dN_dpT_%i",i),Form("dN_dpT_%i;p_{T};dN/dp_{T}",i),100,0,4);
    }

  hNscatt_pT = new TProfile("hNscatt_pT","Profile of Nscatt vs pT",100,0,5,0,8);

  //Read initial parton information file
  ifstream myFileInitialInfo;
  myFileInitialInfo.open(initialInfoFile);

  //Read parton evolution from file
  ifstream myFileEvolInfo;
  myFileEvolInfo.open(evolInfoFile);

  if (!myFileInitialInfo)
    {
      cout << "File does not exist!" << endl;
      return;
    }
  else
    {
      cout << "--> Successfully opened file " << initialInfoFile << endl;
    }

  while (myFileInitialInfo)
    {
      if(evtnumber % 100 == 0)
        {
          cout << "-----> // // Reading Event No. " << evtnumber << endl;
        }
      
      //Read the event header                                                                                                                                                
      int nlist;
      int iterindex;
      int nbaryons;
      int nmesons;
      int nparticles;
      int nparticleszpc;
      
      myFileInitialInfo >> evtnumber >> iterindex >> nlist >> nbaryons >> nmesons >> nparticles >> nparticleszpc;
      
      //Avoid double reading the last entry                                                                                                                                  
      if (!myFileInitialInfo) break;
      
      //Loop over each parton in the event                                                                                                                                   
      for (int i=0;i<nlist;i++)
        {
          int partid;
          float pvec[3];
          float mass;
          double spacetime[4];
	  
          myFileInitialInfo >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

	  parton part;
	  part.px = pvec[0];
	  part.py = pvec[1];
	  part.pz = pvec[2];
	  part.x = spacetime[0];
	  part.y = spacetime[1];
	  part.z = spacetime[2];
	  part.t = spacetime[3];
	  part.m = mass;

	  //Add the parton as the first stage of its own scattering evolution
	  vector<parton> aux;
	  aux.push_back(part);
	  eventParticles.push_back(aux);
	}
      
      //Get evolution information
      string line;
      TString tline;
      string description;
      int junk1;
      int junk2;
      int partonindex1;
      int partonindex2;
      int partid_init_1;
      float pvec_init_1[3];
      float mass_init_1;
      float spacetime_init_1[4];
      int partid_final_1;
      float pvec_final_1[3];
      float mass_final_1;
      float spacetime_final_1[4];
      int partid_init_2;
      float pvec_init_2[3];
      float mass_init_2;
      float spacetime_init_2[4];
      int partid_final_2;
      float pvec_final_2[3];
      float mass_final_2;
      float spacetime_final_2[4];
       
      string regexp = Form(" event,miss,iscat,jscat= %i [0-9]+ [0-9]+ [0-9]+", evtnumber);
      string regexp_break = Form("%i [0-9]+", evtnumber+1);
      while(myFileEvolInfo)
	{
	  TRegexp e(regexp);
	  TRegexp b(regexp_break);
	  std::getline(myFileEvolInfo,line);
	  tline = line;

	  if(tline.Contains(b))
	    {
	      break;
	    }

	  if(!tline.Contains(e))
	    {
	      continue;
	    }

	  stringstream headerstream(line);
	  headerstream >> description >> junk1 >> junk2 >> partonindex1 >> partonindex2;

	  myFileEvolInfo >> partid_init_1 >> pvec_init_1[0] >> pvec_init_1[1] >> pvec_init_1[2] >> mass_init_1 >> spacetime_init_1[0] >> spacetime_init_1[1] >> spacetime_init_1[2] >> spacetime_init_1[3];	
	  
	  myFileEvolInfo >> partid_init_2 >> pvec_init_2[0] >> pvec_init_2[1] >> pvec_init_2[2] >> mass_init_2 >> spacetime_init_2[0] >> spacetime_init_2[1] >> spacetime_init_2[2] >> spacetime_init_2[3];	 
	  
	  myFileEvolInfo >> partid_final_1 >> pvec_final_1[0] >> pvec_final_1[1] >> pvec_final_1[2] >> mass_final_1 >> spacetime_final_1[0] >> spacetime_final_1[1] >> spacetime_final_1[2] >> spacetime_final_1[3];
	  
	  myFileEvolInfo >> partid_final_2 >> pvec_final_2[0] >> pvec_final_2[1] >> pvec_final_2[2] >> mass_final_2 >> spacetime_final_2[0] >> spacetime_final_2[1] >> spacetime_final_2[2] >> spacetime_final_2[3];	 
	    
	  //Create scattered parton object from file information
	  parton part1;
	  part1.px = pvec_final_1[0];
	  part1.py = pvec_final_1[1];
	  part1.pz = pvec_final_1[2];
	  part1.x = spacetime_final_1[0];
	  part1.y = spacetime_final_1[1];
	  part1.z = spacetime_final_1[2];
	  part1.t = spacetime_final_1[3];
	  part1.m = mass_final_1;

	  parton part2;
	  part2.px = pvec_final_2[0];
	  part2.py = pvec_final_2[1];
	  part2.pz = pvec_final_2[2];
	  part2.x = spacetime_final_2[0];
	  part2.y = spacetime_final_2[1];
	  part2.z = spacetime_final_2[2];
	  part2.t = spacetime_final_2[3];
	  part2.m = mass_final_2;

	  eventParticles[partonindex1-1].push_back(part1);
	  eventParticles[partonindex2-1].push_back(part2);
	}

      //Run analysis on particles on an event-by-event basis
      processEvent();
      eventParticles.clear();
    }

  draw();
}
