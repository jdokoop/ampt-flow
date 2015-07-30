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

//Categories for the number of scattering events. NSCAT = 0, 1, 2, 3, >3
const int NSCATT = 5;

//Vector of vectors to contain parton evolution for a given event
vector<vector<parton> > eventParticles;

//Output histograms
TH1F *hDeltaRapidity[NSCATT];
TH1F *hRapidity[NSCATT];

//Counter for the number of events processed
int evtnumber = 0;

//-----------------------------------
// Functions
//-----------------------------------

void draw()
{
  gStyle->SetOptStat(0);

  TCanvas *cDeltaRapidity = new TCanvas("cDeltaRapidity","cDeltaRapidity",500,500);
  hDeltaRapidity[0]->SetLineColor(kViolet-3);
  hDeltaRapidity[1]->SetLineColor(kAzure-3);
  hDeltaRapidity[2]->SetLineColor(kSpring-6);
  hDeltaRapidity[3]->SetLineColor(kOrange-3);
  hDeltaRapidity[4]->SetLineColor(kRed);

  hDeltaRapidity[0]->Draw();
  hDeltaRapidity[1]->Draw("same");
  hDeltaRapidity[2]->Draw("same");
  hDeltaRapidity[3]->Draw("same");
  hDeltaRapidity[4]->Draw("same");

  TCanvas *cRapidity = new TCanvas("cRapidity","cRapidity",500,500);
  hRapidity[0]->SetLineColor(kViolet-3);
  hRapidity[1]->SetLineColor(kAzure-3);
  hRapidity[2]->SetLineColor(kSpring-6);
  hRapidity[3]->SetLineColor(kOrange-3);
  hRapidity[4]->SetLineColor(kRed);

  hRapidity[0]->Draw();
  hRapidity[1]->Draw("same");
  hRapidity[2]->Draw("same");
  hRapidity[3]->Draw("same");
  hRapidity[4]->Draw("same");
}

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

void processEvent()
{
  //Each entry in the 'eventParticles' vector corresponds to a single parton and its history
  for(int i=0; i<eventParticles.size(); i++)
    {   
      vector<parton> v = eventParticles[i];
      int numstages = v.size();
      int numscatterings = v.size()-1;

      if(numscatterings > 3) numscatterings = 4;

      float y = computeRapidity(v[numstages-1]);
      float delta_y = y - computeRapidity(v[0]);

      hRapidity[numscatterings]->Fill(y);
      hDeltaRapidity[numscatterings]->Fill(delta_y);

      //Iterate over successive scatterings of the ith parton in the event
	/*
      cout << endl << "--------------" << endl;
      cout << "PARTON " << i+1 << endl;
      cout << "--------------" << endl << endl;

      for(int j=0; j<v.size(); j++)
	{
	  parton p = v[j];
	  cout << "   p_" << j << " = ("<< p.px << ", " << p.py << ", " << p.pz <<")" << endl;
	}
	*/
    }
}

void parseHistoryFiles(char *initialInfoFile = "parton-initial-afterPropagation.dat", char *evolInfoFile = "parton-collisionsHistory.dat", char *outputFile = "evolution_out.root")
{
  //Initialize histograms
  for(int i=0; i<NSCATT; i++)
    {
      hDeltaRapidity[i] = new TH1F(Form("hDeltaRapidity_%i",i),Form("hDeltaRapidity_%i",i),500,-10,10);
      hRapidity[i] = new TH1F(Form("hRapidity_%i",i),Form("hRapidity_%i",i),500,-5,5);
    }

  //Read initial parton information file
  ifstream myFileInitialInfo;
  myFileInitialInfo.open(initialInfoFile);

  //Read parton evolution information from file
  ifstream myFileEvolInfo;
  myFileEvolInfo.open(evolInfoFile);

  if (!myFileInitialInfo || !myFileEvolInfo)
    {
      cout << "File does not exist!" << endl;
      return;
    }
  else
    {
      cout << "--> Successfully opened file " << initialInfoFile << endl;
      cout << "--> Successfully opened file " << evolInfoFile << endl << endl;
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
      while(myFileEvolInfo)
	{
	  TRegexp e(regexp);
	  std::getline(myFileEvolInfo,line);
	  tline = line;

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
  cout << "Hola Inmundo!" << endl;
}
