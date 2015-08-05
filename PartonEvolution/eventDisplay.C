//----------------------------------------------------
// Code to parse the partonic history files from AMPT
// and trace the scattering evolution of quarks prior
// to hadronization.
//
// Author: J. Orjuela-Koop
// 07-29-2015
//----------------------------------------------------
#include "TLorentzVector.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TPaletteAxis.h"

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

//Verbosity flag
bool verbosity = false;

//Categories for the number of scattering events. NSCAT = 0, 1, 2, >=3
const int NSCATT = 4;

//Vector of vectors to contain parton evolution for a given event
vector<vector<parton> > eventParticles;
vector<vector<TArrow*> > arrows;
vector<vector<TEllipse*> > circles;

//Output histograms
TH2F *hDisplay;

//Counter for the number of events processed
int evtnumber = 0;

//Parton radius for display
const float PARTON_RADIUS = 0.02;

//-----------------------------------
// Functions
//-----------------------------------

/*
 * Determines if a given parton occupies the same location as any other *scattered* parton
 */
bool partonOverlap(float x, float y)
{
  int count = 0;
  for(int i=0; i<eventParticles.size(); i++)
    {
      for(int j=0; j<eventParticles[i].size() && eventParticles[i].size()>1; j++)
	{
	  parton p = eventParticles[i][j];

	  if(p.x == x && p.y == y)
	    {
	      count++;
	    }
	}
    }

  if(count >=2)
    {
      return true;
    }

  return false;
}

/*
 * Prepare parton information for displaying by adusting positions of overlapping partons
 */
void initializeDisplay()
{
  //Displace partons occupying same location, just for ease of visualization
  for(int i=0; i<eventParticles.size(); i++)
    {
      for(int j=0; j<eventParticles[i].size() && eventParticles[i].size()>1; j++)
	{
	  parton p = eventParticles[i][j];

	  if(partonOverlap(p.x,p.y))
	    {
	      eventParticles[i][j].x = eventParticles[i][j].x+0.15*PARTON_RADIUS;
	      eventParticles[i][j].y = eventParticles[i][j].y+0.15*PARTON_RADIUS;
	    }
	}
    }
}

/*
 * Draw event display
 */
void draw()
{
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","Event Display",700,700);
  gPad->SetFillColor(kBlack);
  hDisplay->GetXaxis()->SetAxisColor(kWhite);
  hDisplay->GetYaxis()->SetAxisColor(kWhite);
  hDisplay->GetXaxis()->SetLabelColor(kWhite);
  hDisplay->GetYaxis()->SetLabelColor(kWhite);
  hDisplay->GetXaxis()->SetTitleColor(kWhite);
  hDisplay->GetYaxis()->SetTitleColor(kWhite);
  hDisplay->GetXaxis()->SetRangeUser(-0.2,2.6);
  hDisplay->GetYaxis()->SetRangeUser(-1.8,2.3);

  //Define gray COLZ scale for 'spectator' background
  UInt_t Number = 2;
  Double_t Red[2]   = { 0.00, 1.00};
  Double_t Green[2] = { 0.00, 1.00};
  Double_t Blue[2]  = { 0.00, 1.00};
  Double_t Stops[2] = { 0.00, 1.00};

  Int_t nb=50;
  TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);
  
  hDisplay->SetContour(nb);  
  
  hDisplay->Draw("COLZ");
  gPad->Update();
  TPaletteAxis *palette1 = (TPaletteAxis*) hDisplay->GetListOfFunctions()->FindObject("palette");
  palette1->SetLineColor(kWhite);
  palette1->SetLabelColor(kWhite);

  Color_t colorPartons[8] = {kViolet, kAzure-3, kGreen-3, kOrange-3, kRed, kGray, kPink+9, kSpring-9};
  
  int cindex = 0;
  for(int i=0; i<circles.size(); i++)
    {
      vector<TEllipse*> t = circles[i];
      for(int j=0; j<t.size(); j++)
	{
	  t[j]->SetLineColor(colorPartons[cindex]);
	  t[j]->SetFillColorAlpha(kWhite,0);
	  t[j]->Draw("same");
	}

      if(t.size() > 0) cindex++;
  }
  
  cindex = 0;
  for(int i=0; i<arrows.size(); i++)
    {
      vector<TArrow*> ta = arrows[i];
      for(int j=0; j<ta.size(); j++)
	{
	  ta[j]->SetLineWidth(2);
	  ta[j]->SetLineColor(colorPartons[cindex]);
	  //ta[j]->SetAngle(60);
	  ta[j]->Draw();
	}

      if(ta.size() > 0) cindex++;
    }
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
 * Compute azimuthal angle of a given parton
 */
float computePhi(parton p)
{
  float px = p.px;
  float py = p.py;
  float pz = p.pz;
  float mass = p.m;

  float energy =  sqrt (pow(px,2)+pow(py,2)+pow(pz,2)+pow(mass,2));
  TLorentzVector ev(px,py,pz,energy);

  return ev.Phi();
}

/*
 * Loop over partons in each event, filling desired histograms
 */
void processEvent()
{ 
  //Each entry in the 'eventParticles' vector corresponds to a single parton and its history
  for(int i=0; i<eventParticles.size(); i++)
    {   
      //Determine number of scattering events undergone by the parton
      vector<parton> v = eventParticles[i];
      int numstages = v.size();
      int numscatterings = v.size()-1;

      if(numscatterings > 2) numscatterings = 3;

      //Fill histograms
      float x_initial = v[0].x;
      float y_initial = v[0].y;
      float z_initial = v[0].z;

      //Fill in scattered parton positions
      if(numscatterings > -1)
	{
	  for(int j=0; j<numstages; j++)
	    {
	      hDisplay->Fill(v[j].x,v[j].y);
	    }
	}

      //Fill in arrows
      vector<TArrow*> partonArrows;
      for(int j=1; j<numstages && numstages>1; j++)
	{
	  TArrow *t = new TArrow(v[j-1].x, v[j-1].y, v[j].x, v[j].y,0.05,"->-");
	  partonArrows.push_back(t);
	}
      arrows.push_back(partonArrows);
      partonArrows.clear();

      //Fill in circles
      vector<TEllipse*> partonCircles;
      for(int j=0; j<numstages && numstages>1; j++)
	{
	  TEllipse *t = new TEllipse(v[j].x,v[j].y,PARTON_RADIUS,PARTON_RADIUS);
	  partonCircles.push_back(t);
	}
      circles.push_back(partonCircles);
      partonCircles.clear();

      //If verbosity is enabled, print evolution of every single parton
      if(verbosity)
	{
	  cout << endl << "-----------------------" << endl;
	  cout << "PARTON " << i << endl;
	  cout << "-----------------------" << endl << endl;
	  for(int j=0; j<numstages; j++)
	    {
	      parton p = v[j];
	      cout << "x_" << j << " = (" << p.x << ", " << p.y << ", " << p.z << ")" << endl;
	      cout << "p_" << j << " = (" << p.px << ", " << p.py << ", " << p.pz << ")" << endl << endl;
	    }
	}
    }
}

void eventDisplay(char *initialInfoFile = "parton-initial-afterPropagation.dat", char *evolInfoFile = "parton-collisionsHistory.dat", char *outputFile = "evolution_out.root")
{
  //Initialize event display
  hDisplay = new TH2F("hDisplay","Single AMPT Event;x;y",100,-3,3,100,-3,3);

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
      initializeDisplay();
      processEvent();
      eventParticles.clear();
    }

  draw();
}
