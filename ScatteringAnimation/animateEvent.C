#include "TLorentzVector.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPad.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TLatex.h"
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

//Participant plane for event
float psi2;

//Vector of vectors to contain parton evolution for a given event
vector<vector<parton> > eventParticles;

//Vector of scattering times
vector<parton> scatteringTimes;

//Counter for the number of events processed
int evtnumber = 0;

//Time step for animation, in fm/c
float dt = 0.05;

//Number of iterations in animation
int NITER = 100;

//Amount of time during which we display a collision marker
float collMarkTime = 4*dt;

//-----------------------------------
// Functions
//-----------------------------------

void draw(vector<float> x, vector<float> y, int iteration, vector<float> v2)
{
  TCanvas *c = new TCanvas(Form("c_%i",iteration),Form("c_%i",iteration),600,800);
  gStyle->SetOptStat(0);

  //Divide canvas into a pad for the event display, and another for v2(t)
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0, 0.2, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0, 0.0, 1.0, 0.2);

  pad1->Draw();
  pad2->Draw();

  //Go to pad with v2(t)
  pad2->cd();
  TH1F *hTemplate_v2 = new TH1F(Form("hTemplate_v2_%i",iteration),"",NITER,0,NITER*dt);
  hTemplate_v2->GetXaxis()->SetTitleFont(62);
  hTemplate_v2->GetYaxis()->SetTitleFont(62);
  hTemplate_v2->GetXaxis()->SetLabelFont(62);
  hTemplate_v2->GetYaxis()->SetLabelFont(62);
  hTemplate_v2->GetYaxis()->SetRangeUser(0,0.25);
  hTemplate_v2->SetFillColor(kBlue);

  for(int i=1; i<=v2.size(); i++)
    {
      hTemplate_v2->SetBinContent(i,v2[i-1]);
    }

  hTemplate_v2->Draw();

  //Go to pad with scattering animation
  pad1->cd();
  TH2F *hTemplate = new TH2F(Form("hTemplate_%i",iteration),Form("hTemplate_%i",iteration),100,-10,10,100,-10,10);
  hTemplate->SetTitle("");
  hTemplate->GetXaxis()->SetTitle("x [fm]");
  hTemplate->GetYaxis()->SetTitle("y [fm]");
  hTemplate->GetXaxis()->SetTitleFont(62);
  hTemplate->GetYaxis()->SetTitleFont(62);
  hTemplate->GetXaxis()->SetLabelFont(62);
  hTemplate->GetYaxis()->SetLabelFont(62);
  hTemplate->Draw();

  //Iterate over partons position vectors and draw TEllipse
  for(int i=0; i<x.size(); i++)
    {
      TEllipse *tell = new TEllipse(x[i], y[i], 0.05, 0.05);
      tell->SetFillColor(kBlack);
      tell->Draw("same");
    }

  TLatex *tTime = new TLatex(0.62,0.8,Form("t = %g",iteration*dt));
  tTime->SetNDC(kTRUE);
  tTime->Draw("same");
  TLatex *tTimeUnits = new TLatex(0.78,0.8,"fm/c");
  tTimeUnits->SetNDC(kTRUE);
  tTimeUnits->Draw("same");

  //Check if we are within 'collMarkTime' following a scattering event
  for(int i=0; i<scatteringTimes.size(); i++)
    {
      if(iteration*dt >= scatteringTimes[i].t && iteration*dt <= scatteringTimes[i].t + collMarkTime)
	{
	  TEllipse *tCollMark = new TEllipse(scatteringTimes[i].x, scatteringTimes[i].y, 0.2, 0.2);
	  tCollMark->SetLineColor(kRed);
	  tCollMark->SetFillStyle(0);
	  tCollMark->Draw("same");
	}
    }

  
  if(iteration < 10)
    {
      c->SaveAs(Form("IterationFrame_00%i.gif",iteration));
    }
  else if(iteration >= 10 && iteration < 99)
    {
      c->SaveAs(Form("IterationFrame_0%i.gif",iteration));
    }
  else
    {
      c->SaveAs(Form("IterationFrame_%i.gif",iteration));
    }
  
}

void computePosition(float &x, float &y, vector<parton> evol, float t0)
{
  //Determine if formation has occurred for the parton
  float t_form = evol[0].t;
  if(t0 < t_form)
    {
      x = -999.0;
      y = -999.0;
      return;
    }  

  //If parton has already formed, find which scattering stage it is in
  int stage = 0;
  int numStages = evol.size()-1;

  if(t0 >= evol[numStages].t)
    { 
      stage = numStages; 
    }
  else
    { 
      for(int i=0; i<numStages; i++)
	{
	  if(t0 >= evol[i].t && t0 < evol[i+1].t) stage = i;
	}
    }

  //Compute energy and velocity using the momentum of the particular scattering stage
  //Velocity in units of c
  float px = evol[stage].px;
  float py = evol[stage].py;
  float pz = evol[stage].pz;
  float mass = evol[stage].m;
  float energy = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
  float vel = 1 - 0.5*(mass/energy)*(mass/energy);
  float p = TMath::Sqrt(px*px + py*py + pz*pz);

  //Compute orientation of particle
  float theta = TMath::ACos(pz/p);
  float phi = TMath::ATan2(py,px);

  //Compute new coordinate after time t0
  float distTraveled = vel*t0;
  float xTraveled = distTraveled*TMath::Sin(theta)*TMath::Cos(phi);
  float yTraveled = distTraveled*TMath::Sin(theta)*TMath::Sin(phi);
  float zTraveled = distTraveled*TMath::Cos(theta);

  //Offset by previous location
  float x0 = evol[stage].x;
  float y0 = evol[stage].y;
  x = x0 + xTraveled;
  y = y0 + yTraveled;
}

void loadEventPlane()
{
  ifstream psi_2_file;
  psi_2_file.open("psi_2.txt");

  while(psi_2_file)
    {
      psi_2_file >> psi2;
    }

  psi_2_file.close();
}

float computeEllipticFlow(vector<float> x, vector<float> y)
{
  float v2 = 0.0;
  int nPartons = x.size();

  for(int i=0; i<nPartons; i++)
    {
      float x_pos = x[i];
      float y_pos = y[i];
      float phi = TMath::ATan2(y_pos,x_pos);

      v2 += TMath::Cos(2*(phi-psi2));
    }

  v2 = (float) v2/nPartons;
  return v2;
}

void processEvent()
{
  loadEventPlane();

  int numPartons = eventParticles.size();

  //Determine the x,y,t coordinates of each scattering event
  for(int i=0; i<numPartons; i++)
    {
      vector<parton> v = eventParticles[i];

      for(int j=1; j<v.size(); j++)
	{
	  parton scatt;
	  scatt.x = v[j].x;
	  scatt.y = v[j].y;
	  scatt.t = v[j].t;

	  scatteringTimes.push_back(scatt);
	}
    }

  //Loop over all particles in N iterations to find their positions
  vector<float> xvals;
  vector<float> yvals;
  vector<float> v2vals;

  int iteration = 0;
  while(iteration < NITER)
    {
      for(int i=0; i<numPartons; i++)
	{
	  vector<parton> v = eventParticles[i];
	  
	  float x = 0;
	  float y = 0;
	  computePosition(x, y, v, dt*iteration);
	  xvals.push_back(x);
	  yvals.push_back(y);
	}
 
      float v2 = computeEllipticFlow(xvals, yvals);
      v2vals.push_back(v2);

      draw(xvals, yvals, iteration, v2vals);
      xvals.clear();
      yvals.clear();
      iteration++;
    }
}

void animateEvent(char *initialInfoFile = "parton-initial-afterPropagation.dat", char *evolInfoFile = "parton-collisionsHistory.dat", char *outputFile = "evolution_out.root")
{
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
      //initializeDisplay();
      processEvent();
      eventParticles.clear();
    }
}
