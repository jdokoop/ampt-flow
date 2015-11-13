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

//Collision system to display
int collSystem = 2; //0 = p+Au, 1 = d+Au, 2 = He3+Au

//-----------------------------------
// Functions
//-----------------------------------

void draw(vector<float> x_green, vector<float> y_green, vector<float> x_blue, vector<float> y_blue, vector<float> x_yellow, vector<float> y_yellow, int iteration, vector<float> v2, vector<float> epsilon2)
{
  TCanvas *c = new TCanvas(Form("c_%i",iteration),Form("c_%i",iteration),480,800);
  gStyle->SetOptStat(0);

  //Divide canvas into a pad for the event display, v2(t), and epsilon2(t)
  TPad *pad1 = new TPad("pad1", "The pad 60% of the height",0.0, 0.4, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0, 0.2, 1.0, 0.4);
  TPad *pad3 = new TPad("pad3", "The pad 20% of the height",0.0, 0.0, 1.0, 0.2);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();

  //Go to pad with epsilon2(t)
  pad3->cd();
  pad3->SetFillColor(kBlack);
  pad3->SetBottomMargin(0.3);
  pad3->SetTickx();
  pad3->SetTicky();

  TH1F *hTemplate_ep2 = new TH1F(Form("hTemplate_ep2_%i",iteration),";t [fm/c]; #varepsilon_{2}",NITER,0,NITER*dt);
  hTemplate_ep2->GetXaxis()->SetTitleFont(62);
  hTemplate_ep2->GetYaxis()->SetTitleFont(62);
  hTemplate_ep2->GetXaxis()->SetLabelFont(62);
  hTemplate_ep2->GetYaxis()->SetLabelFont(62);
  hTemplate_ep2->GetYaxis()->SetTitleOffset(0.3);
  hTemplate_ep2->GetYaxis()->SetTitleSize(0.15);
  hTemplate_ep2->GetYaxis()->SetLabelSize(0.09);
  hTemplate_ep2->GetYaxis()->SetRangeUser(0,1);
  hTemplate_ep2->GetYaxis()->SetNdivisions(8);

  hTemplate_ep2->GetXaxis()->SetAxisColor(kWhite);
  hTemplate_ep2->GetYaxis()->SetAxisColor(kWhite);
  hTemplate_ep2->GetXaxis()->SetTitleColor(kWhite);
  hTemplate_ep2->GetYaxis()->SetTitleColor(kWhite);
  hTemplate_ep2->GetXaxis()->SetLabelColor(kWhite);
  hTemplate_ep2->GetYaxis()->SetLabelColor(kWhite);

  hTemplate_ep2->GetXaxis()->SetTitleOffset(0.8);
  hTemplate_ep2->GetXaxis()->SetTitleSize(0.12);
  hTemplate_ep2->GetXaxis()->SetLabelSize(0.09);

  hTemplate_ep2->SetLineColor(kSpring-9);
  hTemplate_ep2->SetLineWidth(2);

  for(int i=1; i<=epsilon2.size(); i++)
    {
      hTemplate_ep2->SetBinContent(i,epsilon2[i-1]);
    }

  hTemplate_ep2->Draw("L");

  //Go to pad with v2(t)
  pad2->cd();
  pad2->SetFillColor(kBlack);
  pad2->SetBottomMargin(0.3);
  pad2->SetTickx();
  pad2->SetTicky();

  TH1F *hTemplate_v2 = new TH1F(Form("hTemplate_v2_%i",iteration),";t [fm/c]; v_{2}",NITER,0,NITER*dt);
  hTemplate_v2->GetXaxis()->SetTitleFont(62);
  hTemplate_v2->GetYaxis()->SetTitleFont(62);
  hTemplate_v2->GetXaxis()->SetLabelFont(62);
  hTemplate_v2->GetYaxis()->SetLabelFont(62);
  hTemplate_v2->GetYaxis()->SetTitleOffset(0.3);
  hTemplate_v2->GetYaxis()->SetTitleSize(0.15);
  hTemplate_v2->GetYaxis()->SetLabelSize(0.09);
  hTemplate_v2->GetYaxis()->SetRangeUser(-0.1,0.1);
  hTemplate_v2->GetYaxis()->SetNdivisions(8);

  hTemplate_v2->GetXaxis()->SetAxisColor(kWhite);
  hTemplate_v2->GetYaxis()->SetAxisColor(kWhite);
  hTemplate_v2->GetXaxis()->SetTitleColor(kWhite);
  hTemplate_v2->GetYaxis()->SetTitleColor(kWhite);
  hTemplate_v2->GetXaxis()->SetLabelColor(kWhite);
  hTemplate_v2->GetYaxis()->SetLabelColor(kWhite);

  hTemplate_v2->GetXaxis()->SetTitleOffset(0.8);
  hTemplate_v2->GetXaxis()->SetTitleSize(0.12);
  hTemplate_v2->GetXaxis()->SetLabelSize(0.09);

  hTemplate_v2->SetLineColor(kOrange+0);
  hTemplate_v2->SetLineWidth(2);

  for(int i=1; i<=v2.size(); i++)
    {
      hTemplate_v2->SetBinContent(i,v2[i-1]);
    }

  hTemplate_v2->Draw("L");

  //Go to pad with scattering animation
  pad1->cd();
  pad1->SetFillColor(kBlack);
  pad1->SetTickx();
  pad1->SetTicky();

  TH2F *hTemplate = new TH2F(Form("hTemplate_%i",iteration),Form("hTemplate_%i",iteration),100,-5,5,100,-5,5);
  hTemplate->SetTitle("");
  hTemplate->GetXaxis()->SetTitle("x [fm]");
  hTemplate->GetYaxis()->SetTitle("y [fm]");
  hTemplate->GetXaxis()->SetTitleFont(62);
  hTemplate->GetYaxis()->SetTitleFont(62);
  hTemplate->GetXaxis()->SetLabelFont(62);
  hTemplate->GetYaxis()->SetLabelFont(62);

  hTemplate->GetXaxis()->SetAxisColor(kWhite);
  hTemplate->GetYaxis()->SetAxisColor(kWhite);
  hTemplate->GetXaxis()->SetTitleColor(kWhite);
  hTemplate->GetYaxis()->SetTitleColor(kWhite);
  hTemplate->GetXaxis()->SetLabelColor(kWhite);
  hTemplate->GetYaxis()->SetLabelColor(kWhite);
  
  hTemplate->Draw();

  //Draw participant plane
  TLine *lPlane = new TLine(0, 0, 10*TMath::Cos(psi2), 10*TMath::Sin(psi2));
  lPlane->SetLineColor(kOrange-3);
  lPlane->SetLineWidth(2);
  lPlane->SetLineStyle(2);
  //lPlane->Draw("same");

  //Iterate over partons position vectors and draw TEllipse
  for(int i=0; i<x_green.size(); i++)
    {
      TEllipse *tell = new TEllipse(x_green[i], y_green[i], 0.09, 0.09);
      tell->SetFillColor(kPink+7);
      tell->Draw("same");
    }

  for(int i=0; i<x_blue.size(); i++)
    {
      TEllipse *tell = new TEllipse(x_blue[i], y_blue[i], 0.09, 0.09);
      tell->SetFillColor(kSpring-8);
      tell->Draw("same");
    }

  for(int i=0; i<x_yellow.size(); i++)
    {
      TEllipse *tell = new TEllipse(x_yellow[i], y_yellow[i], 0.09, 0.09);
      tell->SetFillColor(kAzure+1);
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
	  TEllipse *tCollMark = new TEllipse(scatteringTimes[i].x, scatteringTimes[i].y, 0.25, 0.25);
	  tCollMark->SetLineColor(kWhite);
	  tCollMark->SetLineWidth(3);
	  //tCollMark->SetFillColor(kRed);
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

void computePosition(float &x, float &y, float &phi_angle, vector<parton> evol, float t0)
{
  //Determine if formation has occurred for the parton
  float t_form = evol[0].t;
  if(t0 < t_form)
    {
      x = -999.0;
      y = -999.0;
      phi_angle = -999.0;
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
  TLorentzVector ev(px,py,pz,energy);
  double phi = ev.Phi();
  double theta = ev.Theta();
  phi_angle = phi;

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

  if(collSystem == 0)
    {
      psi_2_file.open("psi_2.txt");
    }
  if(collSystem == 1)
    {
      psi_2_file.open("psi_2_He3Au.txt");
    }
  if(collSystem == 2)
    {
      psi_2_file.open("psi_2_He3Au.txt");
    }

  while(psi_2_file)
    {
      psi_2_file >> psi2;
    }

  psi_2_file.close();
}

float computeEllipticFlow(vector<float> phi)
{
  float v2 = 0.0;
  int nPartons = 0;

  for(int i=0; i<phi.size(); i++)
    {
      //Ignore partons that haven't yet formed
      if(phi[i] < -900)
	{
	  continue;
	}

      float angle = phi[i];
      v2 += TMath::Cos(2*(angle-psi2));
      nPartons++;
    }

  v2 = (float) v2/nPartons;
  if(v2 != v2) return 0; //If we get a NaN (i.e. partons not yet formed), return -999
  return v2;
}

float computeEpsilon2(vector<float> x_green, vector<float> y_green, vector<float> x_yellow, vector<float> y_yellow, vector<float> x_blue, vector<float> y_blue)
{
  //Concatenate all partons in a single vector
  vector<float> x;
  vector<float> y;

  for(int i=0; i<x_green.size(); i++)
    {
      x.push_back(x_green[i]);
      y.push_back(y_green[i]);
    }

  for(int i=0; i<x_blue.size(); i++)
    {
      x.push_back(x_blue[i]);
      y.push_back(y_blue[i]);
    }

  for(int i=0; i<x_yellow.size(); i++)
    {
      x.push_back(x_yellow[i]);
      y.push_back(y_yellow[i]);
    }

  //Start by computing center of mass
  int nPartons = 0;

  float x_cm = 0.0;
  float y_cm = 0.0;

  for(int i=0; i<x.size(); i++)
    {
      //Only use partons that have formed, otherwise their coordinates are (-999,-999)
      if(x[i] < -900 && y[i] < -900)
	{
	  continue;
	}

      x_cm = x_cm + x[i];
      y_cm = y_cm + y[i];
      nPartons++;
    }

  x_cm = x_cm/nPartons;
  y_cm = y_cm/nPartons;

  //Compute eccentricity with respect to center of mass
  float qx = 0.0;
  float qy = 0.0;
  float r2 = 0.0;

  for(int i=0; i<x.size(); i++)
    {
      if(x[i] < -900 && y[i] < -900)
	{
	  continue;
	}

      float x_pos = x[i] - x_cm;
      float y_pos = y[i] - y_cm;
      float r = TMath::Sqrt(x_pos*x_pos + y_pos*y_pos);
      float phi = TMath::ATan2(y_pos,x_pos);

      qx += r*r*TMath::Cos(2*phi);
      qy += r*r*TMath::Sin(2*phi);
      r2 += r*r;
    }

  qx = qx/nPartons;
  qy = qy/nPartons;
  r2 = r2/nPartons;

  float ep2 = TMath::Sqrt(qx*qx + qy*qy)/r2;
  if(ep2 != ep2) return 0; //If NaN (i.e. partons not yet formed), return -999
  return ep2;
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

  //If d+Au, for nucleons with y_init > 0
  //If p+Au, for nucleons that don't scatter
  vector<float> xvals_green;
  vector<float> yvals_green;

  //If d+Au, for nucleons with y_init < 0
  //If p+Au, for nucleons with Nscatt = 1
  vector<float> xvals_yellow;
  vector<float> yvals_yellow;

  //If d+Au, don't use these vectors
  //If p+Au, for nucleons with Nscatt >= 2
  vector<float> xvals_blue;
  vector<float> yvals_blue;

  vector<float> phivals;
  vector<float> v2vals;
  vector<float> epsilon2vals;

  int iteration = 0;
  while(iteration < NITER)
    {
      for(int i=0; i<numPartons; i++)
	{
	  vector<parton> v = eventParticles[i];
	  
	  float x = 0;
	  float y = 0;
	  float phi = 0;

	  computePosition(x, y, phi, v, dt*iteration);

	  //Sort coordinates into different vectors depending on Nscatt or initial hotspot location (for pAu, dAu respectively)
	  if(collSystem == 1)
	    {
	      if(v[0].y < 0)
		{
		  xvals_yellow.push_back(x);
		  yvals_yellow.push_back(y);
		}
	      else
		{
		  xvals_green.push_back(x);
		  yvals_green.push_back(y);
		}
	    }
	  else if(collSystem == 0)
	    {
	      if(v.size() == 1)
		{
		  xvals_green.push_back(x);
		  yvals_green.push_back(y);
		}
	      else if(v.size() == 2)
		{
		  xvals_yellow.push_back(x);
		  yvals_yellow.push_back(y);
		}
	      else if(v.size() >= 3)
		{
		  xvals_blue.push_back(x);
		  yvals_blue.push_back(y);
		}
	    }
	  else if(collSystem == 2)
	    {
	      if(v[0].y < -1)
		{
		  xvals_yellow.push_back(x);
		  yvals_yellow.push_back(y);
		}
	      else if(v[0].y < 0 && v[0].y > -1)
		{
		  xvals_green.push_back(x);
		  yvals_green.push_back(y);
		}
	      else if(v[0].y > 0)
		{
		  xvals_blue.push_back(x);
		  yvals_blue.push_back(y);
		}
	    }

	  phivals.push_back(phi);
	}
 
      float v2 = computeEllipticFlow(phivals);
      v2vals.push_back(v2);

      float epsilon2 = computeEpsilon2(xvals_green, yvals_green, xvals_yellow, yvals_yellow, xvals_blue, yvals_blue);
      epsilon2vals.push_back(epsilon2);

      draw(xvals_green, yvals_green, xvals_blue, yvals_blue, xvals_yellow, yvals_yellow, iteration, v2vals, epsilon2vals);
      xvals_green.clear();
      yvals_green.clear();
      xvals_yellow.clear();
      yvals_yellow.clear();
      xvals_blue.clear();
      yvals_blue.clear();
      phivals.clear();
      iteration++;
    }
}

void animateEvent(char *initialInfoFile = "", char *evolInfoFile = "")
{
  if(collSystem == 0)
    {
      initialInfoFile = "parton-initial-afterPropagation.dat";
      evolInfoFile = "parton-collisionsHistory.dat";
    }
  else if(collSystem == 1)
    {
      initialInfoFile = "parton-initial-afterPropagation_He3Au.dat";
      evolInfoFile = "parton-collisionsHistory_He3Au.dat";
    }
  else if(collSystem == 2)
    {
      initialInfoFile = "parton-initial-afterPropagation_He3Au.dat";
      evolInfoFile = "parton-collisionsHistory_He3Au.dat";
    }

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
      processEvent();
      eventParticles.clear();
    }
}
