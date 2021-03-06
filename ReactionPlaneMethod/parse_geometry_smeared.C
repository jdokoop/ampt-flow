// Macro to parse output geometry files from AMPT
// Apr. 27 2015
// J. Orjuela-Koop

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

#include "TMath.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"

using namespace std;

//---------------------------------------------------
// Variables
//---------------------------------------------------

//Simplified particle object, with just (x,y) position
struct particle
{
	double x;
	double y;
};

//Participants in the target and projectile nuclei
vector<particle> particles_targ;
vector<particle> particles_proj;

double epsilon_2 = 0.0; // jln - really these should be initialized to zero
double epsilon_3 = 0.0;
double psi_2 = 0.0;
double psi_3 = 0.0;
double ncollsum  = 0.0; // jln - add calculation of mean Ncoll
double npartsumT = 0.0; //     - add calculation of mean Npart_targ
double npartsumP = 0.0; //     - add calculation of mean Npart_proj

//Files to output event plane angles
ofstream psi2_file;
ofstream psi3_file;
ofstream ncoll_file;

//---------------------------------------------------
// Functions
//---------------------------------------------------

//Determines if a given target particle has already been taken into account
bool exists_targ(particle p)
{
	for(unsigned int i=0; i<particles_targ.size(); i++)
	{
		particle aux = particles_targ[i];

		if(aux.x == p.x && aux.y == p.y)
		{
			return true;
		}
	}

	return false;
}

//Determines if a given projectile particle has already been taken into account
bool exists_proj(particle p)
{
	for(unsigned int i=0; i<particles_proj.size(); i++)
	{
		particle aux = particles_proj[i];

		if(aux.x == p.x && aux.y == p.y)
		{
			return true;
		}
	}

	return false;
}

//Get phi angle from cartesian coordinates
double get_phi(double x, double y)
{
	return TMath::ATan2(y,x);
}

//Get distance from center in polar coordinates
double get_r(double x, double y)
{
	return TMath::Sqrt(x*x + y*y);
}

void process_event()
{  
	bool WithProj = true;

	double x_cm = 0;
	double y_cm = 0;

	//Compute centroid
	for(unsigned int j=0; j<particles_targ.size(); j++)
	{
		x_cm += particles_targ[j].x;
		y_cm += particles_targ[j].y;
	}

	double Count = (double) particles_targ.size();

	if (WithProj) {
		for(unsigned int j=0; j<particles_proj.size(); j++)
		{
			x_cm += particles_proj[j].x;
			y_cm += particles_proj[j].y;
		}
		Count += (double) particles_proj.size();
	}

	x_cm = x_cm/Count;
	y_cm = y_cm/Count;

	//Shift origin to centroid
	for(unsigned int j=0; j<particles_targ.size(); j++)
	{
		particles_targ[j].x = particles_targ[j].x - x_cm;
		particles_targ[j].y = particles_targ[j].y - y_cm;
	}

	if (WithProj) {
		//Shift origin to centroid
		for(unsigned int j=0; j<particles_proj.size(); j++)
		{
			particles_proj[j].x = particles_proj[j].x - x_cm;
			particles_proj[j].y = particles_proj[j].y - y_cm;
		}
	}

	//Compute eccentricity
	//JLN - Option to do with blurred Gaussian distributions... (does not change above x_cm,y_cm calculation)
	bool NucleonSmear = true;
	double SmearSigma = 0.4; // units of fm
	TF1 *fradius = new TF1("fradius","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradius->SetParameter(0,SmearSigma);
	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());
	int SmearSamplings = 100;
	if (NucleonSmear) Count = 0;

	double qx_2 = 0;
	double qy_2 = 0;
	double qx_3 = 0;
	double qy_3 = 0;
	double phi = 0;
	double rsq = 0;
	double r = 0;

	double e2 = 0;
	double e3 = 0;

	//Loop over all particles in the target nucleus and compute eccentricity
	for(unsigned int i=0; i<particles_targ.size(); i++)
	{
		if (!NucleonSmear) {
			r = get_r(particles_targ[i].x,particles_targ[i].y);
			phi = get_phi(particles_targ[i].x,particles_targ[i].y);

			qx_2 += r*r*TMath::Cos(2*phi);
			qy_2 += r*r*TMath::Sin(2*phi);
			qx_3 += r*r*TMath::Cos(3*phi);
			qy_3 += r*r*TMath::Sin(3*phi);
			rsq += r*r;
		}
		if (NucleonSmear) {
			for (int is=0;is<SmearSamplings;is++) {
				double rtemp   = fradius->GetRandom();
				double phitemp = fphi->GetRandom();
				double xtemp   = particles_targ[i].x + rtemp*TMath::Sin(phitemp);
				double ytemp   = particles_targ[i].y + rtemp*TMath::Cos(phitemp);

				r    = TMath::Sqrt(xtemp*xtemp + ytemp*ytemp);
				phi  = TMath::ATan2(ytemp,xtemp);

				qx_2 += r*r*TMath::Cos(2*phi);
				qy_2 += r*r*TMath::Sin(2*phi);
				qx_3 += r*r*TMath::Cos(3*phi);
				qy_3 += r*r*TMath::Sin(3*phi);
				rsq += r*r;

				Count = Count + 1;  // noting this is a double
			}
		}
	}

	if (WithProj) {
		for(unsigned int i=0; i<particles_proj.size(); i++)
		{
			if (!NucleonSmear) {
				r = get_r(particles_proj[i].x,particles_proj[i].y);
				phi = get_phi(particles_proj[i].x,particles_proj[i].y);

				qx_2 += r*r*TMath::Cos(2*phi);
				qy_2 += r*r*TMath::Sin(2*phi);
				qx_3 += r*r*TMath::Cos(3*phi);
				qy_3 += r*r*TMath::Sin(3*phi);
				rsq += r*r;
			}
			if (NucleonSmear) {
				for (int is=0;is<SmearSamplings;is++) {
					double rtemp   = fradius->GetRandom();
					double phitemp = fphi->GetRandom();
					double xtemp   = particles_proj[i].x + rtemp*TMath::Sin(phitemp);
					double ytemp   = particles_proj[i].y + rtemp*TMath::Cos(phitemp);

					r    = TMath::Sqrt(xtemp*xtemp + ytemp*ytemp);
					phi  = TMath::ATan2(ytemp,xtemp);

					qx_2 += r*r*TMath::Cos(2*phi);
					qy_2 += r*r*TMath::Sin(2*phi);
					qx_3 += r*r*TMath::Cos(3*phi);
					qy_3 += r*r*TMath::Sin(3*phi);
					rsq += r*r;

					Count = Count + 1;  // noting this is a double
				}
			}

		}
	}

	qx_2 = qx_2/Count;
	qy_2 = qy_2/Count;
	qx_3 = qx_3/Count;
	qy_3 = qy_3/Count;
	rsq  = rsq /Count;

	e2 = TMath::Sqrt(qx_2*qx_2 + qy_2*qy_2)/rsq;
	e3 = TMath::Sqrt(qx_3*qx_3 + qy_3*qy_3)/rsq;

	psi_2 = (TMath::ATan2(qy_2,qx_2) + TMath::Pi())/2;
	psi_3 = (TMath::ATan2(qy_3,qx_3) + TMath::Pi())/3;

	psi2_file << psi_2 << "\n";
	psi3_file << psi_3 << "\n";

	//cout << "e2 = " << e2 << endl;
	//cout << "psi_2 = " << psi_2 << endl;
	if (Count == 1) e2 = 0.0;

	//Accumulate event-by-event eccentricity to compute average for all events
	epsilon_2 += e2;
	epsilon_3 += e3;
}

void parse_geometry_smeared()
{
	//Read input file
	//File format for each event:
	//   Ncoll
	//   Projectile participant x y z
	//   Target participant x y z
	//   ...
	string fname = "geometry.dat";
	ifstream infile(fname.c_str());

	string line;
	int ncoll = 0;
	int nevent = 0;

	double a = 0;
	double b = 0;
	double c = 0;

	psi2_file.open ("psi_2.txt");
	psi3_file.open ("psi_3.txt");
	ncoll_file.open ("ncoll_file.txt");

	while(getline(infile, line))
	{
		istringstream iss(line);

		//Do we have a new event? You can tell if there's only one number in the line, corresponding to ncoll
		if(!(iss >> a >> b >> c))
		{
			ncoll = a;
			ncoll_file << ncoll << "\n";
		}
		else
		{
			continue;
		}

		//There is a pair of lines (projectile, target) for each binary collision in the event.
		//Hence, iterate over 2*ncoll lines.
		for(int i=1; i<=2*ncoll; i++)
		{
			particle part;

			getline(infile, line);
			istringstream iss(line);
			iss >> a >> b >> c;

			part.x = a;
			part.y = b;

			if(i%2 == 0) //The second line in the pair corresponds to target participants
			{
				if(!exists_targ(part))
				{
					particles_targ.push_back(part);
				}
			}
			else //The first line in the pair corresponds to projectile participants
			{
				if(!exists_proj(part))
				{
					particles_proj.push_back(part);
				}
			}
			a = 0;
			b = 0;
		}

		//Process the particles in the event to compute epsilon2,3
		process_event();
		ncollsum += ncoll;
		npartsumT += particles_targ.size();
		npartsumP += particles_proj.size();

		//Reset after event
		nevent++;
		ncoll = 0;
		particles_targ.clear();
		particles_proj.clear();
	}

	//Average epsilon over all events
	epsilon_2 = epsilon_2/(double)nevent;  // jln - should not divide by an int, cast it first
	epsilon_3 = epsilon_3/(double)nevent;

	cout << "***RAN OVER " << nevent << " EVENTS" << endl;
	cout << "<epsilon_2> = " << epsilon_2 << endl;
	cout << "<epsilon_3> = " << epsilon_3 << endl;

	cout << "<ncoll> = " << ncollsum/(double)nevent << endl;
	cout << "<npart_targ> = " << npartsumT/(double)nevent << endl;
	cout << "<npart_proj> = " << npartsumP/(double)nevent << endl;

	psi2_file.close();
	psi3_file.close();
	ncoll_file.close();

}
