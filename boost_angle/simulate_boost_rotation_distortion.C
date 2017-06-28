#include <TLorentzVector.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1F.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

// a macro written in C++ with objects in the ROOT framework
// meant to be run with ROOT.

// this function is meant to discern the effect of boost + rotation and just rotation

class particle{// particle position class

public:

	particle(float a, float b, float c){x = a; y = b; z = c;}

	float getX(){return x;}
	float getY(){return y;}
	float getZ(){return z;}

	void setX(float a){x = a;}
	void setY(float a){y = a;}
	void setZ(float a){z = a;}

	float getPhi(){return atan2(y,x);}

protected:
	float x;
	float y;
	float z;

};

void simulate_boost_rotation_distortion()
{

	vector<particle> pions;

	int npions = 10000000;

	double pi = TMath::Pi();

	gRandom->SetSeed(0);

	for(int i = 0; i < npions; i++)
	{

		float phi = gRandom->Uniform(-pi, pi);

		float rand_px = 0.25 * TMath::Cos(phi);
		float rand_py = 0.25 * TMath::Sin(phi);
		float pz = 0.25 * TMath::SinH(-3.4);
		particle tmp(rand_px,rand_py,pz);
		pions.push_back(tmp);
	}

	TH1D *azimuth = new TH1D("azimuth","azimuth",100,-pi,pi);
	TH1D *azimuth2 = new TH1D("azimuth2","azimuth2",100,-pi,pi);
	TH1D *azimuth3 = new TH1D("azimuth3","azimuth3",100,-pi,pi);

	for(int i = 0; i < npions; i++)
	{
		float phi = pions[i].getPhi();
		azimuth->Fill(phi);
	}

	TLegend *leg = new TLegend(0.15,0.75,0.35,0.9);

	TCanvas *c0 = new TCanvas("c0","c0",800,800);
	azimuth->GetYaxis()->SetRangeUser(0.0,azimuth->GetMaximum()*1.5);
	azimuth->SetXTitle("#phi");
	azimuth->SetLineColor(kBlack);
	azimuth->SetLineWidth(2.0);
	azimuth->Draw("");
	leg->AddEntry(azimuth,"no modification");
	for(int i = 0; i < npions; i++)
	{
		float angle = 0.003;

		float x = pions[i].getX();
		float z = pions[i].getZ();

		float new_x = x * TMath::Cos(angle) - z * TMath::Sin(angle);
		float new_z = x * TMath::Sin(angle) + z * TMath::Cos(angle);
	
		pions[i].setX(new_x);
		pions[i].setZ(new_z);

		azimuth2->Fill(pions[i].getPhi());
	}

	azimuth2->SetLineColor(kBlue);
	azimuth2->SetLineWidth(2.0);
	azimuth2->Draw("same");
	leg->AddEntry(azimuth2,"with 3 mRad rotation");

	float px1 = 100*TMath::Sin(0.0016);
	float py1 = 0.0;
	float pz1 = 100*TMath::Cos(0.0016);
	float mass1 = 0.9;

	float energy1 =  sqrt (pow(px1, 2) + pow(py1, 2) + pow(pz1, 2) + pow(mass1, 2));
	TLorentzVector ev1(px1, py1, pz1, energy1);

	float px2 = 100*TMath::Sin(TMath::Pi()+0.0036);
	float py2 = 0.0;
	float pz2 = 100*TMath::Cos(TMath::Pi()+0.0036);
	float mass2 = 197*0.9;

	float energy2 =  sqrt (pow(px2, 2) + pow(py2, 2) + pow(pz2, 2) + pow(mass2, 2));
	TLorentzVector ev2(px2, py2, pz2, energy2);

	TLorentzVector cms = ev1 + ev2;

	for(int i = 0; i < npions; i++)
	{

		float px = pions[i].getX();
		float py = pions[i].getY();
		float pz = pions[i].getZ();

		TLorentzVector vec(px,py,pz,TMath::Sqrt(pow(px,2)+pow(py,2)+pow(pz,2)+pow(0.139,2)));

		vec.Boost(-cms.BoostVector());

		pions[i].setX(vec.Px());
		pions[i].setY(vec.Py());
		pions[i].setZ(vec.Pz());

		azimuth3->Fill(pions[i].getPhi());
	}

	azimuth3->SetLineColor(kRed);
	azimuth3->SetLineWidth(2.0);
	azimuth3->Draw("same");
	leg->AddEntry(azimuth3,"with CMS boost AND rotation");
	leg->Draw("same");
}