#include <TLorentzVector.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1F.h>

#include <iostream>

using namespace std;

// a macro written in C++ with objects in the ROOT framework
// meant to be run with ROOT.

// calculate modified eta coverage with offset and rotation as a function of phi angle

// 1) using FVTX 0-3 layer rmin, rmax, zplane and a zvertex, calculate two vectors
// 2) determine the pseudorapidity coverage of the two vectors (integrate AMPT particles in that range)
// 3) offset, boost and rotate vectors and then recalculate...
// 4) plot etamin, etamax modification versus phi and the weight factor

// start with just one plane and zvertex = 0.0
// UPDATE - loop over all four FVTX planes and z-vertex bins and write array to output file

void eta_modification_grid_original()
{

  // initialization for TStyle object
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.3f");

  // hard coded values for beam offset
  double zvertex = 0.0;
  double xoffset = 0.191;
  double yoffset = 0.057; // these are both positive values
  double beamangle = -0.0036;

  // input distribution in eta from AMPT simulations
  TFile *fin = new TFile("fout_eta_ampt_paub2.root");
  TH1F  *hampt = static_cast <TH1F *> (fin->Get("heta"));

  // layers 0 - 3
  double fvtx_z[5] = {-20.11,-26.14,-32.17,-38.20, -144.4 };  // last one is BBC
  double fvtx_r_min[5] = {4.5,4.5,4.5,4.5,          5.8};
  double fvtx_r_max[5] = {10.5,17.2,17.2,17.2,      13.0};     // update BBC values...

  TGraph *gridgraph[5][60];
  TGraph *etamin  = new TGraph();
  TGraph *etamax  = new TGraph();
  TGraph *etaminM = new TGraph();
  TGraph *etamaxM = new TGraph();

  // index FVTX layers 0,1,2,3 and BBC = 4
  for (int idet=0;idet<1;idet++) {
    // 60 bins in zvertex from -30 to +30
    for (int iz=30;iz<31;iz++) {

      zvertex = -29.5 + ((float) iz) * 1.0;
      double rmin   = fvtx_r_min[idet];
      double rmax   = fvtx_r_max[idet];
      double zplane = fvtx_z[idet];

      char fooname[60];
      sprintf(fooname,"multcorr_layer_%1d_zvertex_%d",idet,iz);
      gridgraph[idet][iz] = new TGraph();
      gridgraph[idet][iz]->SetName(fooname);

      TGraph *particles = new TGraph();
      TGraph *particlesM = new TGraph();
      TGraph *modification = new TGraph();
      
      // step uniformly around in phi angle
      for(int ipart = 0; ipart < 100; ipart++) {

        float phi = -TMath::Pi() + 0.01 * ((float) ipart) * 2.0 * TMath::Pi();

  for (int ii=0;ii<2;ii++) { // min and max pieces

    float r;
    if (ii==0) r = rmin;
    if (ii==1) r = rmax;
    
    // really x, y, z coordinate for the detector element
    float x = r * TMath::Cos(phi);
    float y = r * TMath::Sin(phi);
    float z = zplane - zvertex;
    float eta = TMath::ATanH(z/TMath::Sqrt(x*x+y*y+z*z));
    
    // then x,y,z prime with x and y offsets
    float xprime = x - xoffset;
    float yprime = y - yoffset;
    float zprime = z;  // no additional offset
    
    // then rotate which interchanges a bit of x,z
    float xdoubleprime = TMath::Cos(beamangle) * xprime - TMath::Sin(beamangle) * zprime;
    float ydoubleprime = yprime;
    float zdoubleprime = TMath::Sin(beamangle) * xprime + TMath::Cos(beamangle) * zprime;
    float etadoubleprime = TMath::ATanH(zdoubleprime/TMath::Sqrt(xdoubleprime*xdoubleprime+
     ydoubleprime*ydoubleprime+
     zdoubleprime*zdoubleprime));    
    
    if (ii==0) {
      etamin->SetPoint(ipart,phi,eta);
      etaminM->SetPoint(ipart,phi,etadoubleprime);
    } else {
      etamax->SetPoint(ipart,phi,eta);
      etamaxM->SetPoint(ipart,phi,etadoubleprime);
    }
    
  } // min/max 
  
  double eta0;
  double eta0M;
  double eta1;
  double eta1M;
  double phireal;
  etamin->GetPoint(ipart,phireal,eta0);
  etamax->GetPoint(ipart,phireal,eta1);
  etaminM->GetPoint(ipart,phireal,eta0M);
  etamaxM->GetPoint(ipart,phireal,eta1M);
  
  // integrate particles...
  particles->SetPoint( ipart,phi, hampt->Integral(hampt->FindBin(eta0),hampt->FindBin(eta1)));
  particlesM->SetPoint(ipart,phi, hampt->Integral(hampt->FindBin(eta0M),hampt->FindBin(eta1M))); 
  modification->SetPoint(ipart,phi,
   hampt->Integral(hampt->FindBin(eta0),hampt->FindBin(eta1)) / 
   hampt->Integral(hampt->FindBin(eta0M),hampt->FindBin(eta1M))); 
  gridgraph[idet][iz]->SetPoint(ipart,phi,
   hampt->Integral(hampt->FindBin(eta0),hampt->FindBin(eta1)) / 
   hampt->Integral(hampt->FindBin(eta0M),hampt->FindBin(eta1M))); 
  
      } // end loop over phi angles


    } // end loop over z-vertex values
  } // end loop over detector planes

  //==================================================

  TFile *fout = new TFile("fout_grid.root","RECREATE");
  for (int idet=0;idet<5;idet++) {
   for (int iz=0;iz<60;iz++) {
     gridgraph[idet][iz]->Write();
   }
 }
 fout->Close();
}