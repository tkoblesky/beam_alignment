#include <TLorentzVector.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TColor.h>

#include <iostream>

using namespace std;

// a macro written in C++ with objects in the ROOT framework
// meant to be run with ROOT.

// This macro creates a grid of possible 'yellow' and 'blue' beams orientations
// and computes the Lorentz boost and rotation necessary in order to bring
// the two beams into a colinear reference frame
// once the colinear boost and rotation is found, the detector (BBC) elements are rotated into place
// the detector elements are represented as low momentum particles.

static const int verbosity = 1;

// this function takes in 2 lorentz vectors (energy, momentum) (each representing one of the beams)
// and a vector representing the detector eleements
// this function returns the rotation angle necessary to fix the beams
float boost_and_rotate(TLorentzVector & vec, TLorentzVector input1, TLorentzVector input2)
{
  // takes in 3 vectors in the LAB FRAME

  TVector3 z(0,0,1);//this represents the z-axis (ideal beam axis)

  if(verbosity > 0)// print to screen the 
  {
    cout << "-----------------BEFORE BOOST-----------------" << endl;
    cout << "px1 = " << input1.Px() << ", py1 = " << input1.Py() << ", pz1 = " << input1.Pz() << endl;
    cout << "px2 = " << input2.Px() << ", py2 = " << input2.Py() << ", pz2 = " << input2.Pz() << endl;
    cout << "Angle between beams = " << input1.Angle(input2.Vect()) << endl;
    cout << "Angle between beam1 and z axis = " << input1.Angle(z) << endl;
    cout << "Angle between beam2 and z axis = " << input2.Angle(z) << endl << endl;
  }

  TLorentzVector cms = input1 + input2;
  input1.Boost(-cms.BoostVector());//blue
  input2.Boost(-cms.BoostVector());//yellow
  vec.Boost(-cms.BoostVector());//boost the opposite direction to LAB FRAME
  cms.Boost(-cms.BoostVector());

  if(verbosity > 0)
  {
    cout << "CHECKING THAT TOTAL MOMENTUM IN BOOSTED CoM FRAME BE IDENTICALLY ZERO..." << endl;
    cout << "px_com = " << cms.Px() << ", py_com = " << cms.Py() << ", pz_com = " << cms.Pz() << endl << endl;
  }
  //check to see if input2.Angle(z) is the complement of my default angle

  double rotAngleY = -input1.Angle(z);
  if(verbosity > 0 )
  {
    cout<<"-input1.Angle(z): "<<rotAngleY<<" input2.Angle(z): "<<input2.Angle(z)<<endl;
  }

  input1.RotateY(rotAngleY);
  input2.RotateY(rotAngleY);
  vec.RotateY(rotAngleY);

  if(verbosity > 0)
  {
    cout << "----------------- AFTER BOOST + ROTATION -----------------" << endl;
    cout << "blue: px1' = " << input1.Px() << ", py1' = " << input1.Py() << ", pz1' = " << input1.Pz() << endl;
    cout << "yellow: px2' = " << input2.Px() << ", py2' = " << input2.Py() << ", pz2' = " << input2.Pz() << endl;
    cout << "vec: px3' = " << vec.Px() << ", py3' = " << vec.Py() << ", pz3' = " << vec.Pz() << endl;
  }

  //if rotation was the wrong way (minus sign error), go twice the opposite way
  if(TMath::Abs(input1.Px()) >  0.00001 )
  {
    input1.RotateY(-2*rotAngleY);
    input2.RotateY(-2*rotAngleY);
    vec.RotateY(-2*rotAngleY);

    if(verbosity > 0)
    {
      cout << "----------------- AFTER BOOST + INVERSE ROTATION-----------------" << endl;
      cout << "Rotation angle: "<<-2*rotAngleY<<endl;
      cout << "px1' = " << input1.Px() << ", py1' = " << input1.Py() << ", pz1' = " << input1.Pz() << endl;
      cout << "px2' = " << input2.Px() << ", py2' = " << input2.Py() << ", pz2' = " << input2.Pz() << endl;
      cout << "px3' = " << vec.Px() << ", py3' = " << vec.Py() << ", pz3' = " << vec.Pz() << endl;
    }

    return -2*rotAngleY;
  }
  return rotAngleY;
}

// function which scans over different beam angles and makes histograms of the results
// calls boost_and_rotate
void scan_boost_angle()
{

  //TStyle optimization
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.3f");

  double pi = TMath::Pi();

  // utility histograms
  TH2D *average_phi_diff = new TH2D("average_phi_diff","average_phi_diff",7,-0.0024-0.00005,0.0024-0.00005,9,-0.0048-0.00005,0.0048-0.00005);
  TH1D * back_rap_phi_dist = new TH1D("back_rap_phi_dist","back_rap_phi_dist",140,-pi,pi);
  TH1D * for_rap_phi_dist = new TH1D("for_rap_phi_dist","for_rap_phi_dist",140,-pi,pi);
  TH2D *rot_angles = new TH2D("rot_angles","rot_angles",7,-0.0024-0.00005,0.0024-0.00005,9,-0.0048-0.00005,0.0048-0.00005);

  int npart = 1;//number of particles

  // hard coded beam orientations in mRad
  string blue_labels[7] = {"4.0","3.2","2.4","1.6","0.8","0.0","-0.8"};
  string yellow_labels[9] = {"8.4","7.2","6.0","4.8","3.6","2.4","1.2","0.0","-1.2"};

  for(int ibin = 1; ibin < average_phi_diff->GetNbinsX()+1; ibin++)
  {
  	average_phi_diff->GetXaxis()->SetBinLabel(ibin,blue_labels[ibin - 1].c_str());
    rot_angles->GetXaxis()->SetBinLabel(ibin,blue_labels[ibin - 1].c_str());
  }

  for(int jbin = 1; jbin < average_phi_diff->GetNbinsY()+1; jbin++)
  {
  	average_phi_diff->GetYaxis()->SetBinLabel(jbin,yellow_labels[jbin - 1].c_str());
    rot_angles->GetYaxis()->SetBinLabel(jbin,yellow_labels[jbin - 1].c_str());
  }

  TF1 *gaus = new TF1("gaus","TMath::Exp(-x*x/([0]*[0]))",-7,7);
  gaus->SetParameter(0,6);

  //loop over beam angels
  for(int iangle2 = 0; iangle2 < 9; iangle2++)
  {
   float yellow_angle = TMath::Pi()+0.0036+0.0048-iangle2*2*0.0048/8;
   float yellow_px = 100*TMath::Sin(yellow_angle);
   float yellow_py = 0.0;
   float yellow_pz = 100*TMath::Cos(yellow_angle);
   float proton_mass = 0.938;
   float yellow_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
   TLorentzVector yellow_beam(yellow_px,yellow_py,yellow_pz,yellow_energy);

   for(int iangle1 = 0; iangle1 < 7; iangle1++)
   {
    float blue_angle = 0.0016+0.0024-iangle1*2*0.0024/6;
    float blue_px = 100*TMath::Sin(blue_angle);
    float blue_py = 0.0;
    float blue_pz = 100*TMath::Cos(blue_angle);
    //float proton_mass = 0.938;
    float blue_energy = TMath::Sqrt(100*100+proton_mass*proton_mass);
    TLorentzVector blue_beam(blue_px,blue_py,blue_pz,blue_energy);

    for(int ipart = 0; ipart < npart; ipart++)
    {
      float rand_phi = TMath::Pi()*gRandom->Uniform(-1,1);
      float rand_eta = gaus->GetRandom();//gRandom->Uniform(-3,-1);//gaussian centered at zero with sigma of 6
      float particle_pT = 0.250;//in GeV
      float px = particle_pT*TMath::Cos(rand_phi);
      float py = particle_pT*TMath::Sin(rand_phi);
      float pz = particle_pT*TMath::SinH(rand_eta);
      float energy = TMath::Sqrt(px*px+py*py+pz*pz+0.1396*0.1396);

      TLorentzVector pmt_vec(px,py,pz,energy);
      if(verbosity > 0)
      {
        cout<<"------------------------------------------------------------"<<endl;
        cout<<"for angle configuration blue: "<<blue_angle<<" yellow-pi: "<<yellow_angle-TMath::Pi()<<endl;
      }

      rot_angles->SetBinContent(iangle1+1,iangle2+1,boost_and_rotate(pmt_vec,blue_beam,yellow_beam));
      if(verbosity > 0)
        cout<<endl<<endl;
      float phi = TMath::ATan2(pmt_vec.Py(),pmt_vec.Px());
      float eta = TMath::ASinH(pmt_vec.Pz()/TMath::Sqrt(pmt_vec.Px()*pmt_vec.Px()+pmt_vec.Py()*pmt_vec.Py()));
      if(eta > -4 && eta < -3)
        back_rap_phi_dist->Fill(phi);
      else if(eta > 3 && eta < 4)
        for_rap_phi_dist->Fill(phi);
      if(ipart==0)
        average_phi_diff->SetBinContent(iangle1+1,iangle2+1,TMath::Abs(phi-rand_phi));
      else 
        average_phi_diff->SetBinContent(iangle1+1,iangle2+1,TMath::Abs(phi-rand_phi) + average_phi_diff->GetBinContent(iangle1+1,iangle2+1));
    }// end loop ipart
  }// end loop iangle2
}// end loop iangle1

average_phi_diff->Scale(1.0/npart);//normalize by number of particles

//draw results on canvas
TCanvas *c1 = new TCanvas("c1","c1",800,800);
c1->SetRightMargin(0.15);

rot_angles->Draw("colz");
rot_angles->DrawCopy("textsame");

TCanvas *c0 = new TCanvas("c0","c0",800,800);
c0->SetRightMargin(0.15);
average_phi_diff->SetXTitle("blue angles (mRad)");
average_phi_diff->SetYTitle("yellow angles (mRad)");
average_phi_diff->Draw("colz");
average_phi_diff->DrawCopy("textsame");

TLegend *leg = new TLegend(0.15,0.75,0.3,0.9);
TCanvas *c2 = new TCanvas("c2","c2",800,800);
for_rap_phi_dist->SetLineColor(kRed);
back_rap_phi_dist->SetLineColor(kBlue);
for_rap_phi_dist->Draw();
back_rap_phi_dist->Draw("same");

}
