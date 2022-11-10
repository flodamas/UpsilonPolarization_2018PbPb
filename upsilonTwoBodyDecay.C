//
// This program will simulate a two-body decay
//
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TLorentzVector.h"

#include <iostream>
#include <cmath>
//
// upsilonTwoBodyDecay
// 
// Code to generate a two-body decay,
// and study the kinematics of the daughters.
// This 
//

void upsilonTwoBodyDecay() {

  TFile* file = new TFile("JPsi1STwoBodyDecayLargePhaseSpace.root","RECREATE","Upsilon 1S 2-Body Decay");  // File to store new nTuple
  //  Init ups 'tuple
  TString varlist = "upsM:upsRap:upsP:upsPt:upsPz:upsEta:upsPhi:upsCosTheta:eleCosThetaPrime:eleP:elePt:elePz:eleEta:elePhi:posP:posPt:posPz:posEta:posPhi";

  TNtuple* mTuple = new TNtuple("ups", "Generated Upsilon ntuple", varlist);
  
  TRandom3 rnd(0);//
  //float eleM = 0.511/1000.0;  // 0.511 MeV/c^2
  float eleM = 0.113;  // mu mass = 113 MeV/c^2
  double upsM =  0; // in case we need to set in in the loop
  //
  // Start loop to generate upsilons:
  //
  const int iterate = 1e6;
  //TF1* zBreitWigner = new TF1("zBreitWigner","TMath::BreitWigner(x,91.1876,2.474)",.3,200);
  for (int i = 0;i<iterate;i++){
    // For generating random masses
    upsM = 3.09688;
    //upsM =   9.46; //  10.02; 10.36; (30.*rnd.Rndm());
    //upsM = rnd.BreitWigner(91.1876,2.474); // Z-mass
    if (upsM<0.3) continue;
    
    if ((i%10000)==0) {
      cout << "i = " << i << " -------" << endl;
      cout << "upsM = " << upsM << endl;
    }
    //
    // Upsilons:
    //
    // Generate Upsilons with momentum according to
    // Flat pT (in the range 0-5 Gev/c )
    // Flat y  (in the range -1.5 - 1.5 )
    // Flat phi (in the range -pi - pi )
    //

    float upsPt = 30.*rnd.Rndm(); // should be flat in 0-10
    float upsRap = 8.0*rnd.Rndm()-4.;// 
    float upsPhi = TMath::TwoPi()*rnd.Rndm()-TMath::Pi();
    float upsMt = sqrt(upsM*upsM+upsPt*upsPt); // mt^2 = m^2 + pt^2 
    float upsPz = upsMt*sinh(upsRap); // pz = mt sinh (y)
    float upsE = upsMt*cosh(upsRap);
    TVector3 upsPvec(upsPt*cos(upsPhi),upsPt*sin(upsPhi),upsPz);
    float upsEta = upsPvec.PseudoRapidity();
    if (upsEta>1000) {
      cout << "Large eta!!" << endl;
      cout << "upsPt = " << upsPt << endl;
      cout << "upsPz = " << upsPz << endl;
    }
   //upsPvec vector is the Upsilon 3-momentum in the lab frame
    TLorentzVector upsVec(upsPvec,upsE); // Upsilon 4-mom in the lab frame
 
    //
    // Electron daughter decay in the Parent rest frame:
    //
    // Generating random phi (between -pi and pi)
    // Generate random cos(theta) flat between -1 and 1
    // the electron momentum vector in the rest frame
    // of the parent always has magnitude equal to M_ups/2
    //
    float eleCosThetaPrime = 2.0*rnd.Rndm()-1.0;
    float elePhiGenPrime = TMath::TwoPi()*rnd.Rndm()-TMath::Pi();
    
    // direction vector of the electron daughter with respect to parent direction, in parent rest frame
    TVector3 dir(cos(elePhiGenPrime)*sin(acos(eleCosThetaPrime)),sin(elePhiGenPrime)*sin(acos(eleCosThetaPrime)),eleCosThetaPrime);
    float eleE = upsM/2.; // the energy of each electron is just given by half of the parent mass
    float elePprime = sqrt(eleE*eleE-eleM*eleM); // p = sqrt (E^2 - m^2)  
    TVector3 elePprimeVec(elePprime*dir); // momentum vector of electron = magnitude x direction_vector
    
    // At this point elePprimeVec is given with respect to the parent direction (in other words, elePprime
    // is written as if the parent direction is the z-axis).
    // Since the parent direction is really given by upsPvec in the lab (and not the z-axis in the lab), 
    // we need to rotate the vector so that its coordinates are given in the lab 
    // (albeit still in the parent rest frame)
    // 
    elePprimeVec.RotateY(upsPvec.Theta());
    elePprimeVec.RotateZ(upsPvec.Phi());
    // Positrons thrown opposite electrons
    TVector3 posPprimeVec(-elePprimeVec.X(),-elePprimeVec.Y(),-elePprimeVec.Z());
    
    // check: The dot product (divided by the magnitudes)
    // of the electron vector and the upsilon vector should be equal to cosThetaPrime:
    //cout << "cosTheta* (from rotated vector) = " << elePprimeVec.Dot(upsPvec)/(elePprimeVec.Mag()*upsPvec.Mag()) << endl;
    //cout << "cosTheta* (as thrown)           = " << eleCosThetaPrime << endl;
    //
    // We are still in the rest frame of the Upsilon, need to boost into the lab frame
    //
    // First, make the 4-vectors of the electron and positron
    // we need the energy of the electron (and positron). These are the same, because the momentum
    // magnitudes of the 
    TLorentzVector eleVec(elePprimeVec,eleE);
    TLorentzVector posVec(posPprimeVec,eleE);   
    // Now, boost them into the Lab frame using the Upsilon boost velocity
    eleVec.Boost(upsVec.BoostVector());
    posVec.Boost(upsVec.BoostVector());

    //  eleVec & posVec are now in the lab frame
    
    TLorentzVector upsilon = eleVec + posVec; 
    // compare upsilon to upsVec (they should be identical, modulo machine precision)
    //cout << "eleVec + posVec  = " << upsilon.E() << ' ' << upsilon.Px() << ' ' << upsilon.Py() << ' ' << upsilon.Pz() << endl;
    //cout << "upsVec           = " << upsVec.E() << ' ' << upsVec.Px() << ' ' << upsVec.Py() << ' ' << upsVec.Pz() << endl;
    //  calculate MC opening angle:
    TVector3 eleLab3 = eleVec.Vect();
    TVector3 posLab3 = posVec.Vect();
    float upsilonCosTheta = cos(eleLab3.Angle(posLab3));
    
    float elePGen = eleVec.P(); 
    float posPGen = posVec.P(); 

    float elePtGen = eleVec.Perp(); 
    float posPtGen = posVec.Perp(); 

    float elePzGen = eleVec.Z();  
    float posPzGen = posVec.Z();
    
    float elePhiGen = eleVec.Phi();
    float posPhiGen = posVec.Phi();

    float eleEtaGen = eleVec.Eta();
    float posEtaGen = posVec.Eta();
    if (eleEtaGen<-1000) {
      cout << "Large eta electron!!" << endl;
      cout << "elePtGen = " << elePtGen << endl;
      cout << "elePzGen = " << elePzGen << endl;
    }
    if (posEtaGen<-1000) {
      cout << "Large eta positron!!" << endl;
      cout << "posPtGen = " << posPtGen << endl;
      cout << "posPzGen = " << posPzGen << endl;
    }
    
    float tuple[] = {
      static_cast<float>(upsM),
      upsRap,
      static_cast<float>(upsVec.P()),
      upsPt,
      upsPz,
      upsEta,
      upsPhi,
      upsilonCosTheta,
      eleCosThetaPrime,
      elePGen,
      elePtGen,
      elePzGen,
      eleEtaGen,
      elePhiGen,
      posPGen,
      posPtGen,
      posPzGen,
      posEtaGen,
      posPhiGen
    };
    mTuple->Fill(tuple);
  }// loop over generated upsilons
  
  mTuple->Write();
  file->Close();
  return;
}
