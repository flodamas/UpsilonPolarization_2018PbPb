// when the header is included several times, to avoid the redefinition error
#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

/// Helping functions returning the positive muon's coordinates in a given reference frame based on the TLorentzVectors in the lab frame
#include "TVector3.h"
#include "TLorentzVector.h"

// Lab to Helicity
TVector3 MuPlusVector_Helicity(const TLorentzVector QQLV_Lab, const TLorentzVector MuPlusLV_Lab) {
	// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	TLorentzVector MuPlusLV_QQRestFrame(MuPlusLV_Lab);

	//(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
	//(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
	MuPlusLV_QQRestFrame.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 MuPlusVec_Boosted = MuPlusLV_QQRestFrame.Vect();

	//Note: TVector3.Rotate() rotates the vectors, not the coordinates, so should rotate -phi and -theta

	MuPlusVec_Boosted.RotateZ(-QQVector_Lab.Phi());

	MuPlusVec_Boosted.RotateY(-QQVector_Lab.Theta());

	return MuPlusVec_Boosted;
}

// Lab to Collins-Soper
// requires the beam parameters
TVector3 MuPlusVector_CollinsSoper(const TLorentzVector QQLV_Lab, const TLorentzVector MuPlusLV_Lab) {
	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.02;                    //(Center of mass Energy per nucleon pair in TeV)
	double beamEnergy = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)

	// ******** HX to CS (rotation from HX frame to CS frame) ******** //
	// (1. Boost two beams to upsilon's rest frame)
	// (2. Rotate the coordinates)
	// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
	// (4. Calculate delta (angle btw ZHX and ZCS))

	// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
	TLorentzVector Beam1LV_Boosted(0., 0., beamEnergy, beamEnergy);
	TLorentzVector Beam2LV_Boosted(0., 0., -beamEnergy, beamEnergy); // mind the sign!!

	Beam1LV_Boosted.Boost(-QQLV_Lab.BoostVector());
	Beam2LV_Boosted.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 Beam1Vector_QQRestFrame(Beam1LV_Boosted.Vect());
	TVector3 Beam2Vector_QQRestFrame(Beam2LV_Boosted.Vect());

	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	Beam1Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam1Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	Beam2Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam2Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	// ******** Calculate the angle between z_HX and z_CS ******** //
	TVector3 ZHXunitVec(0, 0, 1.);                                                  //(define z_HX unit vector)
	double Angle_B1ZHX = Beam1Vector_QQRestFrame.Angle(ZHXunitVec);                //(angle between beam1 and z_HX)
	double Angle_B2ZHX = Beam2Vector_QQRestFrame.Angle(-ZHXunitVec);               //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
	double Angle_B1miB2 = Beam1Vector_QQRestFrame.Angle(-Beam2Vector_QQRestFrame); //(angle between beam1 and -beam2)

	double delta = 0; //(define and initialize the angle between z_HX and z_CS)

	// Maths for calculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
	if (Angle_B1ZHX > Angle_B2ZHX)
		delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
	else if (Angle_B1ZHX < Angle_B2ZHX)
		delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
	else
		std::cout << "beam1PvecBoosted.Pz() = 0?" << std::endl;

	// ******** Rotate the coordinates along the y-axis by the angle between z_HX and z_CS ******** //
	TVector3 MuPlusVec_CS(MuPlusVector_Helicity(QQLV_Lab, MuPlusLV_Lab));

	MuPlusVec_CS.RotateY(delta);

	return MuPlusVec_CS;
}

#endif