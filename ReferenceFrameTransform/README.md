# Reference Frame Transformation Codes
This is the code to transform daughter muons in the lab frame to the Helicity, and the Collis-Soper frame.
First, transform daughter muons in the lab frame to the Helicity frame.
Second, transform daughter muons in the Helicity frame to the Collis-Soper frame using the angle between two reference frames.

Oniatree_MC_numEvent1000_1342.root : Input Oniatree File
UpsilonRefFrame2.C : Used the same method in the sample code (upsilonTwoBodyDecay.C). Boosted muons to the upsilon's rest frame and rotated the coordinate by #phi and #theth of upsilon in the lab frame.
UpsilonRefFrame3.C : Rotated the coordinate by #phi and #theth of upsilon in the lab frame and boosted muons to the upsilon's rest frame.

UpsilonRefFrame2.C and UpsilonRefFrame3.C give the same results.
