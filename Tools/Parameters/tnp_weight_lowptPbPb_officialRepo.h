#ifndef tnp_weight_lowptPbPb_h
#define tnp_weight_lowptPbPb_h
#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// +++++++++++++++++++++++++++++++++++++++
// - Trigger: (tnp_weight_trg_pbpb)
//   * filterId = 0: Jpsi L2 filter
//   * filterId = 1: Jpsi L3 filter
//   * filterId = 2: Upsi L2 filter
//   * filterId = 3: Upsi L3 filter
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - MuID: (tnp_weight_muid_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - Inner tracking: (tnp_weight_trk_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// +++++++++++++++++++++++++++++++++++++++

double tnp_weight_muid_pbpb(double pt, double eta, int idx=0);
double tnp_weight_trk_pbpb(double eta, int idx=0);
double tnp_weight_trg_pbpb(double pt, double eta, int filterId=0,int idx=0);

///////////////////////////////////////////////////
//             M u I D    P b P b                //
///////////////////////////////////////////////////

double tnp_weight_muid_pbpb(double pt, double eta, int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  // SF for 0 < |eta| < 1.2
  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    if (x >= 3.5 && x <4) {num = 0.946159; den = 0.951487; statUp = 0.0037907; statDown = 0.00389305;}
    else if (x >= 4 && x <4.5) {num = 0.973653; den = 0.975115; statUp = 0.00284919; statDown = 0.00295999;}
    else if (x >= 4.5 && x <5) {num = 0.984681; den = 0.983948; statUp = 0.00264355; statDown = 0.00278459;}
    else if (x >= 5 && x <5.5) {num = 0.985002; den = 0.989214; statUp = 0.00286271; statDown = 0.00302699;}
    else if (x >= 5.5 && x <6.5) {num = 0.988848; den = 0.990738; statUp = 0.00200687; statDown = 0.00212372;}
    else if (x >= 6.5 && x <8) {num = 0.991114; den = 0.992378; statUp = 0.00174382; statDown = 0.00186571;}
    else if (x >= 8 && x <10.5) {num = 0.986281; den = 0.989496; statUp = 0.00215883; statDown = 0.00230715;}
    else if (x >= 10.5 && x <14) {num = 0.976775; den = 0.968327; statUp = 0.00361383; statDown = 0.00385533;}
    else if (x >= 14 && x <18) {num = 0.980589; den = 0.977123; statUp = 0.00486848; statDown = 0.00537745;}
    else {num = 0.984897; den = 0.982298; statUp = 0.00550232; statDown = 0.00622724;}
  }
  // SF for 1.2 < |eta| < 1.8
  if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    if (x >= 2.07 && x <3) {num = 0.981569; den = 0.981282; statUp = 0.00677595; statDown = 0.00689087;}
    else if (x >= 3 && x <3.5) {num = 0.986502; den = 0.985373; statUp = 0.00410813; statDown = 0.00425456;}
    else if (x >= 3.5 && x <4) {num = 0.98077; den = 0.989287; statUp = 0.00400325; statDown = 0.00400325;}
    else if (x >= 4 && x <4.5) {num = 0.986286; den = 0.992139; statUp = 0.00359094; statDown = 0.00382453;}
    else if (x >= 4.5 && x <5) {num = 0.992842; den = 0.992784; statUp = 0.00351371; statDown = 0.00380977;}
    else if (x >= 5 && x <6) {num = 0.989169; den = 0.995025; statUp = 0.00291819; statDown = 0.00314488;}
    else if (x >= 6 && x <7.5) {num = 0.994996; den = 0.996835; statUp = 0.00241413; statDown = 0.00267229;}
    else if (x >= 7.5 && x <10) {num = 0.998548; den = 0.99638; statUp = 0.00145236; statDown = 0.00254248;}
    else if (x >= 10 && x <15) {num = 0.979249; den = 0.971019; statUp = 0.00526182; statDown = 0.0057076;}
    else {num = 0.994124; den = 0.98437; statUp = 0.0058759; statDown = 0.007951;}
  }
  // SF for 1.8 < |eta| < 2.1
  if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    if (x >= 1.5 && x <2.5) {num = 0.970228; den = 0.953223; statUp = 0.0110656; statDown = 0;}
    else if (x >= 2.5 && x <3) {num = 0.987297; den = 0.983762; statUp = 0.00593338; statDown = 0.00622044;}
    else if (x >= 3 && x <3.5) {num = 0.991708; den = 0.991014; statUp = 0.00482514; statDown = 0.00446851;}
    else if (x >= 3.5 && x <4) {num = 0.982788; den = 0.995327; statUp = 0.00481342; statDown = 0.00528113;}
    else if (x >= 4 && x <4.5) {num = 0.99717; den = 0.997323; statUp = 0.00282996; statDown = 0.00378255;}
    else if (x >= 4.5 && x <5.5) {num = 0.997013; den = 0.998097; statUp = 0.0022717; statDown = 0;}
    else if (x >= 5.5 && x <7) {num = 1; den = 0.998592; statUp = 1.50776e-07; statDown = 0.00122323;}
    else if (x >= 7 && x <9) {num = 0.996501; den = 0.999317; statUp = 0.00297386; statDown = 0.00362572;}
    else if (x >= 9 && x <12) {num = 0.997501; den = 0.996859; statUp = 0.00249868; statDown = 0.00478656;}
    else {num = 0.989345; den = 0.997618; statUp = 0.0051521; statDown = 0.00686294;}
  }
  // SF for 2.1 < |eta| < 2.4
  if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    if (x >= 1.5 && x <2.2) {num = 0.901632; den = 0.910903; statUp = 0.0172668; statDown = 0;}
    else if (x >= 2.2 && x <2.7) {num = 0.939984; den = 0.964945; statUp = 0.0107612; statDown = 0.0109264;}
    else if (x >= 2.7 && x <3.2) {num = 0.973211; den = 0.978488; statUp = 0.00846784; statDown = 0.00879021;}
    else if (x >= 3.2 && x <3.7) {num = 0.977562; den = 0.988794; statUp = 0.006863; statDown = 0.00730481;}
    else if (x >= 3.7 && x <4.7) {num = 0.996762; den = 0.993708; statUp = 0.00323819; statDown = 0.00498314;}
    else if (x >= 4.7 && x <8) {num = 0.996858; den = 0.997702; statUp = 0.00255714; statDown = 0.00289927;}
    else if (x >= 8 && x <11) {num = 1; den = 0.999093; statUp = 3.22096e-09; statDown = 0.00251078;}
    else if (x >= 11 && x <14) {num = 0.978749; den = 0.997809; statUp = 0.0100421; statDown = 0.0137177;}
    else {num = 1; den = 0.995974; statUp = 1.46346e-07; statDown = 0.00540197;}
  }

  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    // syst uncertainties
    if (x >= 3.5 && x < 4) syst = 0.000799144;
    else if (x >= 4 && x < 4.5) syst = 0.000316587;
    else if (x >= 4.5 && x < 5) syst = 0.000900099;
    else if (x >= 5 && x < 5.5) syst = 0.00110215;
    else if (x >= 5.5 && x < 6.5) syst = 0.000271437;
    else if (x >= 6.5 && x < 8) syst = 0.000512112;
    else if (x >= 8 && x < 10.5) syst = 0.00123065;
    else if (x >= 10.5 && x < 14) syst = 0.00263086;
    else if (x >= 14 && x < 18) syst = 0.00137269;
    else syst = 0.00370496;
  }
  else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    // syst uncertainties
    if (x >= 2.07 && x < 3) syst = 0.00212907;
    else if (x >= 3 && x < 3.5) syst = 0.00176869;
    else if (x >= 3.5 && x < 4) syst = 0.0018028;
    else if (x >= 4 && x < 4.5) syst = 0.00161887;
    else if (x >= 4.5 && x < 5) syst = 0.00177587;
    else if (x >= 5 && x < 6) syst = 0.000768592;
    else if (x >= 6 && x < 7.5) syst = 0.00104735;
    else if (x >= 7.5 && x < 10) syst = 0.00121085;
    else if (x >= 10 && x < 15) syst = 0.00181423;
    else syst = 0.00257806;
  }
  else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    // syst uncertainties
    if (x >= 1.5 && x < 2.5) syst = 0.00334452;
    else if (x >= 2.5 && x < 3) syst = 0.00116118;
    else if (x >= 3 && x < 3.5) syst = 0.00483524;
    else if (x >= 3.5 && x < 4) syst = 0.00108025;
    else if (x >= 4 && x < 4.5) syst = 0.00282868;
    else if (x >= 4.5 && x < 5.5) syst = 0.000630971;
    else if (x >= 5.5 && x < 7) syst = 7.68922e-07;
    else if (x >= 7 && x < 9) syst = 0.00185991;
    else if (x >= 9 && x < 12) syst = 0.000838922;
    else syst = 0.000763013;
  }
  else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    // syst uncertainties
    if (x >= 1.5 && x < 2.2) syst = 0.00406598;
    else if (x >= 2.2 && x < 2.7) syst = 0.00430655;
    else if (x >= 2.7 && x < 3.2) syst = 0.00328867;
    else if (x >= 3.2 && x < 3.7) syst = 0.00327706;
    else if (x >= 3.7 && x < 4.7) syst = 0.0018308;
    else if (x >= 4.7 && x < 8) syst = 0.00193408;
    else if (x >= 8 && x < 11) syst = 2.29509e-07;
    else if (x >= 11 && x < 14) syst = 0.00129097;
    else syst = 1.46346e-07;
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  else if (idx == -2) syst_factor = -1*syst;
  else if (idx == +1) stat_factor = statUp;
  else if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//              T R G     P b P b                //
///////////////////////////////////////////////////

double tnp_weight_trg_pbpb(double pt, double eta, int filterId,int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (filterId==0) { //L2 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.66246; den = 0.617044; statUp = 0.00759822; statDown = 0.00764213;}
      else if (x >= 4 && x <4.5) {num = 0.859542; den = 0.838993; statUp = 0.00554625; statDown = 0.00565142;}
      else if (x >= 4.5 && x <5) {num = 0.916898; den = 0.90054; statUp = 0.00471306; statDown = 0.00485668;}
      else if (x >= 5 && x <5.5) {num = 0.934879; den = 0.92338; statUp = 0.00458324; statDown = 0.0047695;}
      else if (x >= 5.5 && x <6.5) {num = 0.947552; den = 0.940172; statUp = 0.00344535; statDown = 0.00357327;}
      else if (x >= 6.5 && x <8) {num = 0.94982; den = 0.956634; statUp = 0.00341385; statDown = 0.00355042;}
      else if (x >= 8 && x <10.5) {num = 0.95068; den = 0.962863; statUp = 0.00364888; statDown = 0.00380937;}
      else if (x >= 10.5 && x <14) {num = 0.947627; den = 0.963749; statUp = 0.00479744; statDown = 0.00506693;}
      else if (x >= 14 && x <18) {num = 0.939736; den = 0.965572; statUp = 0.00768532; statDown = 0.00828185;}
      else {num = 0.943; den = 0.957231; statUp = 0.00880706; statDown = 0.00961831;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.07 && x <3) {num = 0.672197; den = 0.634105; statUp = 0.0116105; statDown = 0.0116585;}
      else if (x >= 3 && x <3.5) {num = 0.830014; den = 0.793033; statUp = 0.00781546; statDown = 0.00793036;}
      else if (x >= 3.5 && x <4) {num = 0.901627; den = 0.867368; statUp = 0.00658615; statDown = 0.00658615;}
      else if (x >= 4 && x <4.5) {num = 0.918707; den = 0.916768; statUp = 0.00628767; statDown = 0.00654209;}
      else if (x >= 4.5 && x <5) {num = 0.935159; den = 0.935819; statUp = 0.0064133; statDown = 0.00675361;}
      else if (x >= 5 && x <6) {num = 0.940538; den = 0.949225; statUp = 0.0052386; statDown = 0.00549113;}
      else if (x >= 6 && x <7.5) {num = 0.943354; den = 0.957099; statUp = 0.00541232; statDown = 0.0057046;}
      else if (x >= 7.5 && x <10) {num = 0.952833; den = 0.959401; statUp = 0.00555608; statDown = 0.00593069;}
      else if (x >= 10 && x <15) {num = 0.950014; den = 0.960094; statUp = 0.00676702; statDown = 0.0073103;}
      else {num = 0.914564; den = 0.953433; statUp = 0.0132629; statDown = 0.0146698;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.5 && x <2.5) {num = 0.69574; den = 0.673526; statUp = 0.0132158; statDown = 0.0132339;}
      else if (x >= 2.5 && x <3) {num = 0.881983; den = 0.892382; statUp = 0.0102195; statDown = 0.0104922;}
      else if (x >= 3 && x <3.5) {num = 0.908943; den = 0.929761; statUp = 0.00966169; statDown = 0.0100646;}
      else if (x >= 3.5 && x <4) {num = 0.920464; den = 0.94472; statUp = 0.00960848; statDown = 0.0100994;}
      else if (x >= 4 && x <4.5) {num = 0.907279; den = 0.948185; statUp = 0.0108863; statDown = 0.0115854;}
      else if (x >= 4.5 && x <5.5) {num = 0.929741; den = 0.951593; statUp = 0.00797835; statDown = 0.00841839;}
      else if (x >= 5.5 && x <6.5) {num = 0.892463; den = 0.953392; statUp = 0.0118058; statDown = 0.0123937;}
      else if (x >= 6.5 && x <8) {num = 0.917335; den = 0.947148; statUp = 0.0100574; statDown = 0.0107488;}
      else if (x >= 8 && x <9.5) {num = 0.912892; den = 0.941109; statUp = 0.0137674; statDown = 0.0149277;}
      else if (x >= 9.5 && x <13) {num = 0.915007; den = 0.941082; statUp = 0.013898; statDown = 0.0151274;}
      else {num = 0.904239; den = 0.925954; statUp = 0.0221633; statDown = 0.0243071;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <2.2) {num = 0.749861; den = 0.728213; statUp = 0.01662; statDown = 0.0165872;}
      else if (x >= 2.2 && x <2.7) {num = 0.854468; den = 0.865648; statUp = 0.0140236; statDown = 0.0142579;}
      else if (x >= 2.7 && x <3.2) {num = 0.878285; den = 0.893429; statUp = 0.0129941; statDown = 0.0133672;}
      else if (x >= 3.2 && x <3.7) {num = 0.8799; den = 0.909; statUp = 0.0120868; statDown = 0.0125586;}
      else if (x >= 3.7 && x <4.7) {num = 0.887508; den = 0.911976; statUp = 0.0108445; statDown = 0.0105222;}
      else if (x >= 4.7 && x <6.5) {num = 0.913489; den = 0.926877; statUp = 0.00954989; statDown = 0.0100281;}
      else if (x >= 6.5 && x <8.5) {num = 0.900517; den = 0.922227; statUp = 0.0144257; statDown = 0.0153605;}
      else if (x >= 8.5 && x <11) {num = 0.868413; den = 0.91872; statUp = 0.0200071; statDown = 0.0215637;}
      else {num = 0.889934; den = 0.912376; statUp = 0.0231784; statDown = 0.0255139;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00124485;
      else if (x >= 4 && x < 4.5) syst = 0.00103897;
      else if (x >= 4.5 && x < 5) syst = 0.00129541;
      else if (x >= 5 && x < 5.5) syst = 0.000789457;
      else if (x >= 5.5 && x < 6.5) syst = 0.000603473;
      else if (x >= 6.5 && x < 8) syst = 0.000282966;
      else if (x >= 8 && x < 10.5) syst = 0.000187072;
      else if (x >= 10.5 && x < 14) syst = 0.000282897;
      else if (x >= 14 && x < 18) syst = 0.00120202;
      else syst = 0.00323389;
    }
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.07 && x < 3) syst = 0.0057738;
      else if (x >= 3 && x < 3.5) syst = 0.00516636;
      else if (x >= 3.5 && x < 4) syst = 0.00227344;
      else if (x >= 4 && x < 4.5) syst = 0.00228932;
      else if (x >= 4.5 && x < 5) syst = 0.00180741;
      else if (x >= 5 && x < 6) syst = 0.00169068;
      else if (x >= 6 && x < 7.5) syst = 0.00053155;
      else if (x >= 7.5 && x < 10) syst = 0.00184744;
      else if (x >= 10 && x < 15) syst = 0.000513966;
      else syst = 0.00398062;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.5 && x < 2.5) syst = 0.00483512;
      else if (x >= 2.5 && x < 3) syst = 0.00268291;
      else if (x >= 3 && x < 3.5) syst = 0.00462034;
      else if (x >= 3.5 && x < 4) syst = 0.00501381;
      else if (x >= 4 && x < 4.5) syst = 0.0104532;
      else if (x >= 4.5 && x < 5.5) syst = 0.00248333;
      else if (x >= 5.5 && x < 6.5) syst = 0.00505641;
      else if (x >= 6.5 && x < 8) syst = 0.00309659;
      else if (x >= 8 && x < 9.5) syst = 0.00373394;
      else if (x >= 9.5 && x < 13) syst = 0.00272605;
      else syst = 0.0174777;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.5 && x < 2.2) syst = 0.0173906;
      else if (x >= 2.2 && x < 2.7) syst = 0.00738026;
      else if (x >= 2.7 && x < 3.2) syst = 0.00966916;
      else if (x >= 3.2 && x < 3.7) syst = 0.00539735;
      else if (x >= 3.7 && x < 4.7) syst = 0.0124751;
      else if (x >= 4.7 && x < 6.5) syst = 0.00722603;
      else if (x >= 6.5 && x < 8.5) syst = 0.00296615;
      else if (x >= 8.5 && x < 11) syst = 0.00969476;
      else syst = 0.0145531;
    }
  }
  else if (filterId==1) { //L3 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0977768; den = 0.0718102; statUp = 0.00477386; statDown = 0.00465538;}
      else if (x >= 4 && x <4.5) {num = 0.310189; den = 0.236612; statUp = 0.00733234; statDown = 0.00727546;}
      else if (x >= 4.5 && x <5) {num = 0.497824; den = 0.427955; statUp = 0.00854784; statDown = 0.00854347;}
      else if (x >= 5 && x <5.5) {num = 0.647029; den = 0.570059; statUp = 0.00899919; statDown = 0.00905681;}
      else if (x >= 5.5 && x <6.5) {num = 0.718238; den = 0.667318; statUp = 0.0069276; statDown = 0.00699004;}
      else if (x >= 6.5 && x <8) {num = 0.770636; den = 0.738796; statUp = 0.00665117; statDown = 0.00673002;}
      else if (x >= 8 && x <10.5) {num = 0.791377; den = 0.780799; statUp = 0.0068868; statDown = 0.00698942;}
      else if (x >= 10.5 && x <14) {num = 0.828454; den = 0.817398; statUp = 0.00832993; statDown = 0.00852906;}
      else if (x >= 14 && x <18) {num = 0.798313; den = 0.82568; statUp = 0.0132365; statDown = 0.0136444;}
      else {num = 0.844425; den = 0.840963; statUp = 0.0140504; statDown = 0.0146774;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.07 && x <3) {num = 0.312297; den = 0.282025; statUp = 0.0108985; statDown = 0.0107675;}
      else if (x >= 3 && x <3.5) {num = 0.428221; den = 0.427322; statUp = 0.00996613; statDown = 0.00990255;}
      else if (x >= 3.5 && x <4) {num = 0.545302; den = 0.523792; statUp = 0.0106388; statDown = 0.0106824;}
      else if (x >= 4 && x <4.5) {num = 0.590907; den = 0.59983; statUp = 0.011421; statDown = 0.0114593;}
      else if (x >= 4.5 && x <5) {num = 0.642516; den = 0.645172; statUp = 0.0127105; statDown = 0.0128095;}
      else if (x >= 5 && x <6) {num = 0.67265; den = 0.678719; statUp = 0.0105191; statDown = 0.0106482;}
      else if (x >= 6 && x <7.5) {num = 0.685196; den = 0.718066; statUp = 0.0110631; statDown = 0.0111844;}
      else if (x >= 7.5 && x <10) {num = 0.745909; den = 0.754042; statUp = 0.0115011; statDown = 0.0117057;}
      else if (x >= 10 && x <15) {num = 0.774291; den = 0.802324; statUp = 0.0133219; statDown = 0.0136598;}
      else {num = 0.766916; den = 0.825355; statUp = 0.0213965; statDown = 0.0220731;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.5 && x <2.5) {num = 0.127351; den = 0.116205; statUp = 0.00812389; statDown = 0.00790635;}
      else if (x >= 2.5 && x <3) {num = 0.469654; den = 0.461475; statUp = 0.0151839; statDown = 0.0150975;}
      else if (x >= 3 && x <3.5) {num = 0.66052; den = 0.657218; statUp = 0.0159007; statDown = 0.0159966;}
      else if (x >= 3.5 && x <4) {num = 0.711318; den = 0.716352; statUp = 0.0159382; statDown = 0.0161569;}
      else if (x >= 4 && x <4.5) {num = 0.72707; den = 0.743816; statUp = 0.0174651; statDown = 0.0178183;}
      else if (x >= 4.5 && x <5.5) {num = 0.760903; den = 0.783352; statUp = 0.0132749; statDown = 0.0135294;}
      else if (x >= 5.5 && x <6.5) {num = 0.76791; den = 0.818708; statUp = 0.0164832; statDown = 0.0169401;}
      else if (x >= 6.5 && x <8) {num = 0.810282; den = 0.834464; statUp = 0.0146744; statDown = 0.015183;}
      else if (x >= 8 && x <9.5) {num = 0.816841; den = 0.841405; statUp = 0.019156; statDown = 0.0201825;}
      else if (x >= 9.5 && x <13) {num = 0.830536; den = 0.857282; statUp = 0.0191548; statDown = 0.0201111;}
      else {num = 0.865685; den = 0.871118; statUp = 0.0261773; statDown = 0.0280197;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.5 && x <2.2) {num = 0.0687152; den = 0.0653009; statUp = 0.00763497; statDown = 0.0073555;}
      else if (x >= 2.2 && x <2.7) {num = 0.267663; den = 0.247036; statUp = 0.0154354; statDown = 0.0150849;}
      else if (x >= 2.7 && x <3.2) {num = 0.417871; den = 0.405866; statUp = 0.0183278; statDown = 0.0180847;}
      else if (x >= 3.2 && x <3.7) {num = 0.491052; den = 0.496108; statUp = 0.0185629; statDown = 0.018464;}
      else if (x >= 3.7 && x <4.7) {num = 0.57329; den = 0.569608; statUp = 0.0170203; statDown = 0.0170343;}
      else if (x >= 4.7 && x <6.5) {num = 0.676741; den = 0.652632; statUp = 0.0161595; statDown = 0.0163126;}
      else if (x >= 6.5 && x <8.5) {num = 0.717631; den = 0.71118; statUp = 0.0218543; statDown = 0.0223093;}
      else if (x >= 8.5 && x <11) {num = 0.713947; den = 0.744621; statUp = 0.0278094; statDown = 0.0286323;}
      else {num = 0.748393; den = 0.778517; statUp = 0.0328678; statDown = 0.0342111;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000736946;
      else if (x >= 4 && x < 4.5) syst = 0.00117449;
      else if (x >= 4.5 && x < 5) syst = 0.0031303;
      else if (x >= 5 && x < 5.5) syst = 0.00296318;
      else if (x >= 5.5 && x < 6.5) syst = 0.000684434;
      else if (x >= 6.5 && x < 8) syst = 0.000814384;
      else if (x >= 8 && x < 10.5) syst = 0.000744841;
      else if (x >= 10.5 && x < 14) syst = 0.000675376;
      else if (x >= 14 && x < 18) syst = 0.00210979;
      else syst = 0.00279935;
    }
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.07 && x < 3) syst = 0.00304772;
      else if (x >= 3 && x < 3.5) syst = 0.00425846;
      else if (x >= 3.5 && x < 4) syst = 0.00447385;
      else if (x >= 4 && x < 4.5) syst = 0.001133;
      else if (x >= 4.5 && x < 5) syst = 0.00368832;
      else if (x >= 5 && x < 6) syst = 0.00158046;
      else if (x >= 6 && x < 7.5) syst = 0.00182771;
      else if (x >= 7.5 && x < 10) syst = 0.00308072;
      else if (x >= 10 && x < 15) syst = 0.00292029;
      else syst = 0.00529299;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.5 && x < 2.5) syst = 0.00136009;
      else if (x >= 2.5 && x < 3) syst = 0.0038152;
      else if (x >= 3 && x < 3.5) syst = 0.00207952;
      else if (x >= 3.5 && x < 4) syst = 0.00964032;
      else if (x >= 4 && x < 4.5) syst = 0.00613411;
      else if (x >= 4.5 && x < 5.5) syst = 0.00141201;
      else if (x >= 5.5 && x < 6.5) syst = 0.00675142;
      else if (x >= 6.5 && x < 8) syst = 0.00255167;
      else if (x >= 8 && x < 9.5) syst = 0.00276745;
      else if (x >= 9.5 && x < 13) syst = 0.00195002;
      else syst = 0.0272819;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.5 && x < 2.2) syst = 0.00249486;
      else if (x >= 2.2 && x < 2.7) syst = 0.00185886;
      else if (x >= 2.7 && x < 3.2) syst = 0.0148211;
      else if (x >= 3.2 && x < 3.7) syst = 0.00395368;
      else if (x >= 3.7 && x < 4.7) syst = 0.00413762;
      else if (x >= 4.7 && x < 6.5) syst = 0.0162388;
      else if (x >= 6.5 && x < 8.5) syst = 0.0037397;
      else if (x >= 8.5 && x < 11) syst = 0.00190913;
      else syst = 0.0175906;
    }
  }
  else if (filterId==2) { //L2 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.687436; den = 0.628026; statUp = 0.00741958; statDown = 0.00747367;}
      else if (x >= 4 && x <4.5) {num = 0.879668; den = 0.855425; statUp = 0.00516048; statDown = 0.0052745;}
      else if (x >= 4.5 && x <5) {num = 0.936813; den = 0.91093; statUp = 0.00414716; statDown = 0.00429608;}
      else if (x >= 5 && x <5.5) {num = 0.948744; den = 0.930006; statUp = 0.00405143; statDown = 0.00424179;}
      else if (x >= 5.5 && x <6.5) {num = 0.964058; den = 0.95091; statUp = 0.00286448; statDown = 0.00299563;}
      else if (x >= 6.5 && x <8) {num = 0.967607; den = 0.967695; statUp = 0.0027257; statDown = 0.00286633;}
      else if (x >= 8 && x <10.5) {num = 0.977261; den = 0.975624; statUp = 0.0024585; statDown = 0.00262625;}
      else if (x >= 10.5 && x <14) {num = 0.975199; den = 0.979933; statUp = 0.00326164; statDown = 0.00354865;}
      else if (x >= 14 && x <18) {num = 0.978416; den = 0.983164; statUp = 0.00450526; statDown = 0.00515026;}
      else {num = 0.986566; den = 0.984365; statUp = 0.00413127; statDown = 0.00501923;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.707798; den = 0.663284; statUp = 0.0115376; statDown = 0.0116219;}
      else if (x >= 3 && x <3.5) {num = 0.841787; den = 0.79966; statUp = 0.00756481; statDown = 0.00765363;}
      else if (x >= 3.5 && x <4) {num = 0.915148; den = 0.874984; statUp = 0.0058153; statDown = 0;}
      else if (x >= 4 && x <4.5) {num = 0.932092; den = 0.925267; statUp = 0.00578591; statDown = 0.0060332;}
      else if (x >= 4.5 && x <5) {num = 0.95135; den = 0.943497; statUp = 0.005532; statDown = 0.00587852;}
      else if (x >= 5 && x <6) {num = 0.968164; den = 0.958147; statUp = 0.00392457; statDown = 0.00419381;}
      else if (x >= 6 && x <7.5) {num = 0.962352; den = 0.96425; statUp = 0.00435955; statDown = 0.00466282;}
      else if (x >= 7.5 && x <10) {num = 0.979864; den = 0.969778; statUp = 0.00356325; statDown = 0.00395684;}
      else if (x >= 10 && x <15) {num = 0.983282; den = 0.973169; statUp = 0.00384401; statDown = 0.00442388;}
      else {num = 0.985422; den = 0.979848; statUp = 0.00507078; statDown = 0.00689321;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.651185; den = 0.562127; statUp = 0.0311006; statDown = 0.0309056;}
      else if (x >= 2 && x <2.5) {num = 0.75568; den = 0.744834; statUp = 0.0139774; statDown = 0.0139217;}
      else if (x >= 2.5 && x <3) {num = 0.895072; den = 0.898086; statUp = 0.00967712; statDown = 0.00997311;}
      else if (x >= 3 && x <3.5) {num = 0.911965; den = 0.932775; statUp = 0.0093472; statDown = 0.00978825;}
      else if (x >= 3.5 && x <4) {num = 0.933577; den = 0.949087; statUp = 0.008866; statDown = 0.00937714;}
      else if (x >= 4 && x <4.5) {num = 0.92098; den = 0.952348; statUp = 0.0100317; statDown = 0.0107328;}
      else if (x >= 4.5 && x <5.5) {num = 0.945015; den = 0.959924; statUp = 0.00698531; statDown = 0.00744697;}
      else if (x >= 5.5 && x <6.5) {num = 0.915225; den = 0.961517; statUp = 0.0105111; statDown = 0.0112266;}
      else if (x >= 6.5 && x <8) {num = 0.942614; den = 0.958878; statUp = 0.00845552; statDown = 0.00920477;}
      else if (x >= 8 && x <9.5) {num = 0.943992; den = 0.95865; statUp = 0.0111118; statDown = 0.0123342;}
      else if (x >= 9.5 && x <13) {num = 0.953176; den = 0.956512; statUp = 0.0100647; statDown = 0.0114277;}
      else {num = 0.9711; den = 0.959808; statUp = 0.0112001; statDown = 0.0137596;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.80688; den = 0.804533; statUp = 0.017851; statDown = 0.0179625;}
      else if (x >= 2.2 && x <2.7) {num = 0.864384; den = 0.87583; statUp = 0.0134953; statDown = 0.0137959;}
      else if (x >= 2.7 && x <3.2) {num = 0.892179; den = 0.904052; statUp = 0.0121752; statDown = 0.0125694;}
      else if (x >= 3.2 && x <3.7) {num = 0.896585; den = 0.926574; statUp = 0.0112721; statDown = 0.0117765;}
      else if (x >= 3.7 && x <4.7) {num = 0.910998; den = 0.928332; statUp = 0.00974335; statDown = 0.010226;}
      else if (x >= 4.7 && x <6.5) {num = 0.938687; den = 0.949857; statUp = 0.00795089; statDown = 0.00845509;}
      else if (x >= 6.5 && x <8.5) {num = 0.941825; den = 0.954442; statUp = 0.0109965; statDown = 0.0143609;}
      else if (x >= 8.5 && x <11) {num = 0.931535; den = 0.959109; statUp = 0.0136929; statDown = 0.0156881;}
      else {num = 0.969056; den = 0.967084; statUp = 0.00977435; statDown = 0.0126053;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00119995;
      else if (x >= 4 && x < 4.5) syst = 0.000801484;
      else if (x >= 4.5 && x < 5) syst = 0.00142786;
      else if (x >= 5 && x < 5.5) syst = 0.000859141;
      else if (x >= 5.5 && x < 6.5) syst = 0.000855793;
      else if (x >= 6.5 && x < 8) syst = 0.000338442;
      else if (x >= 8 && x < 10.5) syst = 0.000905661;
      else if (x >= 10.5 && x < 14) syst = 0.000193737;
      else if (x >= 14 && x < 18) syst = 0.000621028;
      else syst = 0.0029276;
    }
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00301699;
      else if (x >= 3 && x < 3.5) syst = 0.0051637;
      else if (x >= 3.5 && x < 4) syst = 0.00271564;
      else if (x >= 4 && x < 4.5) syst = 0.00128082;
      else if (x >= 4.5 && x < 5) syst = 0.00105614;
      else if (x >= 5 && x < 6) syst = 0.00120191;
      else if (x >= 6 && x < 7.5) syst = 0.000729975;
      else if (x >= 7.5 && x < 10) syst = 0.00139352;
      else if (x >= 10 && x < 15) syst = 0.00151879;
      else syst = 0.00138277;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0234362;
      else if (x >= 2 && x < 2.5) syst = 0.00781699;
      else if (x >= 2.5 && x < 3) syst = 0.0020642;
      else if (x >= 3 && x < 3.5) syst = 0.00494294;
      else if (x >= 3.5 && x < 4) syst = 0.00372959;
      else if (x >= 4 && x < 4.5) syst = 0.0101533;
      else if (x >= 4.5 && x < 5.5) syst = 0.00248577;
      else if (x >= 5.5 && x < 6.5) syst = 0.00480156;
      else if (x >= 6.5 && x < 8) syst = 0.00535204;
      else if (x >= 8 && x < 9.5) syst = 0.00407749;
      else if (x >= 9.5 && x < 13) syst = 0.000734987;
      else syst = 0.0152108;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0547508;
      else if (x >= 2.2 && x < 2.7) syst = 0.0439035;
      else if (x >= 2.7 && x < 3.2) syst = 0.0100721;
      else if (x >= 3.2 && x < 3.7) syst = 0.00486924;
      else if (x >= 3.7 && x < 4.7) syst = 0.0164241;
      else if (x >= 4.7 && x < 6.5) syst = 0.0045128;
      else if (x >= 6.5 && x < 8.5) syst = 0.00615735;
      else if (x >= 8.5 && x < 11) syst = 0.00521994;
      else syst = 0.00496602;
    }
  }
  else if (filterId==3) { //L3 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0981413; den = 0.0714076; statUp = 0.00475984; statDown = 0.00464341;}
      else if (x >= 4 && x <4.5) {num = 0.309591; den = 0.234967; statUp = 0.00731017; statDown = 0.00724988;}
      else if (x >= 4.5 && x <5) {num = 0.49696; den = 0.427491; statUp = 0.00850388; statDown = 0.00850324;}
      else if (x >= 5 && x <5.5) {num = 0.646567; den = 0.569805; statUp = 0.00897182; statDown = 0.00902901;}
      else if (x >= 5.5 && x <6.5) {num = 0.717727; den = 0.665698; statUp = 0.00690496; statDown = 0.00696443;}
      else if (x >= 6.5 && x <8) {num = 0.771046; den = 0.736859; statUp = 0.00662127; statDown = 0.00670002;}
      else if (x >= 8 && x <10.5) {num = 0.792067; den = 0.777534; statUp = 0.00684886; statDown = 0.00695048;}
      else if (x >= 10.5 && x <14) {num = 0.826589; den = 0.814236; statUp = 0.00832558; statDown = 0.00852162;}
      else if (x >= 14 && x <18) {num = 0.800339; den = 0.820918; statUp = 0.0131166; statDown = 0.0135246;}
      else {num = 0.846856; den = 0.837225; statUp = 0.0139208; statDown = 0.0145458;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.307823; den = 0.284114; statUp = 0.0111767; statDown = 0.0110334;}
      else if (x >= 3 && x <3.5) {num = 0.429139; den = 0.424849; statUp = 0.00992916; statDown = 0.00988034;}
      else if (x >= 3.5 && x <4) {num = 0.54449; den = 0.527662; statUp = 0.45551; statDown = 0.00969908;}
      else if (x >= 4 && x <4.5) {num = 0.591156; den = 0.604174; statUp = 0.0113833; statDown = 0.0114203;}
      else if (x >= 4.5 && x <5) {num = 0.639967; den = 0.645913; statUp = 0.0126515; statDown = 0.0127542;}
      else if (x >= 5 && x <6) {num = 0.673449; den = 0.679156; statUp = 0.0104608; statDown = 0.010556;}
      else if (x >= 6 && x <7.5) {num = 0.685635; den = 0.720417; statUp = 0.0109999; statDown = 0.0111263;}
      else if (x >= 7.5 && x <10) {num = 0.749011; den = 0.754465; statUp = 0.0113948; statDown = 0.0115976;}
      else if (x >= 10 && x <15) {num = 0.773253; den = 0.801728; statUp = 0.0133331; statDown = 0.0136666;}
      else {num = 0.773038; den = 0.833024; statUp = 0.0211603; statDown = 0.0220621;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.00155461; den = 0.000463755; statUp = 0.00211758; statDown = 0.00108585;}
      else if (x >= 2 && x <2.5) {num = 0.00555387; den = 0.00858885; statUp = 0.00315843; statDown = 0.00555387;}
      else if (x >= 2.5 && x <3) {num = 0.451703; den = 0.447277; statUp = 0.0150343; statDown = 0.0149314;}
      else if (x >= 3 && x <3.5) {num = 0.655861; den = 0.650893; statUp = 0.0159254; statDown = 0.0160134;}
      else if (x >= 3.5 && x <4) {num = 0.706533; den = 0.710335; statUp = 0.0158931; statDown = 0.0161213;}
      else if (x >= 4 && x <4.5) {num = 0.726282; den = 0.741482; statUp = 0.0174238; statDown = 0.017772;}
      else if (x >= 4.5 && x <5.5) {num = 0.764391; den = 0.796199; statUp = 0.013155; statDown = 0.0134088;}
      else if (x >= 5.5 && x <6.5) {num = 0.769821; den = 0.824468; statUp = 0.016319; statDown = 0.0167672;}
      else if (x >= 6.5 && x <8) {num = 0.811763; den = 0.834174; statUp = 0.0145482; statDown = 0.0150501;}
      else if (x >= 8 && x <9.5) {num = 0.819571; den = 0.841319; statUp = 0.0190366; statDown = 0.01987;}
      else if (x >= 9.5 && x <13) {num = 0.829677; den = 0.857122; statUp = 0.0192692; statDown = 0.0202277;}
      else {num = 0.874981; den = 0.874474; statUp = 0.0258047; statDown = 0.0277093;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.00366548; den = 0.000413051; statUp = 0.00211949; statDown = 0.00167583;}
      else if (x >= 2.2 && x <2.7) {num = 0.116176; den = 0.109916; statUp = 0.0108872; statDown = 0.0105043;}
      else if (x >= 2.7 && x <3.2) {num = 0.413123; den = 0.401827; statUp = 0.0181733; statDown = 0.017947;}
      else if (x >= 3.2 && x <3.7) {num = 0.482035; den = 0.491334; statUp = 0.0184417; statDown = 0.0183373;}
      else if (x >= 3.7 && x <4.7) {num = 0.568894; den = 0.573223; statUp = 0.0169342; statDown = 0.0169658;}
      else if (x >= 4.7 && x <6.5) {num = 0.675048; den = 0.651616; statUp = 0.0160907; statDown = 0.0162558;}
      else if (x >= 6.5 && x <8.5) {num = 0.722882; den = 0.711847; statUp = 0.0215971; statDown = 0.022054;}
      else if (x >= 8.5 && x <11) {num = 0.714358; den = 0.750096; statUp = 0.0275205; statDown = 0.0283679;}
      else {num = 0.753355; den = 0.779864; statUp = 0.0323455; statDown = 0.0336041;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000650706;
      else if (x >= 4 && x < 4.5) syst = 0.0010869;
      else if (x >= 4.5 && x < 5) syst = 0.00298052;
      else if (x >= 5 && x < 5.5) syst = 0.00341277;
      else if (x >= 5.5 && x < 6.5) syst = 0.000613358;
      else if (x >= 6.5 && x < 8) syst = 0.000658119;
      else if (x >= 8 && x < 10.5) syst = 0.000756756;
      else if (x >= 10.5 && x < 14) syst = 0.000662617;
      else if (x >= 14 && x < 18) syst = 0.00220571;
      else syst = 0.00215326;
    }
    else if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00406324;
      else if (x >= 3 && x < 3.5) syst = 0.00422745;
      else if (x >= 3.5 && x < 4) syst = 0.00493964;
      else if (x >= 4 && x < 4.5) syst = 0.0015019;
      else if (x >= 4.5 && x < 5) syst = 0.00349953;
      else if (x >= 5 && x < 6) syst = 0.00165421;
      else if (x >= 6 && x < 7.5) syst = 0.00195686;
      else if (x >= 7.5 && x < 10) syst = 0.00305233;
      else if (x >= 10 && x < 15) syst = 0.00341103;
      else syst = 0.00425449;
    }
    else if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.00154707;
      else if (x >= 2 && x < 2.5) syst = 0.00239348;
      else if (x >= 2.5 && x < 3) syst = 0.00578521;
      else if (x >= 3 && x < 3.5) syst = 0.0019864;
      else if (x >= 3.5 && x < 4) syst = 0.00826595;
      else if (x >= 4 && x < 4.5) syst = 0.00622458;
      else if (x >= 4.5 && x < 5.5) syst = 0.00155048;
      else if (x >= 5.5 && x < 6.5) syst = 0.00738518;
      else if (x >= 6.5 && x < 8) syst = 0.00155169;
      else if (x >= 8 && x < 9.5) syst = 0.00373986;
      else if (x >= 9.5 && x < 13) syst = 0.00445251;
      else syst = 0.028681;
    }
    else if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.000519911;
      else if (x >= 2.2 && x < 2.7) syst = 0.0288676;
      else if (x >= 2.7 && x < 3.2) syst = 0.013137;
      else if (x >= 3.2 && x < 3.7) syst = 0.00153582;
      else if (x >= 3.7 && x < 4.7) syst = 0.00361851;
      else if (x >= 4.7 && x < 6.5) syst = 0.0155374;
      else if (x >= 6.5 && x < 8.5) syst = 0.00321391;
      else if (x >= 8.5 && x < 11) syst = 0.00306389;
      else syst = 0.0199929;
    }
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  else if (idx == -2) syst_factor = -1*syst;
  else if (idx == +1) stat_factor = statUp;
  else if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//              T R K     P b P b                //
///////////////////////////////////////////////////

double tnp_weight_trk_pbpb(double eta, int idx) {
  double x = eta;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  //SF in eta bins
  if (x >= -2.4 && x < -1.6) {num = 0.994498; den = 0.998413; statUp = 0.00287386; statDown = 0.00290888;}
  if (x >= -1.6 && x < -1.2) {num = 0.973539; den = 0.967322; statUp = 0.00399501; statDown = 0.00407843;}
  if (x >= -1.2 && x < -0.9) {num = 0.964465; den = 0.970816; statUp = 0.00861188; statDown = 0.00869453;}
  if (x >= -0.9 && x < -0.6) {num = 0.96081; den = 0.974407; statUp = 0.0223599; statDown = 0.00682405;}
  if (x >= -0.6 && x < -0.3) {num = 0.964464; den = 0.97802; statUp = 0.00612474; statDown = 0.00621291;}
  if (x >= -0.3 && x < 0.3) {num = 0.963862; den = 0.966583; statUp = 0.00496914; statDown = 0.00503378;}
  if (x >= 0.3 && x < 0.6) {num = 0.956897; den = 0.967248; statUp = 0.00672757; statDown = 0.0068202;}
  if (x >= 0.6 && x < 0.9) {num = 0.964172; den = 0.966882; statUp = 0.00735892; statDown = 0.00754429;}
  if (x >= 0.9 && x < 1.2) {num = 0.961874; den = 0.955987; statUp = 0.00987473; statDown = 0.0099638;}
  if (x >= 1.2 && x < 1.6) {num = 0.964754; den = 0.964653; statUp = 0.0042287; statDown = 0.00430601;}
  if (x >= 1.6 && x < 2.4) {num = 0.999937; den = 0.998771; statUp = 6.33084e-05; statDown = 0.00310832;}

  // syst uncertainties
  if (x >= -2.4 && x < -1.6) syst = 0.0015001;
  else if (x >= -1.6 && x < -1.2) syst = 0.00376932;
  else if (x >= -1.2 && x < -0.9) syst = 0.00125496;
  else if (x >= -0.9 && x < -0.6) syst = 0.00190534;
  else if (x >= -0.6 && x < -0.3) syst = 0.00228604;
  else if (x >= -0.3 && x < 0.3) syst = 0.00493996;
  else if (x >= 0.3 && x < 0.6) syst = 0.00527961;
  else if (x >= 0.6 && x < 0.9) syst = 0.00231575;
  else if (x >= 0.9 && x < 1.2) syst = 0.012059;
  else if (x >= 1.2 && x < 1.6) syst = 0.00278996;
  else syst = 0.000876099;

  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  else if (idx == -2) syst_factor = -1*syst;
  else if (idx == +1) stat_factor = statUp;
  else if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

#endif
