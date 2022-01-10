// -*- C++ -*-
// author: afiq anuar
// short: for computing ttbar spin correlation variables

#ifndef FWK_W_HELICITY_H
#define FWK_W_HELICITY_H

#include "misc/container_util.h"
#include "misc/constants.h"
#include "misc/numeric_vector.h"
#include <cmath>

template <typename Number = float>
const std::vector<std::pair<std::string, Number>>&
compute_w_helicity(Number pTop_pt, Number pTop_eta, Number pTop_phi, Number pTop_m,
                         Number aTop_pt, Number aTop_eta, Number aTop_phi, Number aTop_m,
                         Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                         Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m,
                         Number pW_pt, Number pW_eta, Number pW_phi, Number pW_m,
                         Number aW_pt, Number aW_eta, Number aW_phi, Number aW_m,
                         Number pB_pt, Number pB_eta, Number pB_phi, Number pB_m,
                         Number aB_pt, Number aB_eta, Number aB_phi, Number aB_m)

{
  using namespace Framework;
  static std::vector<std::pair<std::string, Number>> m_w_hel;
  static int initialize = 0;
  if (initialize == 0) {
    m_w_hel.reserve(4); // exact size


    // m_w_hel.emplace_back("f1A", -9999.);
    // m_w_hel.emplace_back("f2A", -9999.);
    m_w_hel.emplace_back("f1", -9999.);
    m_w_hel.emplace_back("f2", -9999.);
    // m_w_hel.emplace_back("f1paper", -9999.);
    // m_w_hel.emplace_back("f2paper", -9999.);

    ++initialize;
  }

  static std::array<Number, 32> arg;
  if (arg[0]  == pTop_pt and arg[1]  == pTop_eta and arg[2]  == pTop_phi and arg[3]  == pTop_m and
      arg[4]  == aTop_pt and arg[5]  == aTop_eta and arg[6]  == aTop_phi and arg[7]  == aTop_m and
      arg[8]  == pLep_pt and arg[9]  == pLep_eta and arg[10] == pLep_phi and arg[11] == pLep_m and
      arg[12] == aLep_pt and arg[13] == aLep_eta and arg[14] == aLep_phi and arg[15] == aLep_m and
      arg[16] == pW_pt and arg[17] == pW_eta and arg[18] == pW_phi and arg[19] == pW_m and
      arg[20] == aW_pt and arg[21] == aW_eta and arg[22] == aW_phi and arg[23] == aW_m and
      arg[24] == pB_pt and arg[25] == pB_eta and arg[26] == pB_phi and arg[27] == pB_m and
      arg[28] == aB_pt and arg[29] == aB_eta and arg[30] == aB_phi and arg[31] == aB_m)
    return m_w_hel;

  arg[0]  = pTop_pt; arg[1]  = pTop_eta; arg[2]  = pTop_phi; arg[3]  = pTop_m;
  arg[4]  = aTop_pt; arg[5]  = aTop_eta; arg[6]  = aTop_phi; arg[7]  = aTop_m;
  arg[8]  = pLep_pt; arg[9]  = pLep_eta; arg[10] = pLep_phi; arg[11] = pLep_m;
  arg[12] = aLep_pt; arg[13] = aLep_eta; arg[14] = aLep_phi; arg[15] = aLep_m;
  arg[16] = pW_pt; arg[17] = pW_eta; arg[18] = pW_phi; arg[19] = pW_m;
  arg[20] = aW_pt; arg[21] = aW_eta; arg[22] = aW_phi; arg[23] = aW_m;
  arg[24] = pB_pt; arg[25] = pB_eta; arg[26] = pB_phi; arg[27] = pB_m;
  arg[28] = aB_pt; arg[29] = aB_eta; arg[30] = aB_phi; arg[31] = aB_m;


  static TLorentzVector p4lab_pTop, p4lab_aTop, p4lab_pLep, p4lab_aLep, p4lab_pW, p4lab_aW,p4lab_pB, p4lab_aB;
  p4lab_pTop.SetPtEtaPhiM(arg[0],  arg[1],  arg[2],  arg[3]);
  p4lab_aTop.SetPtEtaPhiM(arg[4],  arg[5],  arg[6],  arg[7]);
  p4lab_pLep.SetPtEtaPhiM(arg[8],  arg[9],  arg[10], arg[11]);
  p4lab_aLep.SetPtEtaPhiM(arg[12], arg[13], arg[14], arg[15]);
  p4lab_pW.SetPtEtaPhiM(arg[16], arg[17], arg[18], arg[19]);
  p4lab_aW.SetPtEtaPhiM(arg[20], arg[21], arg[22], arg[23]);
  p4lab_pB.SetPtEtaPhiM(arg[24], arg[25], arg[26], arg[27]);
  p4lab_aB.SetPtEtaPhiM(arg[28], arg[29], arg[30], arg[31]);

  // the necessary boosting such that lepton ~ top spin vector
  const TLorentzVector p4lab_TT( p4lab_pTop + p4lab_aTop ); //calculates the t-tbar system 4-momenta
  // const TLorentzVector p4lab_eb( p4lab_aLep + p4lab_pB ); //calculates the t-tbar system 4-momenta


  const auto f_zmf_tt = [&p4lab_TT] (TLorentzVector p4) {
    // function that boost a 4-vector p4 to the zmf-tt sytem
    // p4lab_TT.BoostVector() returns the beta vector (relative speed) of p4lab_TT (TTbar system) with respect to the lab system
    // p4.Boost(-1*vec{beta}) peforms the direct lorentz transform from the lab system to the TTbar system, where beta is the relative speed between them
    p4.Boost( -1. * p4lab_TT.BoostVector() ); //Boost() is the inverse Lorentz transformation, therefore the -1 factor to the relative velocity (beta vector) to give the direct Lorentz transform.
    return p4;
  };
  // const TLorentzVector p4hel_pTop = f_zmf_tt(p4lab_pTop);
  // const TLorentzVector p4hel_aTop = f_zmf_tt(p4lab_aTop);

  // const auto f_zmf_top = [&f_zmf_tt] (const TLorentzVector &wboson, const TLorentzVector &top) {
  //   // function that boost a 4-vector W to the direction of the 4-momentum of the TT sytem, and then reboost it to the direction of the 4-vector top
  //   auto p4 = f_zmf_tt(wboson); // boost a 4-vector wboson to the direction of the TT zmf sytem
  //   p4.Boost( -1. * top.BoostVector() ); // inverse of boosting the 4-vec top to the system where p4 is defined, i.e the tt-zmf. I.e. boost p4 from tt_zmf to the top zmf
  //   return p4;
  // };


  const auto f_zmf_w = [&f_zmf_tt] (const TLorentzVector &daugther, const TLorentzVector &wboson) {
    // function that boost a 4-vector u4 to the direction of the 4-momentum of the TT sytem, and then reboost it to the direction of the 4-vector w
    auto p4 = f_zmf_tt(daugther);
    p4.Boost( -1. * wboson.BoostVector() );
    return p4;
  };

  // const TLorentzVector p4top_pW = f_zmf_top(p4lab_pW, p4hel_pTop); //directon of flight of w in the t-ZMF
  // const TLorentzVector p4top_aW = f_zmf_top(p4lab_aW, p4hel_aTop);

  const TLorentzVector p4hel_pW = f_zmf_tt(p4lab_pW); //directon of flight of W+ in the tt-ZMF
  const TLorentzVector p4hel_aW = f_zmf_tt(p4lab_aW);

  const TLorentzVector p4hel_aLep = f_zmf_w(p4lab_aLep, p4hel_pW); //directon of flight of anti-lepton in the w-ZMF
  const TLorentzVector p4hel_pLep = f_zmf_w(p4lab_pLep, p4hel_aW);

  const TLorentzVector p4hel_pB = f_zmf_w(p4lab_pB, p4hel_pW); //directon of flight of b-quark in the w-ZMF
  const TLorentzVector p4hel_aB = f_zmf_w(p4lab_aB, p4hel_aW);

  // const TLorentzVector p4top_pB = f_zmf_top(p4lab_pB, p4hel_pTop); //directon of flight of b-quark in the t-ZMF
  // const TLorentzVector p4top_aB = f_zmf_top(p4lab_aB, p4hel_aTop);
  //

  // opening angle between the lep direction in w frame and w direction in t frame
  // static const int if1A = index_with_key(m_w_hel, "f1");
  // m_w_hel[if1A].second = p4hel_aLep.Vect().Unit().Dot( p4top_pW.Vect().Unit() );
  //
  // static const int if2A = index_with_key(m_w_hel, "f2A");
  // m_w_hel[if2A].second = p4hel_pLep.Vect().Unit().Dot( p4top_aW.Vect().Unit() );

  // opening angle between the lep direction in w frame and the reversed b direction in w frame
  static const int if1 = index_with_key(m_w_hel, "f1");
  m_w_hel[if1].second = p4hel_aLep.Vect().Unit().Dot( -p4hel_pB.Vect().Unit() );

  static const int if2 = index_with_key(m_w_hel, "f2");
  m_w_hel[if2].second = p4hel_pLep.Vect().Unit().Dot( -p4hel_aB.Vect().Unit() );

  // static const int if1paper = index_with_key(m_w_hel, "f1paper");
  // m_w_hel[if1paper].second = 2*p4lab_eb.M()*p4lab_eb.M()  / (p4lab_pTop.M()*p4lab_pTop.M() - p4lab_pW.M()*p4lab_pW.M()) - 1;


  // static const int if1paper = index_with_key(m_w_hel, "f1paper");
  // m_w_hel[if1paper].second = ( p4hel_aLep*p4hel_pB -  p4hel_aLep.E()*p4hel_pB.E() )  / ( sqrt(p4hel_aLep.Vect().Dot(p4hel_aLep.Vect())) * sqrt(p4hel_pB.Vect().Dot( p4hel_pB.Vect() ) ));
  //
  // // static const int if1paper = index_with_key(m_w_hel, "f1paper");
  // // m_w_hel[if1paper].second = p4hel_aLep.Vect().Unit().Dot( -p4top_pB.Vect().Unit() );
  //
  // static const int if2paper = index_with_key(m_w_hel, "f2paper");
  // m_w_hel[if2paper].second = p4hel_pLep.Vect().Unit().Dot( -p4top_aB.Vect().Unit() );


  return m_w_hel;
}



template <typename Number = float>
auto w_helicity(const std::string &var)
{
  const auto &map = compute_w_helicity(Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number(),
                                             Number(), Number(), Number(), Number());

  return [ivar = index_with_key(map, var)] (Number pTop_pt, Number pTop_eta, Number pTop_phi, Number pTop_m,
                                            Number aTop_pt, Number aTop_eta, Number aTop_phi, Number aTop_m,
                                            Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                                            Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m,
                                            Number pW_pt, Number pW_eta, Number pW_phi, Number pW_m,
                                            Number aW_pt, Number aW_eta, Number aW_phi, Number aW_m,
                                            Number pB_pt, Number pB_eta, Number pB_phi, Number pB_m,
                                            Number aB_pt, Number aB_eta, Number aB_phi, Number aB_m)
  {
    const auto &map = compute_w_helicity(pTop_pt, pTop_eta, pTop_phi, pTop_m,
                                               aTop_pt, aTop_eta, aTop_phi, aTop_m,
                                               pLep_pt, pLep_eta, pLep_phi, pLep_m,
                                               aLep_pt, aLep_eta, aLep_phi, aLep_m,
                                               pW_pt, pW_eta, pW_phi, pW_m,
                                               aW_pt, aW_eta, aW_phi, aW_m,
                                               pB_pt, pB_eta, pB_phi, pB_m,
                                               aB_pt, aB_eta, aB_phi, aB_m);

    return map[ivar].second;
  };
}

#endif
