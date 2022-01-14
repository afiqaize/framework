// -*- C++ -*-
// author: Andre Zimermmane Santos
// short: for computing W-helicity angle in the ttbar system

#ifndef FWK_W_HELICITY_H

#include "misc/numeric_vector.h"
#include "misc/container_util.h"
#include <cmath>


template <typename Number = float>
const std::vector<std::pair<std::string, Number>>&

compute_w_helicity(Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                         Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m,
                         Number pW_pt, Number pW_eta, Number pW_phi, Number pW_m,
                         Number aW_pt, Number aW_eta, Number aW_phi, Number aW_m,
                         Number pB_pt, Number pB_eta, Number pB_phi, Number pB_m,
                         Number aB_pt, Number aB_eta, Number aB_phi, Number aB_m)

{
  using namespace Framework;
  thread_local std::vector<std::pair<std::string, Number>> m_w_hel;
  thread_local int initialize = 0;
  if (initialize == 0) {
    m_w_hel.reserve(2); // exact size

    m_w_hel.emplace_back("f1", -9999.);
    m_w_hel.emplace_back("f2", -9999.);

    ++initialize;
  }

  thread_local std::array<Number, 24> arg;
  if (arg[0]  == pLep_pt and arg[1]  == pLep_eta and arg[2] == pLep_phi and arg[3] == pLep_m and
      arg[4] == aLep_pt and arg[5] == aLep_eta and arg[6] == aLep_phi and arg[7] == aLep_m and
      arg[8] == pW_pt and arg[9] == pW_eta and arg[10] == pW_phi and arg[11] == pW_m and
      arg[12] == aW_pt and arg[13] == aW_eta and arg[14] == aW_phi and arg[15] == aW_m and
      arg[16] == pB_pt and arg[17] == pB_eta and arg[18] == pB_phi and arg[19] == pB_m and
      arg[20] == aB_pt and arg[21] == aB_eta and arg[22] == aB_phi and arg[23] == aB_m)
    return m_w_hel;

  arg[0]  = pLep_pt; arg[1]  = pLep_eta; arg[2] = pLep_phi; arg[3] = pLep_m;
  arg[4] = aLep_pt; arg[5] = aLep_eta; arg[6] = aLep_phi; arg[7] = aLep_m;
  arg[8] = pW_pt; arg[9] = pW_eta; arg[10] = pW_phi; arg[11] = pW_m;
  arg[12] = aW_pt; arg[13] = aW_eta; arg[14] = aW_phi; arg[15] = aW_m;
  arg[16] = pB_pt; arg[17] = pB_eta; arg[18] = pB_phi; arg[19] = pB_m;
  arg[20] = aB_pt; arg[21] = aB_eta; arg[22] = aB_phi; arg[23] = aB_m;


  thread_local TLorentzVector p4lab_pLep, p4lab_aLep, p4lab_pW, p4lab_aW,p4lab_pB, p4lab_aB;
  p4lab_pLep.SetPtEtaPhiM(arg[0],  arg[1],  arg[2], arg[3]);
  p4lab_aLep.SetPtEtaPhiM(arg[4], arg[5], arg[6], arg[7]);
  p4lab_pW.SetPtEtaPhiM(arg[8], arg[9], arg[10], arg[11]);
  p4lab_aW.SetPtEtaPhiM(arg[12], arg[13], arg[14], arg[15]);
  p4lab_pB.SetPtEtaPhiM(arg[16], arg[17], arg[18], arg[19]);
  p4lab_aB.SetPtEtaPhiM(arg[20], arg[21], arg[22], arg[23]);

  const auto f_zmf_w = [] (const TLorentzVector &daugther, const TLorentzVector &wboson) {
    // function that boost a 4-vector u4 to the direction of the 4-momentum of the TT sytem, and then reboost it to the direction of the 4-vector w
    auto p4 = daugther;
    auto u4 = wboson;
    p4.Boost( -1. * u4.BoostVector());
    return p4;
  };


  const TLorentzVector p4hel_aLep = f_zmf_w(p4lab_aLep, p4lab_pW); //directon of flight of anti-lepton in the w-ZMF
  const TLorentzVector p4hel_pLep = f_zmf_w(p4lab_pLep, p4lab_aW);

  const TLorentzVector p4hel_pB = f_zmf_w(p4lab_pB, p4lab_pW); //directon of flight of b-quark in the w-ZMF
  const TLorentzVector p4hel_aB = f_zmf_w(p4lab_aB, p4lab_aW);


  // opening angle between the lep direction in w frame and the reversed b direction in w frame
  static const int if1 = index_with_key(m_w_hel, "f1");
  m_w_hel[if1].second = p4hel_aLep.Vect().Unit().Dot( -p4hel_pB.Vect().Unit() );

  static const int if2 = index_with_key(m_w_hel, "f2");
  m_w_hel[if2].second = p4hel_pLep.Vect().Unit().Dot( -p4hel_aB.Vect().Unit() );


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
                                             Number(), Number(), Number(), Number());

  return [ivar = index_with_key(map, var)] (Number pLep_pt, Number pLep_eta, Number pLep_phi, Number pLep_m,
                                            Number aLep_pt, Number aLep_eta, Number aLep_phi, Number aLep_m,
                                            Number pW_pt, Number pW_eta, Number pW_phi, Number pW_m,
                                            Number aW_pt, Number aW_eta, Number aW_phi, Number aW_m,
                                            Number pB_pt, Number pB_eta, Number pB_phi, Number pB_m,
                                            Number aB_pt, Number aB_eta, Number aB_phi, Number aB_m)
  {
    const auto &map = compute_w_helicity(pLep_pt, pLep_eta, pLep_phi, pLep_m,
                                               aLep_pt, aLep_eta, aLep_phi, aLep_m,
                                               pW_pt, pW_eta, pW_phi, pW_m,
                                               aW_pt, aW_eta, aW_phi, aW_m,
                                               pB_pt, pB_eta, pB_phi, pB_m,
                                               aB_pt, aB_eta, aB_phi, aB_m);

    return map[ivar].second;
  };
}

#endif
