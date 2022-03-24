// -*- C++ -*-
// author: afiq anuar
// short: functions pertaining to converting histogram representation used in this project back to ROOT and saving into output

#ifndef FWK_OUTPUT_UTIL_H
#define FWK_OUTPUT_UTIL_H

#include "string_io.h"
#include "array_histogram.h"

#include "TFile.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

/// for restype, meaningful chars are u. n and a for unroll, normal and all
/// see --result-type option in smoother.cc for explanation
std::vector<std::unique_ptr<TH1>> array_to_root(const std::tuple<
                                                std::vector<std::vector<std::string>>,
                                                std::vector<std::vector<double>>,
                                                std::vector<std::vector<double>>,
                                                std::string> &variables,
                                                const std::string &tag,
                                                const std::vector<std::vector<double>> &edges,
                                                const Arrayhist &hist,
                                                char restype = 'u')
{
  std::vector<std::unique_ptr<TH1>> result;
  if (edges.empty())
    return result;

  const int nbin = count_nbin(edges), nvar = edges.size();
  if (hist.size() != nbin)
    return result;

  const std::string vars = [nvar] (const auto& variables) {
    auto str = std::get<0>(variables)[0][0];
    for (int istr = 1; istr < nvar; ++istr)
      str = str + "_" + std::get<0>(variables)[istr][0];

    return str;
  }(variables);
  const std::string ftag = (tag != "") ? "_" + tag : tag;

  std::unique_ptr<TH1> unroll = nullptr, normal = nullptr;

  if (nvar == 1) {
    unroll = std::make_unique<TH1D>((vars + ftag).c_str(), join({"", std::get<0>(variables)[0][0]}, ";").c_str(),
                                    nbin, edges[0].data());
  }
  else
    unroll = std::make_unique<TH1D>((vars + "_unroll" + ftag).c_str(), ";unrolled bin index", nbin, 0., nbin);

  if (restype != 'u' and nvar == 2) {
    normal = std::make_unique<TH2D>((vars + ftag).c_str(), 
                                    join({"", std::get<0>(variables)[0][0], std::get<0>(variables)[1][0]}, ";").c_str(),
                                    edges[0].size() - 1, edges[0].data(), edges[1].size() - 1, edges[1].data());
  }
  if (restype != 'u' and nvar == 3) {
    normal = std::make_unique<TH3D>((vars + ftag).c_str(),
                                    join({"",
                                          std::get<0>(variables)[0][0],
                                          std::get<0>(variables)[1][0],
                                          std::get<0>(variables)[2][0]}, ";").c_str(),
                                    edges[0].size() - 1, edges[0].data(),
                                    edges[1].size() - 1, edges[1].data(),
                                    edges[2].size() - 1, edges[2].data());
  }

  static std::vector<int> idx(nvar, -1);
  if (nvar > idx.capacity())
    idx.reserve(nvar);

  for (int ibin = 0; ibin < nbin; ++ibin) {
    if (std::isnan(hist(ibin, 0)) or std::isnan(hist(ibin, 1)) or hist(ibin, 1) < 0.)
      continue;

    unroll->SetBinContent(ibin + 1, hist(ibin, 0));
    unroll->SetBinError(ibin + 1, std::sqrt(hist(ibin, 1)));

    if (restype != 'u' and (nvar == 2 or nvar == 3)) {
      for (int iv = 0; iv < nvar; ++iv)
        idx[iv] = index_1n(ibin, iv, edges);

      if (nvar == 2) {
        normal->SetBinContent(idx[0] + 1, idx[1] + 1, hist(ibin, 0));
        normal->SetBinError(idx[0] + 1, idx[1] + 1, std::sqrt(hist(ibin, 1)));
      }
      else if (nvar == 3) {
        normal->SetBinContent(idx[0] + 1, idx[1] + 1, idx[2] + 1, hist(ibin, 0));
        normal->SetBinError(idx[0] + 1, idx[1] + 1, idx[2] + 1, std::sqrt(hist(ibin, 1)));
      }
    }
  }

  if (restype == 'u' or restype == 'a')
    result.emplace_back( std::move(unroll) );
  if (normal != nullptr and restype != 'u')
    result.emplace_back( std::move(normal) );

  return result;
}



template <typename ...Hists>
void save_all_as(const std::string &name, const Hists &...hists)
{
  std::array<std::reference_wrapper<const std::vector<std::unique_ptr<TH1>>>, sizeof...(hists)> refs = { std::cref(hists)... };

  // create smallest possible file
  auto file = std::make_unique<TFile>(name.c_str(), "recreate", "", 209);
  file->cd();

  for (const auto &ref : refs) {
    for (const auto &h : ref.get())
      h->Write();
  }
}

#endif
