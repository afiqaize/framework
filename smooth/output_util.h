// -*- C++ -*-
// author: afiq anuar
// short: functions pertaining to converting histogram representation used in this project back to ROOT and saving into output

#ifndef FWK_OUTPUT_UTIL_H
#define FWK_OUTPUT_UTIL_H

#include "misc/string_io.h"
#include "array_histogram.h"

#include "TFile.h"

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

std::vector<std::unique_ptr<TH1>> array_to_root(const std::tuple<
                                                std::vector<std::vector<std::string>>,
                                                std::vector<std::vector<double>>,
                                                std::vector<std::vector<double>>,
                                                std::string> &variables,
                                                const std::string &tag,
                                                const std::vector<std::vector<double>> &edges,
                                                const Arrayhist &hist)
{
  std::vector<std::unique_ptr<TH1>> result;
  if (edges.empty())
    return result;

  const int nbin = count_nbin(edges);
  if (hist.size() != nbin)
    return result;

  const std::string vars = [nvar = edges.size()] (const auto& variables) {
    auto str = std::get<0>(variables)[0][0];
    for (int istr = 1; istr < nvar; ++istr)
      str = str + "_" + std::get<0>(variables)[istr][0];

    return str;
  }(variables);
  const std::string ftag = (tag != "") ? "_" + tag : tag;

  if (edges.size() == 1) {
    result.emplace_back( std::make_unique<TH1D>((vars + ftag).c_str(), join({"", std::get<0>(variables)[0][0]}, ";").c_str(),
                                                nbin, edges[0].data()) );
  }
  else
    result.emplace_back( std::make_unique<TH1D>((vars + "_unroll" + ftag).c_str(), ";unrolled bin index", nbin, 0., nbin) );

  if (edges.size() == 2) {
    result.emplace_back( std::make_unique<TH2D>((vars + ftag).c_str(), 
                                                join({"", std::get<0>(variables)[0][0], std::get<0>(variables)[1][0]}, ";").c_str(),
                                                edges[0].size() - 1, edges[0].data(), edges[1].size() - 1, edges[1].data()) );
  }
  if (edges.size() == 3) {
    result.emplace_back( std::make_unique<TH3D>((vars + ftag).c_str(),
                                                join({"",
                                                      std::get<0>(variables)[0][0],
                                                      std::get<0>(variables)[1][0],
                                                      std::get<0>(variables)[2][0]}, ";").c_str(),
                                                edges[0].size() - 1, edges[0].data(),
                                                edges[1].size() - 1, edges[1].data(),
                                                edges[2].size() - 1, edges[2].data()) );
  }

  static std::vector<int> idx(edges.size(), -1);
  if (edges.size() > idx.capacity())
    idx.reserve(edges.size());

  for (int ibin = 0; ibin < nbin; ++ibin) {
    if (std::isnan(hist(ibin, 0)) or std::isnan(hist(ibin, 1)) or hist(ibin, 1) < 0.)
      continue;

    result[0]->SetBinContent(ibin + 1, hist(ibin, 0));
    result[0]->SetBinError(ibin + 1, std::sqrt(hist(ibin, 1)));

    if (edges.size() == 2 or edges.size() == 3) {
      for (int iv = 0; iv < edges.size(); ++iv)
        idx[iv] = index_1n(ibin, iv, edges);

      if (edges.size() == 2) {
        result[1]->SetBinContent(idx[0] + 1, idx[1] + 1, hist(ibin, 0));
        result[1]->SetBinError(idx[0] + 1, idx[1] + 1, std::sqrt(hist(ibin, 1)));
      }
      else if (edges.size() == 3) {
        result[1]->SetBinContent(idx[0] + 1, idx[1] + 1, idx[2] + 1, hist(ibin, 0));
        result[1]->SetBinError(idx[0] + 1, idx[1] + 1, idx[2] + 1, std::sqrt(hist(ibin, 1)));
      }
    }
  }

  return result;
}



template <typename ...Hists>
void save_all_as(const std::string &name, const Hists &...hists)
{
  std::array<std::reference_wrapper<const std::vector<std::unique_ptr<TH1>>>, sizeof...(hists)> refs = { std::cref(hists)... };

  auto file = std::make_unique<TFile>(name.c_str(), "recreate");
  file->cd();

  for (const auto &ref : refs) {
    for (const auto &h : ref.get())
      h->Write();
  }
}

#endif
