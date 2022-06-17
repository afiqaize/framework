// -*- C++ -*-
// author: afiq anuar
// short: an interface for constructing a lookup table and computing weights due to prescales
// note: the weights require as input four numbers; run, lumi, and two bitsets representing the hlt and l1 bit masks
// note: the code takes as input the json file generated by bparking/bril_prescale.py

#ifndef FWK_PRESCALER_H
#define FWK_PRESCALER_H

#include <bitset>

#include "misc/container_util.h"
#include "misc/string_io.h"
#include "json/json.hpp"

template <size_t NPATH = 128, size_t NSEED = NPATH>
class Prescaler {
public:
  struct Data {
    std::vector<int> path_indices;
    std::vector<int> path_prescales;
    std::vector<std::vector<int>> seed_indices;
    std::vector<std::vector<int>> seed_prescales;
  };

  /// constructor
  Prescaler() = delete;
  Prescaler(const std::string &input, const std::vector<std::string> &keep = {}) : run(0), irun(0), lumi(0), ilumi(0) { initialize(input, keep); check(); }
  Prescaler(const std::vector<std::string> &input, const std::vector<std::string> &keep = {}) : run(0), irun(0), lumi(0), ilumi(0) {
    for (const auto &ii : input)
      initialize(ii, keep);
    check();
  }

  /// find the weight of the event; weight = 1 / effective prescale
  /// effective prescale is the probability of the event to be triggered; 1 - prod_i(1 - 1/N_i)
  /// with events passing a trigger with prescale N = HL (with H/L = hlt/l1 prescale) have a probability of 1/N to be triggered
  /// prescales are statistically independent from each other
  double weight(int run_, int lumi_, const std::bitset<NPATH> &hltbits, const std::bitset<NSEED> &l1bits);

  /// get the hlt paths stored
  const std::vector<std::string>& hlt_paths() const { return paths; };

  /// get the l1 seeds stored
  const std::vector<std::string>& l1_seeds() const { return seeds; };

private:
  /// make the table using the input file
  /// input is the input json containing the info, and keep is the list of paths one wanna keep
  /// if keep isn't empty, then everything else is ignored
  void initialize(const std::string &input, const std::vector<std::string> &keep);

  /// find the index of a given key = path or seed
  /// if not already there, adds the key to the container
  int index(const std::string &key, std::vector<std::string> &container);

  /// something's wrong if this function is called
  /// don't crash the code though
  double scream(int run_, int lumi_);

  /// validate that expected bitset is enough, and print some information
  void check() const;

  /// the lookup table storing the info
  std::vector<std::pair<int, std::vector<std::pair<int, Data>>>> prescales;

  /// list of paths and seeds
  std::vector<std::string> paths;
  std::vector<std::string> seeds;

  /// tagging members so that the query happens every run/lumi change instead of every event
  int run;
  int irun;

  int lumi;
  int ilumi;
};



template <size_t NPATH, size_t NSEED>
double Prescaler<NPATH, NSEED>::weight(int run_, int lumi_, const std::bitset<NPATH> &hltbits, const std::bitset<NSEED> &l1bits)
{
  if (run_ != run) {
    run = run_;
    irun = index_with_key(prescales, run);

    if (irun < 0)
      return scream(run_, lumi_);

    lumi = lumi_;
    ilumi = index_greater_equal(prescales[irun].second, lumi);

    if (ilumi < 0)
      return scream(run_, lumi_);
  }

  if (lumi_ != lumi) {
    lumi = lumi_;
    ilumi = index_greater_equal(prescales[irun].second, lumi);

    if (ilumi < 0)
      return scream(run_, lumi_);
  }

  const auto &data = prescales[irun].second[ilumi].second;
  std::vector<int> eps = {}; // effective prescales
  for (int ipa = 0; ipa < data.path_indices.size(); ++ipa) {
    if (not hltbits[data.path_indices[ipa]])
      continue;

    const int pps = data.path_prescales[ipa];

    std::vector<int> sps = {}; // seed prescales
    for (int ise = 0; ise < data.seed_indices[ipa].size(); ++ise) {
      if (not l1bits[data.seed_indices[ipa][ise]])
        continue;

      if (pps * data.seed_prescales[ipa][ise] == 1)
        return 1.;
      else
        sps.emplace_back(data.seed_prescales[ipa][ise]);
    }

    if (std::count(std::begin(sps), std::end(sps), 1) > 0)
      eps.emplace_back(pps);
    else {
      for (const auto &ps : sps)
        eps.emplace_back(pps * ps);
    }
  }

  return 1. / (1. - std::accumulate(std::begin(eps), std::end(eps), 1., [] (double probability, int prescale) { return probability * (1. - (1. / prescale)); }));
}



template <size_t NPATH, size_t NSEED>
void Prescaler<NPATH, NSEED>::initialize(const std::string &input, const std::vector<std::string> &keep)
{
  using json = nlohmann::json;

  std::ifstream ifile(input);
  if (not ifile)
    throw std::runtime_error( "ERROR: Prescaler fails to open the json file!!" );

  json ii;
  ifile >> ii;

  auto sort_by_key = [] (auto &p1, auto &p2) { return p1.first < p2.first; };
  auto key_is_irun = [this] (const auto &p) { return p.first == irun; };
  auto key_is_ilumi = [this] (const auto &p) { return p.first == ilumi; };

  // get the paths and seeds
  for (const auto &[srun, lumis] : ii.items()) {
    const json ll = lumis;

    // check if the run is already there
    // if yes, extend it, otherwise, make a new entry
    irun = std::stoi(srun);
    auto rps = std::find_if(std::begin(prescales), std::end(prescales), key_is_irun);

    if (rps == std::end(prescales)) {
      prescales.emplace_back(irun, std::vector<std::pair<int, Data>>{});
      rps = std::find_if(std::begin(prescales), std::end(prescales), key_is_irun);
    }

    for (const auto &[slumi, triggers] : ll.items()) {
      const json tt = triggers;

      // check if the lumi is already there
      // if yes, extend it, otherwise, make a new entry
      ilumi = std::stoi(slumi);
      auto lps = std::find_if(std::begin(rps->second), std::end(rps->second), key_is_ilumi);

      if (lps == std::end(rps->second)) {
        rps->second.emplace_back( std::make_pair(ilumi, Prescaler::Data{}) );
        lps = std::find_if(std::begin(rps->second), std::end(rps->second), key_is_ilumi);
      }
      auto &data = lps->second;

      for (const auto &[path, prescale] : tt.items()) {
        bool tokeep = keep.empty();
        for (const auto &k : keep)
          tokeep = tokeep or contain(path, k);

        if (not tokeep)
          continue;

        const int ipa = index(path, paths);
        if (std::count(std::begin(data.path_indices), std::end(data.path_indices), ipa) == 0) {
          data.path_indices.emplace_back(ipa);
          data.path_prescales.emplace_back(prescale["hlt_prescale"]);

          data.seed_indices.emplace_back(std::vector<int>{});
          for (const auto &seed : prescale["seeds"])
            data.seed_indices.back().emplace_back(index(seed, seeds));

          data.seed_prescales.emplace_back(std::vector<int>{});
          for (const auto &seed : prescale["seed_prescales"])
            data.seed_prescales.back().emplace_back(seed);
        }
      }
    }

    // lumi and run needs to be sorted, as they are arranged lexically
    // likely due to it being stored as string
    // the data within must NOT be sorted, otherwise indices will break
    std::sort(std::begin(rps->second), std::end(rps->second), sort_by_key);
  }
  std::sort(std::begin(prescales), std::end(prescales), sort_by_key);
}



template <size_t NPATH, size_t NSEED>
int Prescaler<NPATH, NSEED>::index(const std::string &key, std::vector<std::string> &container)
{
  auto ite = std::find(std::begin(container), std::end(container), key);

  if (ite == std::end(container)) {
    container.emplace_back(key);
    return container.size() - 1;
  }
  else
    return std::distance(std::begin(container), ite);
}



template <size_t NPATH, size_t NSEED>
double Prescaler<NPATH, NSEED>::scream(int run_, int lumi_)
{
  std::cout << "WARNING: Prescaler is queried with an unexpected run " << run_ << " or lumi section " << lumi_ << ". Returning a weight of 0!" << std::endl;

  run = 0;
  lumi = 0;

  irun = 0;
  ilumi = 0;

  return 0.;
}



template <size_t NPATH, size_t NSEED>
void Prescaler<NPATH, NSEED>::check() const
{
  if (paths.size() > NPATH or seeds.size() > NSEED) {
    std::cout << "Number of paths: " << paths.size() << "\n";
    std::cout << "Number of seeds: " << seeds.size() << std::endl;
    throw std::out_of_range( "ERROR: Prescaler is too small for the amount of paths and/or seeds! Use a larger Prescaler<npath, nseed> or constrain the paths list!!" );
  }

  std::cout << "Prescaler: prescale table initialized considering the following HLT paths: (" << paths.size() << ")\n";
  for (const auto &path : paths)
    std::cout << path << "\n";
  std::cout << "\nand the following L1 seeds: (" << seeds.size() << ")\n";
  for (const auto &seed : seeds)
    std::cout << seed << "\n";
  std::cout << "\nEnsure that the bit masks provided to weight() are in this order, with the first path/seed decision at the 2^0 bit." << std::endl;
}

#endif
