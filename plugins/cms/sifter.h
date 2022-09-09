// -*- C++ -*-
// author: afiq anuar
// short: an interface for sifting through multiple triggers based on tiers
// note: events may appear in both double and single muon datasets, and this code allows defining some priority ordering
// note: such that if double muon is considered the higher priority one, then single muon events are accepted only if it's not a double muon event

#ifndef FWK_SIFTER_H
#define FWK_SIFTER_H

#include <bitset>
#include <numeric>

#include "misc/container_util.h"
#include "misc/string_io.h"

template <std::size_t NPATH = 128>
class Sifter {
public:
  /// add trigger paths
  /// first call adds to the first tier (highest priority) and so on
  /// remembers only unique paths - the rest are ignored
  void add(const std::vector<std::string> &paths_);

  /// returns how many tiers are there
  int ntier() const { return tiers.size(); };

  /// returns paths at a given tier
  std::vector<std::string> paths_at_tier(int tier) const;

  /// decision whether or not to accept the path, based on the tier
  bool accept(int tier, const std::bitset<NPATH> &bits) const;

private:
  /// list of paths and number of paths corresponding to each tier
  std::vector<std::string> paths;
  std::vector<int> tiers;
};



template <std::size_t NPATH>
void Sifter<NPATH>::add(const std::vector<std::string> &paths_)
{
  int iadd = 0;
  for (const auto &path : paths_) {
    if (std::count(std::begin(paths), std::end(paths), path) == 0) {
      paths.emplace_back(path);
      ++iadd;
    }
  }

  if (paths.size() > NPATH)
    throw std::out_of_range( "ERROR: Sifter is too small for the amount of paths! Use a larger Sifter<npath> or constrain the paths list!!" );

  if (iadd > 0)
    tiers.emplace_back(iadd);
}



template <std::size_t NPATH>
std::vector<std::string> Sifter<NPATH>::paths_at_tier(int tier) const
{
  if (tier >= tiers.size())
    return {};

  const int idx = std::accumulate(std::begin(tiers), std::next(std::begin(tiers), tier), 0);
  return std::vector<std::string>(std::next(std::begin(paths), idx), std::next(std::begin(paths), idx + tiers[tier]));
}



template <std::size_t NPATH>
bool Sifter<NPATH>::accept(int tier, const std::bitset<NPATH> &bits) const
{
  if (tier >= tiers.size())
    return false;

  const int idx = std::accumulate(std::begin(tiers), std::next(std::begin(tiers), tier), 0);

  for (int ii = 0; ii < idx + tiers[tier]; ++ii) {
    if (bits[ii])
      return ii >= idx;
  }

  return false;
}

#endif
