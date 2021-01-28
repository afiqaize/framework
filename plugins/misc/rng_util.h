// -*- C++ -*-
// author: afiq anuar
// short: functions related to rng and their use to generate stuff

#ifndef FWK_RNG_UTIL_H
#define FWK_RNG_UTIL_H

#include <algorithm>
#include <array>
#include <functional>
#include <random>

#include "string_io.h"

/// make a seeded rng
/// seeding strategy is mt19937-specific so maybe don't use with anything else
/// credit https://stackoverflow.com/a/444614/13007174
template <typename T = std::mt19937_64>
T random_generator() {
  auto constexpr seed_bytes = sizeof(typename T::result_type) * T::state_size;
  auto constexpr seed_len = seed_bytes / sizeof(std::seed_seq::result_type);
  auto seed = std::array<std::seed_seq::result_type, seed_len>();
  auto dev = std::random_device();
  std::generate_n(std::begin(seed), seed_len, std::ref(dev));
  auto seed_seq = std::seed_seq(begin(seed), end(seed));
  return T{seed_seq};
}



/// make a random string from an allowed character list
/// called variable_name because mostly intended to be used as attribute name generator
/// default length is the maximum length that would fit within SSO (tested at NAF with g++ 831 and clang++ 801)
/// credit https://stackoverflow.com/a/444614/13007174
std::string random_variable_name(std::size_t len = 15)
{
  static constexpr auto chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_";
  static auto rng = random_generator<>();
  auto dist = std::uniform_int_distribution{{}, std::strlen(chars) - 1};
  auto result = std::string(len, '\0');
  std::generate_n(std::begin(result), len, [&dist] () { return chars[dist(rng)]; });
  return result;
}


#endif
