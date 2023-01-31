// -*- C++ -*-
// author: afiq anuar
// short: a listing of free convenience functions for use in the framework
// note: this is basically an extention of Heap, with no clear distinction on what goes where as of now

#ifndef FWK_FUNCTION_UTIL_H
#define FWK_FUNCTION_UTIL_H

#include "../src/Heap.h"
#include <limits>
#include <bitset>

/// a function that performs a simple copy of the argument
template <typename T = float>
T identity(T t)
{
  return t;
}



/// a or b or c or ... in function form
/// call any_of<N> to get a function that takes N bools and return the OR of them all
template <typename ...Bools>
bool any_of_impl(Bools ...bools)
{
  return (bools or ...);
}



template <std::size_t ...N>
auto any_of_helper(std::index_sequence<N...>) -> bool(*)(typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...)
{
  return any_of_impl<typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...>;
}



template <std::size_t N = 1>
auto any_of() -> decltype(any_of_helper(std::make_index_sequence<N>{}))
{
  return any_of_helper(std::make_index_sequence<N>{});
}



/// a and b and c and ... in function form
/// call all_of<N> to get a function that takes N bools and return the AND of them all
template <typename ...Bools>
bool all_of_impl(Bools ...bools)
{
  return (bools and ...);
}



template <std::size_t ...N>
auto all_of_helper(std::index_sequence<N...>) -> bool(*)(typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...)
{
  return all_of_impl<typename std::tuple_element_t<N, std::array<boolean, sizeof...(N)>>...>;
}



template <std::size_t N = 1>
auto all_of() -> decltype(all_of_helper(std::make_index_sequence<N>{}))
{
  return all_of_helper(std::make_index_sequence<N>{});
}



/// convert the bools into a bitset
template <std::size_t N, typename ...Bools>
std::bitset<N> to_bitset_impl(Bools ...bools)
{
  static_assert(N >= sizeof...(bools), "ERROR: to_bitset needs to be at least as large as the provided bools!!");

  std::array<decltype((bools and ...)), sizeof...(bools)> vals = { bools... };
  std::bitset<N> bits;

  for (int ib = 0; ib < vals.size(); ++ib)
    if (vals[ib])
      bits.set(ib);

  return bits;
}



template <typename T, std::size_t N, std::size_t ...A>
auto to_bitset_helper(std::index_sequence<A...>) -> std::bitset<N>(*)(typename std::tuple_element_t<A, std::array<T, sizeof...(A)>>...)
{
  return to_bitset_impl<N, typename std::tuple_element_t<A, std::array<T, sizeof...(A)>>...>;
}



template <typename T, std::size_t A, std::size_t N = 128>
auto to_bitset() -> decltype(to_bitset_helper<T, N>(std::make_index_sequence<A>{}))
{
  return to_bitset_helper<T, N>(std::make_index_sequence<A>{});
}



/// implementation of the apply_to<N, F>, refer to that for more info
template <typename Number, std::size_t ...I, typename Function, std::size_t ...N>
auto apply_to_impl(std::index_sequence<I...>, Function function, std::index_sequence<N...>)
{
  using Traits = function_traits<decltype(function)>;
  auto f_apply = [function] (typename std::tuple_element_t<N, std::array<Number, sizeof...(N)>> ...args) -> typename Traits::result_type { 
    return function( std::get<sizeof...(N) - Traits::arity + I>(std::array<Number, sizeof...(N)>{args...})... );
  };

  return f_apply;
}



/// helper of apply_to<N, F>, refer to that for more info
template <typename Number, typename Function, std::size_t ...N>
auto apply_to_helper(Function function, std::index_sequence<N...>)
{
  using Traits = function_traits<decltype(function)>;
  return apply_to_impl<typename Traits::result_type>(std::make_index_sequence<Traits::arity>{}, function, std::make_index_sequence<sizeof...(N)>{});
}



/// a function that outputs a function that applies another function on its last F arguments
/// let ff be a function taking F arguments
/// then a call of apply_to<N>(ff) returns a function that takes (N + 1) * F arguments
/// whose result is the same as calling ff on the last F arguments, ignoring the first NF arguments
/// there is no reason why would apply_to<0> not work; it's just prevented as no use case is envisioned for it
/// being that this is written mainly to make index masking a touch easier
template <std::size_t N, typename Function>
auto apply_to(Function function)
{
  static_assert(N > 0, "ERROR: apply_to is not callable with N = 0, as in this case no masking is necessary!!");
  using Traits = function_traits<decltype(function)>;
  return apply_to_helper<typename Traits::result_type>(function, std::make_index_sequence<(N + 1) * Traits::arity>{});
}



/// variadic min and max
/// credit https://stackoverflow.com/a/63330289/13007174
template <typename Arg1, typename Arg2, typename... Args>
constexpr auto min(Arg1 &&arg1, Arg2 &&arg2, Args &&...args)
{
  if constexpr (sizeof...(args) == 0)
                 return arg1 < arg2 ? arg1 : arg2;
  else
    return min(min(arg1, arg2), args...);
}



template <typename Arg1, typename Arg2, typename... Args>
constexpr auto max(Arg1 &&arg1, Arg2 &&arg2, Args &&...args)
{
  if constexpr (sizeof...(args) == 0)
                 return arg1 > arg2 ? arg1 : arg2;
  else
    return max(max(arg1, arg2), args...);
}



template <std::size_t N>
constexpr std::size_t last_n(std::size_t bits)
{
  static_assert(N < sizeof(std::size_t), "ERROR: this function makes no sense for too large N!");
  return bits & std::bitset<N>(std::numeric_limits<std::size_t>::max()).to_ullong();
}

#endif
