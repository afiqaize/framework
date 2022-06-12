// -*- C++ -*-
// author: afiq anuar
// short: collection of tools that are needed by group and friends

#ifndef FWK_HEAP_H
#define FWK_HEAP_H

#include <string>
#include <string_view>

#include <memory>
#include <utility>

#include <vector>
#include <array>
#include <tuple>
#include <variant>

#include <algorithm>
#include <numeric>

#include <type_traits>
#include <functional>
#include <stdexcept>



// credit https://stackoverflow.com/questions/18986560/check-variadic-templates-parameters-for-uniqueness
template <typename T> 
struct Type {};

template <typename ...Ts>
struct Types : Type<Ts>... {
  template <typename T>
  constexpr auto operator+(Type<T>)
  {
    if constexpr (std::is_base_of_v<Type<T>, Types>)
      return Types{};
    else
      return Types<Ts..., T>{};
  }

  constexpr std::size_t size() const { return sizeof...(Ts); }
  constexpr std::size_t sum() const { return (sizeof(Ts) + ...); }
};

template <typename ...Ts>
constexpr bool unique_types = ( (Types<>{} + ... + Type<Ts>{}).size() == sizeof...(Ts) );



// credit https://stackoverflow.com/questions/34099597/check-if-a-type-is-passed-in-variadic-template-parameter-pack
template <typename T, typename ...Ts>
constexpr bool contained_in = std::disjunction_v<std::is_same<T, Ts>...>;



// credit https://stackoverflow.com/questions/42580997/check-if-one-set-of-types-is-a-subset-of-the-other/
template <typename T, typename U>
constexpr bool is_subset_of = false;

template <template <typename, typename...> typename T, typename ...Ts, template <typename, typename...> typename U, typename ...Us>
constexpr bool is_subset_of<T<Ts...>, U<Us...>> 
= (contained_in<Ts, Us...> and ...);

template <typename T, typename U>
constexpr bool mutual_overlap = false;

template <template <typename, typename...> typename T, typename ...Ts, template <typename, typename...> typename U, typename ...Us>
constexpr bool mutual_overlap<T<Ts...>, U<Us...>> 
= (is_subset_of<T<Ts...>, U<Us...>> or is_subset_of<U<Us...>, T<Ts...>>);



// sorting type lists, first by size, and second by name
// needed so that Group<T0, ..., Tn> types all have the same data regardless of permutation
// credits https://codereview.stackexchange.com/questions/131194/selection-sorting-a-type-list-compile-time
// credits https://stackoverflow.com/questions/48723974/how-to-order-types-at-compile-time
// 1- swap types at index i and j within a type pack
template <std::size_t I, std::size_t J, typename Tuple>
struct tuple_element_swap {
  template <typename T>
  struct tuple_element_swap_impl {};

  template <std::size_t... Is>
  struct tuple_element_swap_impl<std::index_sequence<Is...>> {
    using type = std::tuple<std::tuple_element_t<Is == I ? J : Is == J ? I : Is, Tuple>...>;
  };

  using type = typename tuple_element_swap_impl<std::make_index_sequence<std::tuple_size_v<Tuple>>>::type;
};

// 2- sort template argument pack types with selection sort
template <typename Tuple, template <typename, typename> typename Compare>
struct tuple_selection_sort {
  // selection sort's "loop"
  template <std::size_t I, std::size_t J, std::size_t N, typename Pack>
  struct tuple_selection_sort_impl {
    using tuple_type = std::conditional_t<Compare<std::tuple_element_t<I, Pack>, std::tuple_element_t<J, Pack>>::value,
                                          Pack,
                                          typename tuple_element_swap<I, J, Pack>::type
                                          >;

    using type = typename tuple_selection_sort_impl<I, J + 1, N, tuple_type>::type;
  };

  template <std::size_t I, std::size_t N, typename Pack>
  struct tuple_selection_sort_impl<I, N, N, Pack> {
    using type = typename tuple_selection_sort_impl<I + 1, I + 2, N, Pack>::type;
  };

  template <std::size_t J, std::size_t N, typename Pack>
  struct tuple_selection_sort_impl<N, J, N, Pack> {
    using type = Pack;
  };

  using type = typename tuple_selection_sort_impl<0, 1, std::tuple_size_v<Tuple>, Tuple>::type;
};

// 3- a mechanism to obtain type names, using macros
template <typename T> constexpr std::string_view type_name();
template <> constexpr std::string_view type_name<void>() { return "void"; }

namespace detail 
{
  using type_name_prober = void;

  template <typename T>
  constexpr std::string_view wrapped_type_name() 
  {
#ifdef __clang__
    return __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
    return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
    return __FUNCSIG__;
#else
#error "Unsupported compiler"
#endif
  }

  constexpr std::size_t wrapped_type_name_prefix_length() { 
    return wrapped_type_name<type_name_prober>().find(type_name<type_name_prober>()); 
  }

  constexpr std::size_t wrapped_type_name_suffix_length() { 
    return wrapped_type_name<type_name_prober>().length() 
      - wrapped_type_name_prefix_length() 
      - type_name<type_name_prober>().length();
  }
}

template <typename T>
constexpr std::string_view type_name() {
  constexpr auto wrapped_name = detail::wrapped_type_name<T>();
  constexpr auto prefix_length = detail::wrapped_type_name_prefix_length();
  constexpr auto suffix_length = detail::wrapped_type_name_suffix_length();
  constexpr auto type_name_length = wrapped_name.length() - prefix_length - suffix_length;
  return wrapped_name.substr(prefix_length, type_name_length);
}

// 4- the comparator type, in this case we're sorting by size and if they're equal, by name
// other comparators can be similarly defined when there's a use case
template <typename T, typename U>
struct greater_size_or_name : std::bool_constant<(sizeof(T) == sizeof(U) ? (type_name<T>() > type_name<U>()) : sizeof(T) > sizeof(U))> {};



// for std::visit
// credit https://en.cppreference.com/w/cpp/utility/variant/visit
template <typename ...Ts> struct overload : Ts... { using Ts::operator()...; };
template <typename ...Ts> overload(Ts...) -> overload<Ts...>;
template <typename> inline constexpr bool always_false_v = false;



// credit https://stackoverflow.com/questions/7943525/is-it-possible-to-figure-out-the-parameter-type-and-return-type-of-a-lambda
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {};

template <typename Ret, typename ...Args>
struct function_traits<Ret(*)(Args...)> : public function_traits<Ret(Args...)> {};

// specialization for normal lambda
template <typename Cls, typename Ret, typename ...Args>
struct function_traits<Ret(Cls::*)(Args...) const> {
  static constexpr std::size_t arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <std::size_t I>
  using arg = typename std::tuple_element_t<I, tuple_arg_types>;

  template <std::size_t I>
  using bare_arg = typename std::tuple_element_t<I, tuple_arg_bare_types>;
};

// specialization for mutable lambda
template <typename Cls, typename Ret, typename ...Args>
struct function_traits<Ret(Cls::*)(Args...)> {
  static constexpr std::size_t arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <std::size_t I>
  using arg = typename std::tuple_element_t<I, tuple_arg_types>;

  template <std::size_t I>
  using bare_arg = typename std::tuple_element_t<I, tuple_arg_bare_types>;
};

// specialization for regular functions
template <typename Ret, typename ...Args>
struct function_traits<Ret(Args...)> {
  static constexpr std::size_t arity = sizeof...(Args);
  using result_type = Ret;
  using tuple_arg_types = std::tuple<Args...>;
  using tuple_arg_bare_types = std::tuple<std::remove_cv_t<std::remove_reference_t<Args>>...>;

  template <std::size_t I>
  using arg = typename std::tuple_element_t<I, tuple_arg_types>;

  template <std::size_t I>
  using bare_arg = typename std::tuple_element_t<I, tuple_arg_bare_types>;
};




// merging of two packs, keeping the types unique
template <typename ...Ts>
struct tuple_union {};

template <template <typename, typename...> typename Tuple, typename ...Ts, typename ...Us>
struct tuple_union<Tuple<Ts...>, Tuple<Us...>> {
  template <std::size_t I, std::size_t N, typename ...Vs>
  struct tuple_union_impl {};

  template <std::size_t I, std::size_t N, template <typename, typename...> typename T, typename ...Vs, typename ...Ws>
  struct tuple_union_impl<I, N, T<Vs...>, T<Ws...>> {
    using ith = std::conditional_t<contained_in<std::tuple_element_t<I, T<Ws...>>, Vs...>,
                                   T<Vs...>,
                                   T<std::tuple_element_t<I, T<Ws...>>, Vs...>>;

    using type = typename tuple_union_impl<I + 1, N, ith, T<Ws...>>::type;
  };

  template <std::size_t N, template <typename, typename...> typename T, typename ...Vs, typename ...Ws>
  struct tuple_union_impl<N, N, T<Vs...>, T<Ws...>> {
    using type = T<Vs...>;
  };

  using type = typename tuple_union_impl<0, std::tuple_size_v<Tuple<Us...>>, Tuple<Ts...>, Tuple<Us...>>::type;
};

template <typename T, typename U, typename ...Ts>
struct tuple_union<T, U, Ts...> {
  using type = typename tuple_union<typename tuple_union<T, U>::type, Ts...>::type;
};



// make a tuple to references to Group data per function arg types
template <typename Tuple, typename Traits, std::size_t ...Is>
auto tuple_of_ref(const Tuple &tuple, Traits, std::index_sequence<Is...>)
{
  return std::make_tuple( std::ref(std::get<std::vector<typename Traits::template bare_arg<Is>>>( std::get<Is>(tuple) ))... );
}



// credit https://stackoverflow.com/questions/11322095/how-to-make-a-function-that-zips-two-tuples-in-c11-stl
// given tuples vec = {vec1, vec2, ... vecN} and idx = {i1, i2, ... iN}, where vec is anything that has the subscript operator []
// zip_index returns another tuple of references tup = {vec1[i1], vec2[i2], ... vecN[iN]}
// which can be passed to functions with std::apply
// _nn_ in name refers to "mapping n to n"
template <template <typename ...> typename T1, typename ...T1s, template <typename, std::size_t> typename A2, typename T2, std::size_t ...N>
constexpr auto zip_nn_helper(const T1<T1s...> &t1, const A2<T2, sizeof...(N)> &a2, std::index_sequence<N...>) 
-> decltype(std::forward_as_tuple( std::get<N>(t1)[ std::get<N>(a2) ]... ))
{
  return std::forward_as_tuple( std::get<N>(t1)[ std::get<N>(a2) ]... );
}

template <template <typename ...> typename T1, typename ...T1s, template <typename, std::size_t> typename A2, typename T2, std::size_t N>
constexpr auto zip_nn(const T1<T1s...> &t1, const A2<T2, N> &a2) 
-> decltype( zip_nn_helper(t1, a2, std::make_index_sequence<N>{}) ) 
{
  static_assert(sizeof...(T1s) == N, "ERROR: the tuple sizes must be the compatible!");
  return zip_nn_helper(t1, a2, std::make_index_sequence<N>{});
}

// a version implementing tup = {vec[i1], vec[i2], ... vec[N]}
// FIXME funny that _nn_ works in Aggregate:add_attribute swapping the forward_as_tuple to make_tuple
// FIXME but the same is not true for _1n_ in Group::transform_attribute
// FIXME make_tuple( std::ref(...)... ) works, but forward_as_tuple is chosen as it's more homogenous
template <typename T1, template <typename, std::size_t> typename A2, typename T2, std::size_t ...N>
constexpr auto zip_1n_helper(T1 &t1, const A2<T2, sizeof...(N)> &a2, std::index_sequence<N...>) 
-> decltype(std::forward_as_tuple( t1[ std::get<N>(a2) ]... ))
{
  return std::forward_as_tuple( t1[ std::get<N>(a2) ]... );
}

template <typename T1, template <typename, std::size_t> typename A2, typename T2, std::size_t N>
constexpr auto zip_1n(T1 &t1, const A2<T2, N> &a2) 
-> decltype( zip_1n_helper(t1, a2, std::make_index_sequence<N>{}) ) 
{
  return zip_1n_helper(t1, a2, std::make_index_sequence<N>{});
}



// https://stackoverflow.com/questions/670308/alternative-to-vectorbool
class boolean {
public:
  boolean(): value() {}
  boolean(bool value_) : value(value_) {}

  operator bool() const {return value;}

  /// the following operators are to allow bool *b = &v[0]; (v is a vector here)
  bool* operator& () { return &value; }
  const bool* operator& () const { return &value; }

private:
  bool value;
};



// literals
using namespace std::string_literals;
using namespace std::string_view_literals;

#endif
