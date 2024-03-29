// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename ...Ts>
Framework::Group<Ts...>::Group(const std::string &name_, int counter_) :
name(name_),
counter(counter_),
selected(counter_),
v_index(this),
ordered(false)
{
  if (std::get<int>(counter) > 0) {
    v_index.resize(std::get<int>(counter));
    std::iota(std::begin(v_index), std::end(v_index), 0);
  }
}



template <typename ...Ts>
int Framework::Group<Ts...>::n_elements() const noexcept
{
  return selected;
}



template <typename ...Ts>
int Framework::Group<Ts...>::size() const noexcept
{
  return selected;
}



template <typename ...Ts>
const int& Framework::Group<Ts...>::ref_to_n_elements() const noexcept
{
  return selected;
}



template <typename ...Ts>
const int& Framework::Group<Ts...>::ref_to_size() const noexcept
{
  return selected;
}



template <typename ...Ts>
int& Framework::Group<Ts...>::mref_to_n_elements() noexcept
{
  return selected;
}



template <typename ...Ts>
int& Framework::Group<Ts...>::mref_to_size() noexcept
{
  return selected;
}



template <typename ...Ts>
int Framework::Group<Ts...>::n_attributes() const noexcept
{
  return v_attr.size() + v_alias.size();
}



template <typename ...Ts>
bool Framework::Group<Ts...>::has_attribute(const std::string &name) const noexcept
{
  return inquire(name) != -1;
}



template <typename ...Ts>
void Framework::Group<Ts...>::reserve(int attr)
{
  v_attr.reserve(attr);
  v_data.reserve(attr);
}



template <typename ...Ts>
template <typename Function, typename ...Attributes>
bool Framework::Group<Ts...>::transform_attribute(const std::string &attr, Function function, Attributes &&...attrs)
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::transform_attribute requires some attributes to be provided!!");

  using Traits = function_traits<decltype(function)>;
  static_assert(contained_in<typename Traits::result_type, Ts...>, 
                "ERROR: Group::transform_attribute: the function return type is not among the types expected by the Group!!");
  static_assert(mutual_overlap<typename Traits::tuple_arg_bare_types, Group<Ts...>>, 
                "ERROR: Group::transform_attribute: the function argument types do not match the types expected by the Group!!");

  if (has_attribute(attr))
    return false;

  if (auto iA = (has_attribute(attrs) and ...); not iA)
    throw std::invalid_argument( "ERROR: Group::transform_attribute: some of the requested attributes are not within the group!!" );

  // note on the structure
  // the multistage function execution is needed because
  // f_loop runs on actual data and produces the result
  // f_apply matches the attribute indices and forward the relevant refs to f_loop
  // all the functions need to be copied instead of referred 
  // due to scoping and/or lambda vs function pointer support

  auto f_loop = [function, this] (auto &vec, const auto &...vecs) -> void {
    // pack capture not in c++17, so unsure how to do this with std::visit
    // types are known, so good old get_if works fine
    auto ic = std::get_if<int>(&this->counter);
    auto uc = std::get_if<uint>(&this->counter);
    int current_counter = (ic != nullptr) ? *ic : *uc;

    for (int iE = 0; iE < current_counter; ++iE)
      vec[iE] = function(vecs[iE]...);
  };

  const std::array<int, sizeof...(attrs)> iattrs = {inquire(attrs)...};
  retype_per_function(zip_1n(v_data, iattrs), Traits{}, std::make_index_sequence<Traits::arity>{});

  auto f_apply = [f_loop, this, iattr = v_data.size(), iattrs] () -> void {
    auto refs = std::tuple_cat(std::make_tuple(std::ref( std::get<std::vector<typename Traits::result_type>>(v_data[iattr]) )), 
                               tuple_of_ref( zip_1n(v_data, iattrs), Traits{}, std::make_index_sequence<Traits::arity>{}) );

    std::apply(f_loop, refs);
  };

  v_attr.emplace_back(std::make_pair(attr, std::function<void()>(f_apply)));
  v_data.emplace_back(std::vector<typename Traits::result_type>());

  return true;
}



template <typename ...Ts>
bool Framework::Group<Ts...>::alias_attribute(const std::string &alias, const std::string &attr)
{
  if (has_attribute(alias))
    return false;

  if (auto iA = inquire(attr); iA == -1)
    throw std::invalid_argument( "ERROR: Group::alias_attribute: the requested attribute to be aliased is not present in the Group!!" );
  else
    v_alias.emplace_back(alias, iA);

  return true;
}



template <typename ...Ts>
std::vector<std::string> Framework::Group<Ts...>::attributes() const
{
  std::vector<std::string> v_attr_name;
  for (const auto &attr : v_attr)
    v_attr_name.emplace_back(attr.first);
  for (const auto &attr : v_alias)
    v_attr_name.emplace_back(attr.first);
  return v_attr_name;
}



template <typename ...Ts>
const std::vector<typename Framework::Group<Ts...>::data_type>& Framework::Group<Ts...>::data() const
{
  return v_data;
}



template <typename ...Ts>
const typename Framework::Group<Ts...>::data_type& Framework::Group<Ts...>::operator()(const std::string &name) const
{
  auto iA = inquire(name);
  if (iA == -1)
    throw std::invalid_argument( "ERROR: Group::get: requested attribute " + name + " is not within the group!!" );

  return v_data[iA];
}



template <typename ...Ts>
const typename Framework::Group<Ts...>::data_type& Framework::Group<Ts...>::operator()(int iattr) const
{
  return v_data[iattr];
}



template <typename ...Ts>
typename Framework::Group<Ts...>::data_type& Framework::Group<Ts...>::mref_to_attribute(const std::string &name)
{
  return const_cast<typename Framework::Group<Ts...>::data_type &>( (*const_cast<const Framework::Group<Ts...>*>(this))(name) );
}



template <typename ...Ts>
typename Framework::Group<Ts...>::data_type& Framework::Group<Ts...>::mref_to_attribute(int iattr)
{
  return const_cast<typename Framework::Group<Ts...>::data_type &>( (*const_cast<const Framework::Group<Ts...>*>(this))(iattr) );
}



template <typename ...Ts>
template <typename T>
const std::vector<T>& Framework::Group<Ts...>::get(const std::string &name) const
{
  static_assert(contained_in<T, Ts...>, "ERROR: Group::get: called with a type not among the types of by the Group!!");
  return std::get<std::vector<T>>((*this)(name));
}



template <typename ...Ts>
template <typename T>
const std::vector<T>& Framework::Group<Ts...>::get(int iattr) const
{
  static_assert(contained_in<T, Ts...>, "ERROR: Group::get: called with a type not among the types of by the Group!!");
  return std::get<std::vector<T>>((*this)(iattr));
}



template <typename ...Ts>
template <typename T>
const std::vector<T>* Framework::Group<Ts...>::get_if(const std::string &name) const noexcept
{
  auto iA = inquire(name);
  if (iA == -1)
    return nullptr;

  return std::get_if<std::vector<T>>(&v_data[iA]);
}



template <typename ...Ts>
template <typename T>
const std::vector<T>* Framework::Group<Ts...>::get_if(int iattr) const noexcept
{
  if (iattr >= 0 and iattr < n_attributes())
    return std::get_if<std::vector<T>>(&v_data[iattr]);

  return nullptr;
}



template <typename ...Ts>
template <typename T>
const T& Framework::Group<Ts...>::get(const std::string &name, int index) const
{
  return this->get<T>(name)[index];
}



template <typename ...Ts>
template <typename T>
const T& Framework::Group<Ts...>::get(int iattr, int index) const
{
  return this->get<T>(iattr)[index];
}



template <typename ...Ts>
template <typename T>
T Framework::Group<Ts...>::convert(const std::string &name, int index) const
{
  return std::visit([index] (auto &&arg) -> T { return arg[index]; }, (*this)(name));
}



template <typename ...Ts>
template <typename T>
T Framework::Group<Ts...>::convert(int iattr, int index) const
{
  return std::visit([index] (auto &&arg) -> T { return arg[index]; }, (*this)(iattr));
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::indices() const
{
  return v_index;
}



template <typename ...Ts>
const typename Framework::Group<Ts...>::idxs& Framework::Group<Ts...>::ref_to_indices() const
{
  return v_index;
}



template <typename ...Ts>
int& Framework::Group<Ts...>::operator[](int idx)
{
  return v_index[idx];
}



template <typename ...Ts>
const int& Framework::Group<Ts...>::operator[](int idx) const
{
  return v_index[idx];
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::iter Framework::Group<Ts...>::begin() noexcept
{
  return v_index.begin();
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::citer Framework::Group<Ts...>::begin() const noexcept
{
  return v_index.begin();
}


template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::citer Framework::Group<Ts...>::cbegin() const noexcept
{
  return v_index.cbegin();
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::iter Framework::Group<Ts...>::end() noexcept
{
  return v_index.end();
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::citer Framework::Group<Ts...>::end() const noexcept
{
  return v_index.end();
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs::citer Framework::Group<Ts...>::cend() const noexcept
{
  return v_index.cend();
}



template <typename ...Ts>
template <typename Idx>
bool Framework::Group<Ts...>::update_indices(Idx &&v_idx)
{
  static_assert(std::is_same_v<std::decay_t<Idx>, Framework::Group<Ts...>::idxs>, "ERROR: Group::update_indices can only be called with indices type!!");

  v_index = std::forward<Idx>(v_idx);
  selected = v_index.size();
  ordered = false;

  return selected > 0;
}



template <typename ...Ts>
template <typename Function, typename IdxAttr, typename ...Attributes>
void Framework::Group<Ts...>::iterate(Function function, const IdxAttr &idxs_or_attr, Attributes &&...attrs) const
{
  if constexpr (std::is_same_v<IdxAttr, Framework::Group<Ts...>::idxs>) {
      static_assert(sizeof...(attrs) > 0, "ERROR: Group::iterate makes no sense without specifying attributes!!");

      auto iA = (has_attribute(attrs) and ...);
      if (!iA)
        throw std::invalid_argument( "ERROR: Group::iterate: some of the requested attributes are not within the group!!" );

      if (idxs_or_attr.ref != this)
        throw std::invalid_argument( "ERROR: Group::iterate: requesting iteration over index set of another Group is nonsensical!!" );

      std::visit([&function, &idxs_or_attr, this] (const auto &...vec) {
          for (auto &&idx : idxs_or_attr)
            function(vec[idx]...);
        }, v_data[inquire(attrs)]...);
    }
  else {
    if (auto iA = ((has_attribute(idxs_or_attr) and has_attribute(attrs)) and ...); not iA)
      throw std::invalid_argument( "ERROR: Group::iterate: some of the requested attributes are not within the group!!" );

    std::visit([&function, this] (const auto &...vec) {
        for (auto &&idx : v_index)
          function(vec[idx]...);
      }, v_data[inquire(idxs_or_attr)], v_data[inquire(attrs)]...);
  }
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                       Compare compare, Attributes &&...attrs) const
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Group::filter makes no sense without specifying attributes!!");

  if (auto iA = (has_attribute(attrs) and ...); not iA)
    throw std::invalid_argument( "ERROR: Group::filter some of the requested attributes are not within the group!!" );

  if (v_idx.ref == this and v_idx.empty())
    return v_idx;

  return filter_helper( (v_idx.ref != this) ? v_index : v_idx, compare, inquire(attrs)... );
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_less(const std::string &name, Number value, 
                                                                            const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data < value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_less_equal(const std::string &name, Number value, 
                                                                                  const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data <= value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_greater(const std::string &name, Number value, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data > value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_greater_equal(const std::string &name, Number value, 
                                                                                     const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data >= value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_equal(const std::string &name, Number value, 
                                                                             const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data == value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_not(const std::string &name, Number value, 
                                                                           const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return data != value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_absolute_equal(const std::string &name, Number value, 
                                                                                      const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return std::abs(data) == value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_absolute_not(const std::string &name, Number value, 
                                                                                    const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {return std::abs(data) != value;}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_bit_and(const std::string &name, Number value, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&value] (auto &data) {
      if constexpr(std::is_integral_v<std::decay_t<decltype(data)>> and 
                   std::is_integral_v<std::decay_t<Number>>)
                    return (data & value);
      else
        return false;
    }, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_in(const std::string &name, Number min, Number max, 
                                                                          const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&min, &max] (auto &data) {return (data > min and data < max);}, name);
}



template <typename ...Ts>
template <typename Number>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_out(const std::string &name, Number min, Number max, 
                                                                           const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return filter(v_idx, [&min, &max] (auto &data) {return (data < min and data > max);}, name);
}



template <typename ...Ts>
template <typename Compare>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                     Compare compare, const std::string &name) const
{
  auto iA = inquire(name);
  if (iA == -1)
    throw std::invalid_argument( "ERROR: Group::sort: some of the requested attributes are not within the group!!" );

  if (v_idx.ref == this and v_idx.size() < 2)
    return v_idx;

  return sort_helper( (v_idx.ref != this) ? v_index : v_idx, compare, iA );
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_ascending(const std::string &name, 
                                                                               const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (p1.second < p2.second); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_descending(const std::string &name, 
                                                                                const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (p1.second > p2.second); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_absolute_ascending(const std::string &name, 
                                                                                        const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (std::abs(p1.second) < std::abs(p2.second)); }, name);
}



template <typename ...Ts>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_absolute_descending(const std::string &name, 
                                                                                         const typename Framework::Group<Ts...>::idxs &v_idx) const
{
  return sort(v_idx, [] (const auto &p1, const auto &p2) { return (std::abs(p1.second) > std::abs(p2.second)); }, name);
}



template <typename ...Ts>
int Framework::Group<Ts...>::inquire(const std::string &name) const noexcept
{
  for (int iA = 0; iA < v_alias.size(); ++iA) {
    if (auto &[alias, idx] = v_alias[iA]; alias == name)
      return idx;
  }

  for (int iA = 0; iA < v_attr.size(); ++iA) {
    if (v_attr[iA].first == name)
      return iA;
  }

  return -1;
}



template <typename ...Ts>
void Framework::Group<Ts...>::reorder()
{
  if (not ordered) {
    for (int iS = 0; iS < selected; ++iS) {
      if (iS != v_index[iS]) {
        for (auto &dat : v_data)
          std::visit([iS, iI = v_index[iS]] (auto &vec) {std::swap(vec[iS], vec[iI]);}, dat);

        v_index[iS] = iS;
      }
    }

    ordered = true;
  }
}



template <typename ...Ts>
template <template <typename, typename...> typename Other, typename ...Us>
Framework::Group<Ts...>::operator Other<Us...>() const
{
  static_assert(is_subset_of<Group<Ts...>, Other<Us...>>, "ERROR: Group conversion must not be narrowing!!!");
  return reinterpret_cast<Other<Us...> &>(*this);
}



template <typename ...Ts>
void Framework::Group<Ts...>::initialize(int init)
{
  v_index.reserve(init);

  for (auto &dat : v_data)
    std::visit([init] (auto &vec) {vec.reserve(init); vec.clear();}, dat);
}



template <typename ...Ts>
template <typename Number>
void Framework::Group<Ts...>::retype(typename Framework::Group<Ts...>::data_type &dat)
{
  if constexpr (contained_in<Number, Ts...>) {
      if (std::get_if<std::vector<Number>>(&dat) == nullptr) {
        dat = std::vector<Number>();
        std::visit([init = v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, dat);
      }
    }
}



template <typename ...Ts>
template <typename Tuple, typename Traits, std::size_t ...Is>
void Framework::Group<Ts...>::retype_per_function(const Tuple &tuple, Traits, std::index_sequence<Is...>)
{
  (retype<typename Traits::template bare_arg<Is>>(std::get<Is>(tuple)), ...);
}



template <typename ...Ts>
template <typename Compare, typename ...Attributes>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::filter_helper(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                              Compare &compare, Attributes &&...attrs) const
{
  typename Framework::Group<Ts...>::idxs v_pass(this);
  std::visit([&v_pass, &v_idx, &compare] (const auto &...vec) {
      for (auto &&index : v_idx) {
        if (compare(vec[index]...))
          v_pass.emplace_back(index);
      }
    }, v_data[attrs]...);

  return v_pass;
}



template <typename ...Ts>
template <typename Compare>
typename Framework::Group<Ts...>::idxs Framework::Group<Ts...>::sort_helper(const typename Framework::Group<Ts...>::idxs &v_idx, 
                                                                            Compare &compare, int attr) const
{
  typename Framework::Group<Ts...>::idxs v_sort(this);
  std::visit([&v_sort, &v_idx, &compare] (const auto &vec) {
      using VT = typename std::decay_t<decltype(vec)>::value_type;
      std::vector<std::pair<int, VT>> v_zip;
      for (auto &&index : v_idx)
        v_zip.emplace_back(index, vec[index]);

      std::sort(std::begin(v_zip), std::end(v_zip), compare);

      for (auto &&zip : v_zip)
        v_sort.emplace_back(zip.first);
    }, v_data[attr]);

  return v_sort;
}
