// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <std::size_t N, typename ...Ts>
template <typename ...Groups>
Framework::Aggregate<N, Ts...>::Aggregate(const std::string &name_, int reserve_, int init, Groups &...groups) : 
Framework::Group<Ts...>::Group(name_, 0),
v_group{ (Group<Ts...>&) groups... }
{
  static_assert(N > 1, "ERROR: Aggregate must be made out of two or more (not necessarily unique) Groups!!");
  static_assert(sizeof...(groups) == N, "ERROR: non-matching template argument count and number of Groups in Aggregate ctor!!");
  static_assert((is_subset_of<Groups, typename Framework::Group<Ts...>::base> and ...), 
                "ERROR: the Aggregate must contain all the types contained by its underlying Groups!!");

  reserve(reserve_);
  if (init > 0)
    this->initialize(init);
  else
    this->initialize(1);
}



template <std::size_t N, typename ...Ts>
template <typename Indexer>
void Framework::Aggregate<N, Ts...>::set_indexer(Indexer indexer_)
{
  if (not indexer) {
    auto get_indexer = [indexer_] (const auto &...grps) { return indexer_(grps.get()...); };

    auto f_index = [this, get_indexer] () -> void {
      v_indices.clear();
      v_indices = std::apply(get_indexer, v_group);
    };

    indexer = std::function<void()>(f_index);
  }
}



template <std::size_t N, typename ...Ts>
void Framework::Aggregate<N, Ts...>::reserve(int attr)
{
  Group<Ts...>::reserve(attr);
}



template <std::size_t N, typename ...Ts>
template <typename Function, typename ...Attributes>
bool Framework::Aggregate<N, Ts...>::add_attribute(const std::string &attr, Function function, Attributes &&...attrs)
{
  static_assert(sizeof...(attrs) > 0, "ERROR: Aggregate::add_attribute: Aggregate attribute must be made out of at least one underlying attribute!!");

  using Traits = function_traits<decltype(function)>;
  static_assert(contained_in<typename Traits::result_type, Ts...>, 
                "ERROR: Aggregate::add_attribute: the function return type is not among the types expected by the Aggregate!!");
  static_assert(mutual_overlap<typename Traits::tuple_arg_bare_types, Group<Ts...>>, 
                "ERROR: Aggregate::add_attribute: the function argument types do not match the types expected by the Aggregate!!");

  if (this->has_attribute(attr))
    return false;

  static constexpr std::array<int, 2> no_attribute = {-1, -1};
  const bool has_all_attribute = ((inquire_group(attrs) != no_attribute) and ...);
  if (not has_all_attribute)
    return false;

  auto f_calculate = [function] (auto &val, const auto &...vals) -> void {
    val = function(vals...);
  };

  auto f_bump_duplicate = [this] (Attributes &&...attrs) -> std::array<std::array<int, 2>, sizeof...(attrs)> {
    std::array<std::array<int, 2>, sizeof...(attrs)> grp_inq = { inquire_group(attrs)... };

    // check for duplicate underlying attributes in the above
    // happens when the same groups are passed multiple times to the Aggregate
    // eg a dielectron Aggregate, with one attribute being dphi
    // and indexer taking the two eles of the highest pt
    // in this case the latter is bumped to the correct placement based on the indexer
    for (int iA1 = 0; iA1 < grp_inq.size() - 1; ++iA1) {
      for (int iA2 = iA1 + 1; iA2 < grp_inq.size(); ++iA2) {
        if (grp_inq[iA1] == grp_inq[iA2]) {
          const auto grp_name = v_group[ grp_inq[iA1][0] ].get().name;

          for (int iA3 = grp_inq[iA1][0] + 1; iA3 < v_group.size(); ++iA3) {
            if (v_group[ iA3 ].get().name == grp_name) {
              grp_inq[iA2][0] = iA3;
              break;
            }
          }
        }
      }
    }

    return grp_inq;
  };

  const auto grps = tuple_of_ref( std::make_tuple( std::ref(underlying_attribute(attrs))... ), 
                                  Traits{}, std::make_index_sequence<Traits::arity>{} );
  const std::array<std::array<int, 2>, sizeof...(attrs)> grp_inq = f_bump_duplicate(attrs...);

  // grps and grp_inq must be captured by value, since they die outside add_attribute scope
  auto f_apply = [f_calculate, this, iattr = this->v_data.size(), grps, grp_inq] (Attributes &&...attrs) -> void {
    static std::array<int, sizeof...(attrs)> attr_idx;
    for (int iE = 0; iE < std::get<int>(this->counter); ++iE) {
      attr_idx.fill(-1);
      auto single_idx = v_indices[iE];

      for (int iI = 0; iI < grp_inq.size(); ++iI)
        attr_idx[iI] = single_idx[ grp_inq[iI][0] ];

      auto args = zip_nn(grps, attr_idx);
      auto refs = std::tuple_cat(std::make_tuple(std::ref( std::get<std::vector<typename Traits::result_type>>(this->v_data[iattr])[iE] )), 
                                 args );

      std::apply(f_calculate, refs);
    }
  };

  auto f_add = [f_apply, attrs...] () -> void {
    f_apply(attrs...);
  };

  this->v_attr.emplace_back(std::make_pair(attr, std::function<void()>(f_add)));
  this->v_data.emplace_back(std::vector<typename Traits::result_type>());

  std::visit([init = this->v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, this->v_data.back());
  return true;
}



template <std::size_t N, typename ...Ts>
template <typename Function, typename ...Attributes>
bool Framework::Aggregate<N, Ts...>::transform_attribute(const std::string &attr, Function function, Attributes &&...attrs)
{
  if (Group<Ts...>::transform_attribute(attr, function, std::forward<Attributes>(attrs)...)) {
    std::visit([init = this->v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, this->v_data.back());
    return true;
  }

  return false;
}



template <std::size_t N, typename ...Ts>
void Framework::Aggregate<N, Ts...>::populate(long long)
{
  indexer();
  this->counter = int(v_indices.size());
  this->selected = std::get<int>(this->counter);

  if (std::get<int>(this->counter) > this->v_index.capacity())
    this->initialize(std::get<int>(this->counter));

  this->v_index.clear();
  for (int iD = 0; iD < std::get<int>(this->counter); ++iD)
    this->v_index.emplace_back(iD);
  this->ordered = true;

  // then run on the attributes
  for (int iD = 0; iD < this->v_data.size(); ++iD)
      this->v_attr[iD].second();
}



template <std::size_t N, typename ...Ts>
std::array<int, 2> Framework::Aggregate<N, Ts...>::inquire_group(const std::string &name)
{
  const auto iTkn = name.find("::");
  const std::string grp_name = name.substr(0, iTkn), grp_attr = name.substr(iTkn + 2, name.size());

  int iGrp = -1;
  for (int iG = 0; iG < N; ++iG) {
    if (grp_name == v_group[iG].get().name) {
      iGrp = iG;
      break;
    }
  }

  int iAttr = (iGrp == -1) ? -1 : v_group[iGrp].get().inquire(grp_attr);
  return {iGrp, iAttr};
}



template <std::size_t N, typename ...Ts>
const typename Framework::Group<Ts...>::data_type& Framework::Aggregate<N, Ts...>::underlying_attribute(const std::string &name)
{
  const auto iGA = inquire_group(name);
  return v_group[ iGA[0] ].get().data()[ iGA[1] ];
}



template <typename ...Groups>
typename Framework::aggregate_helper<sizeof...(Groups), typename tuple_union<Groups...>::type>::type
Framework::make_aggregate(const std::string &name, int reserve, int init, Groups &...groups)
{
  return typename Framework::aggregate_helper<sizeof...(Groups), typename tuple_union<Groups...>::type>::type(name, reserve, init, groups...);
}
