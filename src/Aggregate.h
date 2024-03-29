// -*- C++ -*-
// author: afiq anuar
// short: aggregate of groups, with arbitrary indexing rule and attributes that are transformations of groups attributes

#ifndef FWK_AGGREGATE_H
#define FWK_AGGREGATE_H

#include "Group.h"

// todo: something that allows Group<A, B, C> and Group<B, A> to go into the same Aggregate
// todo: simple conversion operators is not the way to achieve this given the current model
// todo: std::visit gives the wrong result for variant<A, B> when it's casted to variant<B, A>
// todo: also wrong when casting variant<A, B> to variant<A, B, C>
// todo: i.e. there is something in variant implementation that is sensitive to type order
// todo: as for non-template base of Group with virtual functions, the issue is in the need to cast to original type before calling Group::data()
// todo: and it's not clear how to make Aggregate remember the original types without making it a class member
// todo: get rid of index masking

namespace Framework {
  template <std::size_t N, typename ...Ts>
  class Aggregate : public Group<Ts...> {
  public:
    /// constructor
    Aggregate() = delete;

    /// give directly the collections that one will be using
    /// not adding add_collection in order to simplify things a bit
    /// in that one can fix already the signature to cover only all the constructed collections
    template <typename ...Groups>
    Aggregate(const std::string &name_, int reserve_, int init, Groups &...groups);

    /// provide the indexing function
    /// ie how to go from indices in each group to an index in the aggregate
    /// which common to the entire aggregate
    /// signature: args are the refs to the N groups, and returns one std::vector<std::array<int, N>> output
    template <typename Indexer>
    void set_indexer(Indexer indexer_);

    /// reserve the space for expected number of attributes
    void reserve(int attr);

    /// add an attribute into the collection
    /// returns false upon failure to add the attribute
    /// this can happen if some data types are inconsistent
    /// returns true upon a successful addition
    /// the attributes refer to the attributes of the underlying groups
    /// syntax needs to be as expected by inquire_collection
    template <typename Function, typename ...Attributes>
    bool add_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// transform a group of internal attributes into another attribute
    /// the transformation is done element-wise on every element of held data
    /// as this adds a new attribute, its name has to be unique
    template <typename Function, typename ...Attributes>
    bool transform_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// populate the data by calling the functions provided
    void populate(long long) override;

  protected:
    /// inquire attribute of the underlying groups
    /// similar to the above, but now an array: 
    /// first collection index in v_group, second the attribute index within it
    /// assumes syntax of group::attribute
    std::array<int, 2> inquire_group(const std::string &name);

    /// ie the actual reference to the element returned by inquire_group
    /// necessarily implemented without safety...
    const typename Group<Ts...>::data_type& underlying_attribute(const std::string &name);

    /// indexing function - how to go from indices in each group to an index in the aggregate
    std::function<void()> indexer;

    /// and the indices that are made by the above
    std::vector<std::array<int, N>> v_indices;

    /// references to the groups
    std::array<std::reference_wrapper<Group<Ts...>>, N> v_group;
  };

  template <std::size_t N, typename T>
  struct aggregate_helper {};

  template <std::size_t N, template <typename, typename...> typename G, typename ...Ts>
  struct aggregate_helper<N, G<Ts...>> {
    using type = Aggregate<N, Ts...>;
  };

  // factory function, deduction guide ain't cuttin it anymore
  template <typename ...Groups>
  typename aggregate_helper<sizeof...(Groups), typename tuple_union<Groups...>::type>::type
  make_aggregate(const std::string &name, int reserve, int init, Groups &...groups);
}

#include "Aggregate.cc"

// so that many of the tuple shenanigans work
namespace std {
  template <std::size_t N, typename ...Ts>
  struct tuple_size<Framework::Aggregate<N, Ts...>> : tuple_size<Framework::Group<Ts...>> {};

  template <std::size_t I, std::size_t N, typename ...Ts>
  struct tuple_element<I, Framework::Aggregate<N, Ts...>> : tuple_element<I, typename Framework::Group<Ts...>> {};
}

#endif
