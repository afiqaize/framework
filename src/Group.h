// -*- C++ -*-
// author: afiq anuar
// short: group of objects i.e. sets of attributes e.g. a particle has pt eta phi...

#ifndef FWK_GROUP_H
#define FWK_GROUP_H

// note: attributes that are self-referencing e.g. GenPart_motherIdx can't be handled by iterate(); for this one needs to use operator()

#include "Heap.h"
#include "Indices.h"

namespace Framework {
  template <typename ...Ts>
  class Group {
    static_assert(unique_types<Ts...>, "ERROR: a Group must be initialized with unique types!");
    static_assert(not contained_in<bool, Ts...>, "ERROR: Group relies on a few features that are incompatible with standard C++ "
                  "when bool is among the included types. For boolean attributes please use the 'boolean' type instead, which is a drop-in "
                  "replacement provided precisely to avoid this quirk of the standard.");

  public:
    using base = Group<Ts...>;
    using idxs = Indices<base>;
    using sorted = typename tuple_selection_sort<std::tuple<Ts...>, greater_size_or_name>::type;

    template <typename T>
    struct data_impl {};

    template <std::size_t... Is>
    struct data_impl<std::index_sequence<Is...>> {
      using type = std::variant<std::vector<std::tuple_element_t<Is, sorted>>...>;
    };
    using data_type = typename data_impl<std::make_index_sequence<sizeof...(Ts)>>::type;

    /// constructor
    Group() = delete;
    Group(const std::string &name_, int counter_);

    /// such that if (Group) is a well-defined expression
    /// and to no longer need to do if (Group.n_elements)
    explicit operator bool() const { return selected > 0; };

    /// number of currently held elements
    /// these two are absolutely identical
    int n_elements() const noexcept;
    int size() const noexcept;

    /// ref instead of copy of the above
    const int& ref_to_n_elements() const noexcept;
    const int& ref_to_size() const noexcept;

    /// a mutable ref version
    /// can't be const if it's to be used to write TTree...
    /// might be worth considering to write TTree using copies rather than in-place references?
    int& mref_to_n_elements() noexcept;
    int& mref_to_size() noexcept;

    /// number of currently held attributes
    int n_attributes() const noexcept;

    /// as it says on the tin
    bool has_attribute(const std::string &name) const noexcept;

    /// reserve the space for expected number of attributes
    void reserve(int attr);

    /// transform a group of internal attributes into another attribute
    /// the transformation is done element-wise on every element of held data
    /// as this adds a new attribute, its name has to be unique
    template <typename Function, typename ...Attributes>
    bool transform_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// list of attributes
    std::vector<std::string> attributes() const;

    /// reference to container of elements
    const std::vector<data_type>& data() const;

    /// reference to single attribute array - variant version
    const data_type& operator()(const std::string &name) const;

    /// overload the above for when the attribute index is known e,g. from inquire()
    /// deliberately written without checks for invalid index, so use with care
    const data_type& operator()(int iattr) const;

    /// mutable version of operator()
    /// intended for use by Tree only
    data_type& mref_to_attribute(const std::string &name);

    /// known attribute index overload
    data_type& mref_to_attribute(int iattr);

    /// reference to single attribute array - typed version
    template <typename T>
    const std::vector<T>& get(const std::string &name) const;

    /// known attribute index overload
    template <typename T>
    const std::vector<T>& get(int iattr) const;

    /// like the above, but returns a nullptr if attribute doesn't exist, or type is wrong, etc
    template <typename T>
    const std::vector<T>* get_if(const std::string &name) const noexcept;

    /// known attribute index overload
    /// similar to other known index overloads
    template <typename T>
    const std::vector<T>* get_if(int iattr) const noexcept;

    /// reference to single element in an attribute - typed version
    /// equivalent to get<T>(attr)[index]
    /// i.e. gives the nth element as per the current Group state accounting for previous update_indices calls 
    template <typename T>
    const T& get(const std::string &name, int index) const;

    /// known attribute index overload
    template <typename T>
    const T& get(int iattr, int index) const;

    /// the associated indices to be used with the above
    idxs indices() const;

    /// ref instead of copy of the above
    const idxs& ref_to_indices() const;

    /// a few index access utilities ala Indices
    int& operator[](int idx);

    const int& operator[](int idx) const;

    typename idxs::iter begin() noexcept;

    typename idxs::citer begin() const noexcept;

    typename idxs::citer cbegin() const noexcept;

    typename idxs::iter end() noexcept;

    typename idxs::citer end() const noexcept;

    typename idxs::citer cend() const noexcept;

    /// update the indices with another set e.g. output of filter/sort
    /// no checking if the indices are actually legit
    bool update_indices(const idxs &v_idx);

    /// element iterator taking a function and runs the visitor over it
    /// can be type-dependent or otherwise
    /// second arg can either be the indices or the name of the first attribute to iterate over
    /// if indices, then iterate over those indices only
    /// if an attribute name, then iteration is done over all currently held elements
    /// third args onwards are optional extra attributes to iterate over
    template <typename Function, typename IdxAttr, typename ...Attributes>
    void iterate(Function function, const IdxAttr &idxs_or_attr, Attributes &&...attrs) const;

    /// filter the elements in the collection by some criteria on a given attribute
    /// custom filter needs a function returning a bool and taking two args of type Number
    /// the first arg is the indices to restrict the filtering over
    /// if it equals idxs(), no restriction is imposed i.e. the current Group indices are used
    /// second arg is the function running on (single elements of) the attributes
    /// following args are the attributes themselves
    /// returns the indices after filtering
    /// this method DOES NOT update the index in-place
    template <typename Compare, typename ...Attributes>
    idxs filter(const idxs &v_idx, Compare compare, Attributes &&...attrs) const;

    /// common filters
    template <typename Number>
    idxs filter_less(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_less_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_greater(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_greater_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_not(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_absolute_equal(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_absolute_not(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_bit_and(const std::string &name, Number value, const idxs &v_idx = idxs()) const;

    /// both are min and max exclusive
    template <typename Number>
    idxs filter_in(const std::string &name, Number min, Number max, const idxs &v_idx = idxs()) const;

    template <typename Number>
    idxs filter_out(const std::string &name, Number min, Number max, const idxs &v_idx = idxs()) const;

    /// sort the elements in the collection by a given attribute
    /// custom sorter needs a function returning a bool and taking two args, both of std::pair<int, decltype(data)>
    /// FIXME prepare a more convenient implementation
    /// returns the sorted indices
    template <typename Compare>
    idxs sort(const idxs &v_idx, Compare compare, const std::string &name) const;

    /// common sorts
    idxs sort_ascending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_descending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_absolute_ascending(const std::string &name, const idxs &v_idx = idxs()) const;

    idxs sort_absolute_descending(const std::string &name, const idxs &v_idx = idxs()) const;

    /// returns the index where an attribute occurs
    int inquire(const std::string &name) const noexcept;

    /// populate the Group data
    virtual void populate(long long entry) = 0;

    /// reorder the group data such that selected elements occur in front
    /// selected elements are those whose index is in v_index
    void reorder();

    /// conversion to other Group types, provided the conversion is not narrowing
    /// should be safe... right?
    template <template <typename, typename...> typename Other, typename ...Us>
    operator Other<Us...>() const;

    /// name of the group
    std::string name;

  protected:
    /// this method ensures that all attributes have the proper capacity
    void initialize(int init);

    /// update the type of the attribute
    template <typename Number>
    void retype(data_type &dat);

    /// helper for transform_attribute
    /// where we need to call different instantiations of retype
    template <typename Tuple, typename Traits, std::size_t ...Is>
    void retype_per_function(const Tuple &tuple, Traits, std::index_sequence<Is...>);

    /// helper that actually does the filtering
    template <typename Compare, typename ...Attributes>
    idxs filter_helper(const idxs &v_idx, Compare &compare, Attributes &&...attrs) const;

    /// helper that actually does the sorting
    template <typename Compare>
    idxs sort_helper(const idxs &v_idx, Compare &compare, int attr) const;

    /// element counter before prefiltering
    /// a variant purely to suppress root's error message
    /// counter is still not expected to go beyond int, however
    std::variant<int, uint> counter;

    /// element counter after prefiltering i.e. v_index.size()
    int selected;

    /// element indices in the group
    idxs v_index;

    /// a flag indicating that the data is ordered i.e. indexed as 0, 1, ..., N in memory
    /// needed to optimize reorder(), which is called by Tree for all branches
    bool ordered;

    /// register of the currently available attributes
    /// first string is attribute alias
    /// second function is for the element-wise transformation from other attributes
    std::vector<std::pair<std::string, std::function<void()>>> v_attr;

    /// attribute storage
    std::vector<data_type> v_data;
  };
}

#include "Group.cc"

// so that many of the tuple shenanigans work
namespace std {
  template <typename ...Ts>
  struct tuple_size<Framework::Group<Ts...>> : std::integral_constant<size_t, sizeof...(Ts)> {};

  template <size_t I, typename ...Ts>
  struct tuple_element<I, Framework::Group<Ts...>> : tuple_element<I, typename Framework::Group<Ts...>::sorted> {};
}

#endif
