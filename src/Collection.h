// -*- C++ -*-
// author: afiq anuar
// short: collections are groups whose input come from ROOT files - i.e. it is the group-dataset interface

#ifndef FWK_COLLECTION_H
#define FWK_COLLECTION_H

#include "Group.h"
#include "Dataset.h"
#include "TBranch.h"
#include "TLeaf.h"

namespace Framework {
  template <typename ...Ts>
  class Collection : public Group<Ts...> {

  public:
    /// constructor
    Collection() = delete;

    /// this is for single element collections
    /// reserve_ is an estimate of the number of attributes the collection will hold
    Collection(const std::string &name_, int reserve_);

    /// this is for array collections
    /// counter_branch is the presumed branch in the dataset that specifies element count
    /// reserve_ is an estimate of the number of attributes the collection will hold
    /// init is an estimate for the number of elements the collection will hold
    /// 4 is just a "reasonable default"; the code reallocates if capacity is hit
    /// that said, reallocating too often is not good for performance 
    Collection(const std::string &name_, const std::string &counter_name_, int reserve_, int init = 4);

    Collection(const Collection &coll) = delete;
    Collection(Collection &&coll) = default;

    /// destructor
    ~Collection() { detach(); }

    /// assignment
    Collection& operator=(const Collection &dat) = default;
    Collection& operator=(Collection &&dat) = default;

    /// reserve the space for expected number of attributes
    void reserve(int attr);

    /// add an attribute into the collection
    /// returns false upon failure to add the attribute
    /// this can happen if the data type is inconsistent with the collection
    /// returns true upon a successful registration
    bool add_attribute(const std::string &attr, const std::string &branch);

    /// overload for cases when manually specifying the branch type is desirable
    /// will be overridden by associate or transform_attribute if a wrong type is given
    /// the last argument is also used to specify a default value to be used when the branch is not present
    /// necessary when e.g. the default is not the same as default-constructed value of the first type
    template <typename Number>
    bool add_attribute(const std::string &attr, const std::string &branch, Number number);

    /// transform a group of internal attributes into another attribute
    /// the transformation is done element-wise on every element of held data
    /// as this adds a new attribute, its name has to be unique
    template <typename Function, typename ...Attributes>
    bool transform_attribute(const std::string &attr, Function function, Attributes &&...attrs);

    /// associate the attributes to relevant branches in a Dataset
    template <typename Tree>
    void associate(Dataset<Tree> &dataset);

    /// reassociate the attributes to relevant branches in a Dataset
    /// ie to be ran when a reallocation is needed due to capacity hit
    void reassociate();

    /// populate the data with information read from the branches
    void populate(long long entry) override;

  protected:
    /// detach the branches
    void detach();

    /// TBranch isn't enough for realloc, so Collection needs to remember its tree
    TTree *tree;

    /// branch name which gets the counter
    std::string counter_name;

    /// branch associated to the counter
    TBranch *counter_branch;

    /// attribute branches
    std::vector<std::pair<std::string, TBranch *>> v_branch;

    /// default values in case the branch is not present in the current tree
    std::vector<std::variant<Ts...>> v_default;
  };
}

#include "Collection.cc"

// so that many of the tuple shenanigans work
namespace std {
  template <typename ...Ts>
  struct tuple_size<Framework::Collection<Ts...>> : tuple_size<Framework::Group<Ts...>> {};

  template <size_t I, typename ...Ts>
  struct tuple_element<I, Framework::Collection<Ts...>> : tuple_element<I, typename Framework::Group<Ts...>> {};
}

#endif
