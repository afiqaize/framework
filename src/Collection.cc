// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename ...Ts>
Framework::Collection<Ts...>::Collection(const std::string &name_, int reserve_) : 
Framework::Group<Ts...>::Group(name_, 0),
tree(nullptr),
counter_name(""),
counter_branch(nullptr)
{
  reserve(reserve_);
  this->initialize(1);
}



template <typename ...Ts>
Framework::Collection<Ts...>::Collection(const std::string &name_, const std::string &counter_name_, int reserve_, int init /*= 4*/) : 
Framework::Group<Ts...>::Group(name_, 0),
tree(nullptr),
counter_name(counter_name_),
counter_branch(nullptr)
{
  reserve(reserve_);
  if (counter_name != "")
    this->initialize(init);
  else
    this->initialize(1);
}



template <typename ...Ts>
Framework::Collection<Ts...>::~Collection()
{
  detach();
}



template <typename ...Ts>
void Framework::Collection<Ts...>::reserve(int attr)
{
  v_branch.reserve(attr);
  Group<Ts...>::reserve(attr);
}



template <typename ...Ts>
bool Framework::Collection<Ts...>::add_attribute(const std::string &attr, const std::string &branch)
{
  if (this->has_attribute(attr))
    return false;

  this->v_attr.emplace_back(attr, nullptr);
  v_branch.emplace_back(branch, nullptr);
  this->v_data.emplace_back(std::vector<typename std::tuple_element_t<0, std::tuple<Ts...>>>());

  std::visit([init = this->v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, this->v_data.back());
  return true;
}



template <typename ...Ts>
template <typename Number>
bool Framework::Collection<Ts...>::add_attribute(const std::string &attr, const std::string &branch, Number)
{
  static_assert(contained_in<Number, Ts...> or (std::is_same_v<Number, bool> and contained_in<boolean, Ts...>), 
                "ERROR: Collection::add_attribute: the initializing Number type is not among the types expected by the Collection!!");

  using Attribute = typename std::conditional<std::is_same_v<Number, bool>, boolean, Number>::type;

  if (add_attribute(attr, branch)) {
    Group<Ts...>::template retype<Attribute>(this->v_data.back());
    return true;
  }

  return false;
}



template <typename ...Ts>
template <typename Function, typename ...Attributes>
bool Framework::Collection<Ts...>::transform_attribute(const std::string &attr, Function function, Attributes &&...attrs)
{
  if (Group<Ts...>::transform_attribute(attr, function, std::forward<Attributes>(attrs)...)) {
    v_branch.emplace_back("", nullptr);
    std::visit([init = this->v_index.capacity()] (auto &vec) {vec.reserve(init); vec.clear();}, this->v_data.back());
    return true;
  }

  return false;
}



template <typename ...Ts>
template <typename Tree>
void Framework::Collection<Ts...>::associate(Dataset<Tree> &dataset)
{
  tree = dataset.tree().get();

  if (counter_name != "") {
    tree->SetBranchStatus(counter_name.c_str(), 1);
    tree->SetBranchAddress(counter_name.c_str(), &(this->counter), &counter_branch);
    counter_branch->SetAutoDelete(false);
  }

  for (int iB = 0; iB < v_branch.size(); ++iB) {
    auto &[branch_name, branch] = v_branch[iB];
    if (branch_name == "")
      continue;

    // experimenting with automatic type readout
    auto leaf = tree->GetLeaf(branch_name.c_str());
    if (leaf == nullptr)
      leaf = tree->FindLeaf(branch_name.c_str());
    if (leaf == nullptr)
      throw std::runtime_error( std::string("ERROR: Collection::associate: unable to find the requested branch ") + branch_name + "!!!" );

    std::string btype = leaf->GetTypeName();
    if (btype == "Int_t")
      Group<Ts...>::template retype<int>(this->v_data[iB]);
    else if (btype == "UInt_t")
      Group<Ts...>::template retype<uint>(this->v_data[iB]);
    else if (btype == "Float_t")
      Group<Ts...>::template retype<float>(this->v_data[iB]);
    else if (btype == "Double_t")
      Group<Ts...>::template retype<double>(this->v_data[iB]);
    else if (btype == "Long64_t")
      Group<Ts...>::template retype<long long>(this->v_data[iB]);
    else if (btype == "ULong64_t")
      Group<Ts...>::template retype<unsigned long long>(this->v_data[iB]);
    else if (btype == "Bool_t")
      Group<Ts...>::template retype<boolean>(this->v_data[iB]);
    else if (btype == "Char_t")
      Group<Ts...>::template retype<char>(this->v_data[iB]);
    else if (btype == "UChar_t")
      Group<Ts...>::template retype<unsigned char>(this->v_data[iB]);
    else
      throw std::runtime_error( std::string("ERROR: Collection::associate: branch ") + branch_name + "has an unsupported type " + btype
                                + "!!! If it should have been supported, please add it and/or contact the developer.");

    tree->SetBranchStatus(branch_name.c_str(), 1);
    std::visit([this, &branch = branch, &branch_name = branch_name] (auto &vec) { 
        tree->SetBranchAddress(branch_name.c_str(), vec.data(), &branch);
      }, this->v_data[iB]);
    branch->SetAutoDelete(false);
  }
}



template <typename ...Ts>
void Framework::Collection<Ts...>::reassociate()
{
  if (tree == nullptr)
    throw std::runtime_error( "ERROR: Collection::reassociate: the associated tree is null." 
                              "Perhaps Collection::associate has not been called? Aborting!!" );

  if (counter_name == "")
    return;

  TLeaf *leaf = tree->GetLeaf(counter_name.c_str());
  auto current_max = leaf->GetMaximum();
  if (current_max > this->v_index.capacity()) {
    this->initialize(current_max);

    for (int iD = 0; iD < this->v_data.size(); ++iD) {
      auto &[branch_name, branch] = v_branch[iD];
      if (branch == nullptr)
        continue;

      std::visit([this, &branch = branch, &branch_name = branch_name] (auto &vec) {
          tree->SetBranchAddress(branch_name.c_str(), vec.data(), &branch);
        }, this->v_data[iD]);
    }
  }
}



template <typename ...Ts>
void Framework::Collection<Ts...>::populate(long long entry)
{
  // get the number of elements and fill up indices
  if (counter_branch != nullptr) {
    counter_branch->GetEntry(entry);
    this->selected = this->counter;

    this->v_index.clear();
    for (int iD = 0; iD < this->counter; ++iD)
      this->v_index.emplace_back(iD);
  }
  else if (!this->v_attr.empty()) {
    this->v_index.clear();
    this->v_index.emplace_back(0);
    this->counter = 1;
    this->selected = this->counter;
  }

  // and then get the data of all the branches
  for (int iD = 0; iD < this->v_data.size(); ++iD) {
    if (v_branch[iD].second == nullptr)
      continue;

    v_branch[iD].second->GetEntry(entry);
  }

  // functional transformations can only run after everything else is populated
  for (int iD = 0; iD < this->v_data.size(); ++iD) {
    if (v_branch[iD].second == nullptr and this->v_attr[iD].second)
      this->v_attr[iD].second();
  }
}



template <typename ...Ts>
void Framework::Collection<Ts...>::detach()
{
  for (int iD = 0; iD < this->v_data.size(); ++iD) {
    if (v_branch[iD].second != nullptr)
      v_branch[iD].second->ResetAddress();
  }
}
