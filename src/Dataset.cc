// -*- C++ -*-
// author: afiq anuar
// short: please refer to header for information

template <typename Tree>
Framework::Dataset<Tree>::Dataset(const std::string &name_, const std::string &tree_name_, 
                                  const std::string &tree_struct_ /*= ""*/, char tree_delim_ /*= ' '*/,
                                  const std::vector<std::string> &v_file_ /*= {}*/) :
  name(name_),
  tree_name(tree_name_),
  tree_struct(tree_struct_),
  tree_delim(tree_delim_),
  v_file(v_file_),
  v_entry({}),
  evaluated(false),
  tree_ptr(nullptr),
  allocator(Allocator{}),
  analyzer(nullptr),
  v_weight({})
{ gErrorIgnoreLevel = kBreak; }



template <typename Tree>
void Framework::Dataset<Tree>::set_tree(const std::string &tree_name_)
{
  tree_name = tree_name_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_structure(const std::string &tree_struct_)
{
  tree_struct = tree_struct_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_delimiter(char tree_delim_)
{
  tree_delim = tree_delim_;
}



template <typename Tree>
void Framework::Dataset<Tree>::set_files(const std::vector<std::string> &v_file_, int nfile /*= -1*/, bool force_replace /*= false*/)
{
  if (v_file_.empty())
    return;

  if (force_replace)
    reset();

  if (not v_file.empty())
    return;

  nfile = (nfile > 0 and nfile <= v_file_.size()) ? nfile : v_file_.size();
  if (nfile != v_file_.size()) {
    for (int ifile = 0; ifile < nfile; ++ifile)
      v_file.emplace_back(v_file_[ifile]);
  }
  else 
    v_file = v_file_;
}



template <typename Tree>
bool Framework::Dataset<Tree>::add_file(const std::string &file)
{
  auto iF = std::find(std::begin(v_file), std::end(v_file), file);
  if (iF != std::end(v_file))
    return false;

  v_file.emplace_back(file);
  return true;
}



template <typename Tree>
int Framework::Dataset<Tree>::n_files() const
{
  return v_file.size();
}



template <typename Tree>
bool Framework::Dataset<Tree>::add_weight(const std::string &wgt_name, double wgt)
{
  auto iW = std::find_if(std::begin(v_weight), std::end(v_weight), [&wgt_name] (const auto &weight) {return weight.first == wgt_name;});
  if (iW != std::end(v_weight))
    return false;

  v_weight.emplace_back(wgt_name, wgt);
  return true;
}



template <typename Tree>
Framework::Dataset<Tree> Framework::Dataset<Tree>::split(int npartition /*= 2*/, int ipartition /*= 0*/)
{
  if (npartition < 0 or npartition > v_file.size())
    npartition = 2;

  if (ipartition < 0 or ipartition >= npartition)
    ipartition = 0;

  const int nfile = (v_file.size() % npartition) ? int(v_file.size() / npartition) + 1 : v_file.size() / npartition;
  std::vector<std::string> vf;

  for (int iF = ipartition * nfile, iC = 0; iF < v_file.size() and iC < nfile; ++iF, ++iC)
    vf.emplace_back(v_file[iF]);

  auto dat = Dataset(name + "_n" + to_str(npartition) + "_i" + to_str(ipartition), tree_name, tree_struct, tree_delim, vf);
  dat.v_weight = v_weight;

  return dat;
}



template <typename Tree>
long long Framework::Dataset<Tree>::current_entry(long long entry) const
{
  if (tree_ptr == nullptr)
    return -1;

  return tree_ptr->LoadTree(entry);
}



template <typename Tree>
const std::unique_ptr<Tree>& Framework::Dataset<Tree>::tree(int ispecify /*= -1*/) const
{
  static int ifile = 0;
  ifile = (ispecify > -1 and ispecify < n_files()) ? ispecify : ifile;
  if (ifile >= n_files())
    throw std::runtime_error( "ERROR: Dataset::tree already loaded the last available tree!! Ensure that all requested branches exist in at least one file in the Dataset!!" );
  current_entry(v_entry[ifile++]);
  return tree_ptr;
}



template <typename Tree>
double Framework::Dataset<Tree>::get_weight(const std::string &wgt_name) const
{
  auto iW = std::find_if(std::begin(v_weight), std::end(v_weight), [&wgt_name] (const auto &weight) {return weight.first == wgt_name;});
  if (iW == std::end(v_weight)) {
    // TODO something about logger telling absence of weight
    return 0.;
  }
  else 
    return iW->second;
}



template <typename Tree>
template <typename ...Collections>
void Framework::Dataset<Tree>::associate(Collections &...colls)
{
  static_assert(sizeof...(colls) > 0, "ERROR: Dataset::allocate makes no sense when called without arguments!!");
  evaluate();

  if (tree_ptr == nullptr)
    throw std::runtime_error( "ERROR: Dataset::associate should not be called before assigning the files to be analyzed!!" );

  // apparently the tree is loaded only after this call
  // so associate will fail without it
  tree_ptr->GetEntries();

  (colls.associate(*this), ...);
  allocator.set_allocator([&colls...] () { (colls.reassociate(), ...); });
  tree_ptr->SetNotify(&allocator);
}



template <typename Tree>
template <typename Analyzer>
void Framework::Dataset<Tree>::set_analyzer(Analyzer analyzer_, bool force_replace /*= false*/)
{
  using Traits = function_traits<decltype(analyzer_)>;
  static_assert(Traits::arity == 1 and std::is_convertible_v<typename Traits::template bare_arg<0>, long long>, 
                "ERROR: Dataset::set_analyzer only takes one argument that is convertible to entry number!!");
  static_assert(std::is_same_v<typename Traits::result_type, void>, 
                "ERROR: Dataset::set_analyzer: currently non-void return type is not supported!!");

  if (force_replace or not analyzer)
    analyzer = std::function<void(long long)>(analyzer_);
}



template <typename Tree>
void Framework::Dataset<Tree>::analyze(long long total /*= -1LL*/, long long skip /*= -1LL*/, bool silent /*=false*/) const
{
  if (tree_ptr == nullptr)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before assigning the files to be analyzed!!" );

  if (not allocator)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before calling Dataset::associate!!" );

  if (not analyzer)
    throw std::runtime_error( "ERROR: Dataset::analyze should not be called before calling Dataset::set_analyzer!!" );

  const auto dEvt = (total > 0LL and total <= tree_ptr->GetEntries()) ? total : tree_ptr->GetEntries();
  if (not silent)
    std::cout << "Processing " << dEvt << " events..." << std::endl;

  const bool doskip = skip > 0LL and skip < total;
  if (not silent and doskip)
    std::cout << "Skipping " << skip << " events..." << std::endl;

  for (auto cEvt = (doskip) ? skip : 0LL; cEvt < dEvt; ++cEvt)
    analyzer(current_entry(cEvt));
  if (not silent)
    std::cout << "Processed " << dEvt << " events!" << std::endl;
}



template <typename Tree>
void Framework::Dataset<Tree>::reset()
{
  v_file.clear();
  v_file.shrink_to_fit();

  v_entry.clear();
  v_entry.shrink_to_fit();

  evaluated = false;

  if (tree_ptr != nullptr) {
    tree_ptr->ResetBranchAddresses();
    tree_ptr->Reset();
    tree_ptr.reset();
  }

  allocator = Allocator{};
  analyzer = nullptr;

  v_weight.clear();
  v_weight.shrink_to_fit();
}



template <>
void Framework::Dataset<TChain>::evaluate()
{
  if (v_file.empty())
    return;

  if (evaluated)
    return;

  if (tree_ptr == nullptr) {
    tree_ptr = std::make_unique<TChain>(tree_name.c_str());
    //tree_ptr->SetImplicitMT(false);
    tree_ptr->SetBranchStatus("*", 0);
  }

  for (const auto &file : v_file) {
    if (std::filesystem::is_regular_file(file)) {
      tree_ptr->Add(file.c_str());
      v_entry.emplace_back(tree_ptr->GetEntries() - 1LL);
    }
    else
      throw std::runtime_error( "ERROR: Dataset::evaluate: file " + file + " is not accessible by the program!!" );
  }

  evaluated = true;
}



template <>
void Framework::Dataset<TTree>::evaluate()
{
  static bool flag_struct = false;
  if (not flag_struct and tree_struct == "") {
    // TODO something about logging the error
    return;
  }

  if (v_file.empty())
    return;

  if (evaluated)
    return;

  if (tree_ptr == nullptr) {
    tree_ptr = std::make_unique<TTree>(tree_name.c_str(), "");
    //tree_ptr->SetImplicitMT(false);
    tree_ptr->SetBranchStatus("*", 0);
  }

  for (const auto &file : v_file) {
    if (not std::filesystem::is_regular_file(file))
      throw std::runtime_error( "ERROR: Dataset::evaluate: file " + file + " is not accessible by the program!!" );

    if (not flag_struct) {
      tree_ptr->ReadFile(file.c_str(), tree_struct.c_str(), tree_delim);
      flag_struct = true;
    }
    else
      tree_ptr->ReadFile(file.c_str());
  }

  evaluated = true;
}
