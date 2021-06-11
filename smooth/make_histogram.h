// -*- C++ -*-
// author: afiq anuar
// short: provides function definitions needed to make histograms from flat trees

#ifndef FWK_MAKE_HISTOGRAM_H
#define FWK_MAKE_HISTOGRAM_H

#include "Dataset.h"
#include "Collection.h"

#include <cmath>
#include "misc/string_io.h"
#include "misc/container_util.h"
#include "misc/function_util.h"
#include "misc/rng_util.h"

#include "array_histogram.h"
#include "output_util.h"

bool valid_name(const std::string &str)
{
  const auto is_valid = [] (unsigned char c) { return std::isalnum(c) or c == '_'; };
  return std::all_of(std::begin(str), std::end(str), is_valid);
}



/// is a string elementary or not; if yes, returns a vector of string containing its name, the operation involved and the operands
/// op1(a) and a op2 b are both elementary, anything else are not
/// op1/op2 can be any of the supported binary/unary operations
/// and a, b are valid names i.e. containing only alphanumeric + underscores
std::vector<std::string> is_elementary(const std::string &expression, const std::vector<std::string> &unaries, const std::vector<std::string> &binaries)
{
  const auto iop = expression.find("("), icl = expression.find(")");

  if (iop == std::string::npos and icl == std::string::npos) {
    for (const auto &binary : binaries) {
      auto ebi = split(expression, binary);
      if (ebi.size() == 2) {

        if (valid_name(ebi[0]) and valid_name(ebi[1]))
          return {random_variable_name(3), binary, ebi[0], ebi[1]};
      }
    }
  }
  else if (iop != std::string::npos and icl == expression.size() - 1) {
    for (const auto &unary : unaries) {
      if (count_substring(expression, unary) == 1) {
        std::string eun = expression.substr(iop + 1, icl - iop - 1);

        if (valid_name(eun) and expression.substr(0, unary.size()) == unary)
          return {random_variable_name(3), unary, eun};
      }
    }
  }

  return {};
}



void unique_emplace(std::vector<std::vector<std::string>> &expressions, std::vector<std::string> &expression, 
                    std::string &full, const std::string &toreplace)
{
  if (expression.empty())
    return;

  // first check if we already have the same elementary expression in the list
  auto same_all_but_first = [&expression] (auto &var) {
    bool all_equal = (var.size() == expression.size());

    if (all_equal) {
      for (int iv = 1; iv < var.size(); ++iv) {
        all_equal = all_equal and (var[iv] == expression[iv]);
      }
    }

    return all_equal;
  };

  auto isame = std::find_if(std::begin(expressions), std::end(expressions), same_all_but_first);
  if (isame != std::end(expressions)) {
    while( replace(full, toreplace, (*isame)[0]) )
      ;
    return;
  }

  // ok no, go ahead and dump this in after ensuring its name is unique
  while (std::find_if(std::begin(expressions), std::end(expressions), [&expression] (auto &var) { return var[0] == expression[0]; }) != std::end(expressions))
    expression[0] = random_variable_name(3);

  while( replace(full, toreplace, expression[0]) )
    ;
  expressions.emplace_back(std::move(expression));
}



/// parse the expression list into constituent function and arguments
/// given that each expresions is of the form "c", "a op b" (read: c = a + b) or "c", "op(a)", where a and b are branches
/// the parser transforms this into "c", "op", "a", "b", which are then used by transform_attribute
/// in cases where the expression is a more complicated combination of the above, it first reduces it to the above
/// with the branch name automatically generated
void parse_expressions(std::vector<std::vector<std::string>> &expressions)
{
  static const std::vector<std::string> unaries = {"exp(", "log(", "log10(",
                                                   "sin(", "cos(", "tan(", "asin(", "acos(", "atan(",
                                                   "sqrt(", "abs(", "negate(", "invert("};
  static const std::vector<std::string> binaries = {"*", "/", "+", "-"};

  // decompose expressions like a + b + c + d + e -> f + c + d + e -> g + d + e -> h + e...
  // the refs are what needs to be kept track of for the inplace edits within loop body
  // iop and icl are opening and closing indices of the part of string to analyze 
  auto decompose_binaries = [] (std::vector<std::vector<std::string>> &expressions, std::string &exp, std::string &part,
                                std::vector<std::string> &elem, std::string::size_type iop, std::string::size_type icl) {
    auto imd = min(exp.find("*", iop), exp.find("/", iop));
    if (imd > icl)
      imd = std::string::npos;

    if (imd != std::string::npos) {
      auto ipm1 = min(exp.find("+", iop), exp.find("-", iop));
      auto ibi = min(exp.find("+", imd + 1), exp.find("-", imd + 1), exp.find("*", imd + 1), exp.find("/", imd + 1));

      if (ipm1 > imd)
        part = exp.substr(iop, (ibi < icl) ? ibi - iop : icl - iop);
      else {
        auto ipm2 = min(exp.find("+", ipm1 + 1), exp.find("-", ipm1 + 1));
        while (ipm2 < imd) {
          ipm1 = ipm2;
          ipm2 = min(exp.find("+", ipm1 + 1), exp.find("-", ipm1 + 1));
        }

        part = exp.substr(ipm1 + 1, (ibi < icl) ? ibi - ipm1 - 1 : icl - ipm1 - 1);
      }

      auto belem = is_elementary(part, unaries, binaries);
      if (!belem.empty()) {
        unique_emplace(expressions, belem, exp, part);
        elem = is_elementary(exp, unaries, binaries);
        return true;
      }
    }
    else {
      auto ipm1 = min(exp.find("+", iop), exp.find("-", iop));
      if (ipm1 > icl)
        imd = std::string::npos;

      if (ipm1 != std::string::npos) {
        auto ipm2 = min(exp.find("+", ipm1 + 1), exp.find("-", ipm1 + 1));
        part = exp.substr(iop, (ipm2 < icl) ? ipm2 - iop : icl - iop);

        auto belem = is_elementary(part, unaries, binaries);
        if (!belem.empty()) {
          unique_emplace(expressions, belem, exp, part);
          elem = is_elementary(exp, unaries, binaries);
          return true;
        }
      }
    }

    return false;
  };

  const int nvar = expressions.size();
  for (int ivar = 0; ivar < nvar; ++ivar) {
    auto &exp = expressions[ivar][1];
    if (exp == "")
      continue;

    // get all the 'raw' attribute names and dump them in
    // need to first get rid of the unaries
    auto cexp = exp;
    for (const auto &unary : unaries)
      strip(cexp, unary);

    auto branches = split_if(cexp, [] (unsigned char c) { return std::isalnum(c) or c == '_'; });
    for (auto &branch : branches) {
      if (std::find_if(std::begin(expressions), std::end(expressions), [&branch] (auto &var) { return var[0] == branch; }) == std::end(expressions))
        expressions.emplace_back(std::vector<std::string>{branch, ""});
    }

    auto elem = is_elementary(exp, unaries, binaries);
    std::string part = "";
    while (elem.empty()) {
      // a flag to indicate that an elementary expression is found and dumped into expression list
      // due to the inplace editing, we want to be analyzing it in a fresh iteration
      bool dumped = false;

      // find innermost opening and closing brackets
      auto nop = count_substring(exp, "("), ncl = count_substring(exp, ")");
      if (nop != ncl) {
        std::cerr << "Expression " << exp << " is invalid due to non-matching parenthesis. skipping..." << std::endl;
        return;
      }

      if (nop > 0) {
        auto iop1 = exp.find("("), icl = exp.find(")");
        auto iop2 = exp.find("(", iop1 + 1);
        while (iop2 < icl) {
          iop1 = iop2;
          iop2 = exp.find("(", iop1 + 1);
        }

        // check whether it is an elementary unary
        for (const auto &unary : unaries) {
          int iun = iop1 - unary.size() + 1;

          if (exp.substr((iun < 0) ? 0 : iun, unary.size()) == unary) {
            part = exp.substr((iun < 0) ? 0 : iun, icl - iop1 + unary.size());
            auto uelem = is_elementary(part, unaries, binaries);

            if (!uelem.empty()) {
              unique_emplace(expressions, uelem, exp, part);
              dumped = true;
              break;
            }
          }
        }

        if (dumped) {
          elem = is_elementary(exp, unaries, binaries);
          continue;
        }

        // nope, carry on analyzing the insides for binaries
        if (!decompose_binaries(expressions, exp, part, elem, iop1 + 1, icl)) {
          // neither elementary unary or binary; can only be something like (variable)
          // just delete the brackets and move on
          exp.erase(icl, 1);
          exp.erase(iop1, 1);
          elem = is_elementary(exp, unaries, binaries);
        }
      }
      else
        decompose_binaries(expressions, exp, part, elem, 0, std::string::npos);
    }

    exp = elem[1];
    for (int ielem = 2; ielem < elem.size(); ++ielem)
      expressions[ivar].emplace_back(elem[ielem]);
  }
}



/// arguments: vector of strings as expected by --variable and --weight CL options, and input file names for source branch promotion to double
/// returns a tuple containing three vectors
/// first contains variables and (parsed) expressions information
/// second is the coarse binning information i.e. exactly as requested
/// third is fine binning, i.e. coarse split into n bins along each dimension, n being dimension-specific
/// in any case where anything is invalid, returned output is blank
/// segfaults at dataset.analyze() if the source promotion is done within a make_collection for yet to be understood reasons
auto variables_and_binning(const std::vector<std::string> &variables, const std::string &weight, const std::vector<std::string> &files, const std::string &tree)
{
  using namespace Framework;

  static const std::vector<int> fsplits = {12, 5, 3, 2};
  auto fsplit = (variables.size() > 2) ? fsplits.back() : fsplits[variables.size() - 1];

  std::vector<std::vector<std::string>> expressions;
  std::vector<std::vector<double>> cedges, fedges;

  for (const auto &var : variables) {
    auto veb = split(var, ":");
    if (veb.size() != 2) {
      std::cerr << "Variable or expression " << var << " is invalid" << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{}, 
                             std::string(""));
    }

    auto ve = split(veb[0], "=");
    if (ve.size() != 1 and ve.size() != 2) {
      std::cerr << "Variable or expression " << veb[0] << " is invalid" << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{},
                             std::string(""));
    }
    strip(ve[0]);
    if (!valid_name(ve[0])) {
      std::cerr << "Invalid variable name " << ve[0] << ". Aborting. Current version considers only names containing " 
        "alphanumeric characters or underscores"<< std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{},
                             std::string(""));
    }

    if (ve.size() == 2)
      strip(ve[1]);
    else
      ve.emplace_back("");
    expressions.emplace_back(std::move(ve));

    auto binstr = split(veb[1], ";");
    int nbinf = std::stoi(binstr[0]);
    double minf = (binstr.size() == 3) ? std::stod(binstr[1]) : 0., maxf = (binstr.size() == 3) ? std::stod(binstr[2]) : 0.;
    const auto cedge = [&binstr] {
      std::vector<double> cedge;
      for (const auto &edge : binstr)
        cedge.emplace_back(std::stod(edge));
      return cedge;
    }();

    if (nbinf > 1 and minf < maxf) {
      cedges.emplace_back( make_interval(minf, maxf, (maxf - minf) / nbinf) );
      fedges.emplace_back( make_interval(minf, maxf, (maxf - minf) / (fsplit * nbinf)) );
    }
    else if (binstr.size() > 2 and std::is_sorted(std::begin(cedge), std::end(cedge))) {
      fedges.emplace_back( [&cedge, fsplit] {
        std::vector<double> fedge = {cedge[0]};
        for (int ic = 0; ic < cedge.size() - 1; ++ic)
          fill_interval(fedge, cedge[ic + 1], (cedge[ic + 1] - cedge[ic]) / fsplit);

        return fedge;
      }() );
      cedges.emplace_back( std::move(cedge) );
    }
    else {
      std::cerr << "Binning information for variable  " << ve[0] << " is invalid. Note that variables with only one bin is not accepted" << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{},
                             std::string(""));
    }
  }

  auto we = split(weight, "=");
  strip(we[0]);
  if (weight != "") {
    if (valid_name(we[0])) {
      if (we.size() == 2)
        strip(we[1]);
      else
        we.emplace_back("");
      expressions.emplace_back(we);
    }
    else {
      std::cerr << "Invalid weight name " << we[0] << ". Ignoring the weight. Current version considers only names containing " 
        "alphanumeric characters or underscores" << std::endl;
    }
  }

  /**** start testing of expression parser ****
  for (auto &exp : expressions) {
    for (auto &var : exp)
      std::cout << var << " -- ";
    std::cout << std::endl;
  }
  std::cout << "     ----------------     " << std::endl;
  */
  parse_expressions(expressions);
  /*
  for (auto &exp : expressions) {
    for (auto &var : exp)
      std::cout << var << " -- ";
    std::cout << std::endl;
  }
  std::cout << "     ----------------     " << std::endl;
  ****  end testing of expression parser  ****/

  Dataset<TChain> dataset("dataset", tree);
  dataset.set_files(files);

  // tag source branches i.e. those coming from source root file, which may be float
  std::vector<std::string> branches;
  Collection<float, double> src("src", expressions.size());
  for (const auto &var : expressions) {
    if (var.size() == 2 and var[1] == "") {
      src.add_attribute(var[0], var[0]);
      branches.emplace_back(var[0]);
    }
  }
  dataset.associate(src);

  // promote source attributes to double if float
  for (const auto &branch : branches) {
    const bool isfloat = src.get_if<float>(branch) != nullptr;

    if (isfloat) {
      auto name = random_variable_name(3);
      while (std::find_if(std::begin(expressions), std::end(expressions), [&name] (auto &var) { return var[0] == name; }) != std::end(expressions))
        name = random_variable_name(3);

      auto ib = std::find_if(std::begin(expressions), std::end(expressions), [&branch] (auto &var) { return var[0] == branch; });
      (*ib)[0] = name;
      (*ib)[1] = "__source__(";
      ib->emplace_back(branch);

      expressions.emplace_back(std::vector<std::string>{branch, "__promote__(", name});
    }
  }

  return std::make_tuple(std::move(expressions), std::move(cedges), std::move(fedges), we[0]);
}



/// prepare the collection to read the variables
auto make_collection(Framework::Dataset<TChain> &dataset,
                     const std::tuple<
                     std::vector<std::vector<std::string>>,
                     std::vector<std::vector<double>>,
                     std::vector<std::vector<double>>,
                     std::string> &variables_bins)
{
  using namespace Framework;

  const auto &variables = std::get<0>(variables_bins);
  if (variables.empty())
    throw std::runtime_error( "make_collection() :: an error is encountered. aborting." );

  Collection<float, double> coll("coll", variables.size());
  for (const auto &var : variables) {
    if (var.size() == 2 and var[1] == "")
      coll.add_attribute(var[0], var[0], 1.);
    else if (var.size() > 2 and var[1] == "__source__(")
      coll.add_attribute(var[0], var[2], 1.f);
  }
  dataset.associate(coll);

  while (coll.n_attributes() != variables.size()) {
    for (const auto &var : variables) {
      if (var.size() > 2) {
        bool isthere = true;
        for (int iarg = 2; iarg < var.size(); ++iarg)
          isthere = isthere and coll.has_attribute(var[iarg]);

        if (isthere) {
          if (var[1] == "exp(")
            coll.transform_attribute(var[0], [] (double x) {return std::exp(x);}, var[2]);
          else if (var[1] == "log(")
            coll.transform_attribute(var[0], [] (double x) {return std::log(x);}, var[2]);
          else if (var[1] == "log10(")
            coll.transform_attribute(var[0], [] (double x) {return std::log10(x);}, var[2]);
          else if (var[1] == "sin(")
            coll.transform_attribute(var[0], [] (double x) {return std::sin(x);}, var[2]);
          else if (var[1] == "cos(")
            coll.transform_attribute(var[0], [] (double x) {return std::cos(x);}, var[2]);
          else if (var[1] == "tan(")
            coll.transform_attribute(var[0], [] (double x) {return std::tan(x);}, var[2]);
          else if (var[1] == "asin(")
            coll.transform_attribute(var[0], [] (double x) {return std::asin(x);}, var[2]);
          else if (var[1] == "acos(")
            coll.transform_attribute(var[0], [] (double x) {return std::acos(x);}, var[2]);
          else if (var[1] == "atan(")
            coll.transform_attribute(var[0], [] (double x) {return std::atan(x);}, var[2]);
          else if (var[1] == "sqrt(")
            coll.transform_attribute(var[0], [] (double x) {return std::sqrt(x);}, var[2]);
          else if (var[1] == "abs(")
            coll.transform_attribute(var[0], [] (double x) {return std::abs(x);}, var[2]);
          else if (var[1] == "negate(")
            coll.transform_attribute(var[0], [] (double x) {return -x;}, var[2]);
          else if (var[1] == "invert(")
            coll.transform_attribute(var[0], [] (double x) {return 1. / x;}, var[2]);
          else if (var[1] == "__promote__(")
            coll.transform_attribute(var[0], [] (float x) {return double(x);}, var[2]);
          else if (var[1] == "+")
            coll.transform_attribute(var[0], [] (double x, double y) {return x + y;}, var[2], var[3]);
          else if (var[1] == "-")
            coll.transform_attribute(var[0], [] (double x, double y) {return x - y;}, var[2], var[3]);
          else if (var[1] == "*")
            coll.transform_attribute(var[0], [] (double x, double y) {return x * y;}, var[2], var[3]);
          else if (var[1] == "/")
            coll.transform_attribute(var[0], [] (double x, double y) {return x / y;}, var[2], var[3]);
        }
      }
    }
  }

  return coll;
}



/// returns a set of histograms constructed from the provided file
/// the set contains npartition histograms, each containing roughly 1/N of events from the dataset
/// histogram is represented as a vector of bins, which is in turn represented by 2 doubles: sum of weight and weight squared respectively
/// under- and overflows are automatically added to the first and last bins along the relevant axis respectively
/// iedge denote which binning to be used, 1 for the first binning held by variables_bins, 2 is second
std::vector<Arrayhist> count_and_bin(Framework::Dataset<TChain> &dataset,
                                     Framework::Collection<float, double> &coll,
                                     const std::tuple<
                                     std::vector<std::vector<std::string>>,
                                     std::vector<std::vector<double>>,
                                     std::vector<std::vector<double>>,
                                     std::string> &variables_bins,
                                     int npartition = 1,
                                     int iedge = 1,
                                     bool weighted = true,
                                     bool silent = false)
{
  using namespace Framework;

  const auto &variables = std::get<0>(variables_bins);
  const auto &bins = (iedge == 1) ? std::get<1>(variables_bins) : std::get<2>(variables_bins);
  const auto &weight = std::get<3>(variables_bins);

  dataset.associate(coll);

  const int nvar = bins.size();

  const auto ivars = [&variables, &nvar, &coll] () {
    std::vector<int> ivars(nvar);
    for (int ivar = 0; ivar < nvar; ++ivar)
      ivars[ivar] = coll.inquire(variables[ivar][0]);

    return ivars;
  }();

  const int nbin = count_nbin(bins), iweight = coll.inquire(weight);
  std::vector<Arrayhist> hists;
  hists.reserve(npartition);
  for (int ipart = 0; ipart < npartition; ++ipart)
    hists.emplace_back(Arrayhist(nbin));

  thread_local auto rng = random_generator<>();
  auto dist = std::uniform_int_distribution{0, npartition - 1};
  std::vector<double> values(nvar);

  auto fwgt = [npartition, &dist] (std::vector<Arrayhist> &hists, int ibin, const auto &wgt) {
    int ipart = (npartition == 1) ? 0 : dist(rng);
    hists[ipart](ibin, 0) += wgt;
    hists[ipart](ibin, 1) += wgt * wgt;
  };

  auto fone = [npartition, &dist] (std::vector<Arrayhist> &hists, int ibin) {
    int ipart = (npartition == 1) ? 0 : dist(rng);
    hists[ipart](ibin, 0) += 1.;
    hists[ipart](ibin, 1) += 1.;
  };

  auto f_analyze_wgt = [&coll, &hists, &bins, &values, &ivars, nvar, iweight, &fwgt] (long long entry) {
    coll.populate(entry);
    double wgt = coll.get<double>(iweight, 0);

    for (int ivar = 0; ivar < nvar; ++ivar)
      values[ivar] = coll.get<double>(ivars[ivar], 0);

    auto ibin = find_bin(values, bins);
    fwgt(hists, ibin, wgt);
  };

  auto f_analyze_one = [&coll, &hists, &bins, &values, &ivars, nvar, &fone] (long long entry) {
    coll.populate(entry);

    for (int ivar = 0; ivar < nvar; ++ivar)
      values[ivar] = coll.get<double>(ivars[ivar], 0);

    auto ibin = find_bin(values, bins);
    fone(hists, ibin);
  };

  if (weighted and iweight != -1)
    dataset.set_analyzer(f_analyze_wgt, true);
  else
    dataset.set_analyzer(f_analyze_one, true);

  dataset.analyze(-1, -1, silent);

  return hists;
}



void make_histogram_set(const std::vector<std::string> &files,
                        const std::string &tree,
                        const std::vector<std::string> &variables, const std::string &weight,
                        const std::string &output)
{
  auto variables_bins = variables_and_binning(variables, weight, files, tree);
  Framework::Dataset<TChain> dataset("dataset", tree);
  dataset.set_files(files);
  auto coll = make_collection(dataset, variables_bins);
  auto histogram = count_and_bin(dataset, coll, variables_bins, 1, 1);

  const auto nbvars = std::make_tuple(std::vector<std::vector<std::string>>{{"variables"}},
                                      std::vector<std::vector<double>>{ make_interval(0., double(std::get<1>(variables_bins).size()), 1.) },
                                      std::vector<std::vector<double>>{}, std::string{});

  save_all_as(output, array_to_root(variables_bins, "", std::get<1>(variables_bins), histogram[0]), array_to_root(nbvars, "nbin", std::get<1>(nbvars), nbin_hist(variables_bins)));
}

#endif
