// -*- C++ -*-
// author: afiq anuar
// short: provides function definitions needed to make histograms from flat trees

#ifndef FWK_MAKE_HISTOGRAM_H
#define FWK_MAKE_HISTOGRAM_H

#include "Dataset.h"
#include "Collection.h"

#include <cmath>
#include "string_io.h"
#include "container_util.h"
#include "function_util.h"
#include "rng_util.h"
#include "constants.h"

#include "array_histogram.h"
#include "output_util.h"

using FwkColl = Framework::Collection<boolean,
                                      int, unsigned int, float,
                                      long long, unsigned long long, double>;

bool valid_name(const std::string &str)
{
  const auto is_valid = [] (unsigned char c) { return std::isalnum(c) or c == '_'; };
  return std::all_of(std::begin(str), std::end(str), is_valid);
}



/// to translate the argument of constant argument to unary constant(...) into something that valid_name() returns true for
/// but looks awful enough that presumably no one uses it as an actual branch name
std::string to_constant_str(const std::string &str)
{
  // not that these strings have to be explainable, but human readability is nice
  auto copy = "__cv__"s + str; // constant value
  replace(copy, "-"s, "__ns__"s); // negative sign
  replace(copy, "+"s, ""s);
  replace(copy, "."s, "__fp__"s); // floating point
  return copy;
}



/// the reverse operation of the above
/// no check is made whether it's actually a valid constant string - this is merely an impl detail
double constant_str_value(const std::string &str)
{
  auto copy = str;
  replace(copy, "__fp__"s, "."s);
  replace(copy, "__ns__"s, "-"s);
  replace(copy, "__cv__"s, ""s);
  return std::stod(copy);
}



/// is a string elementary or not; if yes, returns a vector of string containing its name, the operation involved and the operands
/// op1(a) and a op2 b are both elementary, anything else are not
/// op1/op2 can be any of the supported unary/binary operations
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
  static const std::vector<std::string> unaries = {"constant("s, "alias("s,
                                                   "exp("s, "log("s, "log10("s,
                                                   "sin("s, "cos("s, "tan("s, "asin("s, "acos("s, "atan("s, "sqrt("s, "abs("s,
                                                   "negate("s, "relu("s, "step("s, "invert("s, "not("s};
  static const std::vector<std::string> binaries = {"*"s, "/"s, "+"s, "-"s};

  // decompose expressions like a + b + c + d + e -> f + c + d + e -> g + d + e -> h + e...
  // the refs are what needs to be kept track of for the inplace edits within loop body
  // iop and icl are opening and closing indices of the part of string to analyze 
  auto decompose_binaries = [] (std::vector<std::vector<std::string>> &expressions, std::string &exp, std::string &part,
                                std::vector<std::string> &elem, std::string::size_type iop, std::string::size_type icl) {
    auto imd = min(exp.find("*"s, iop), exp.find("/"s, iop));
    if (imd > icl)
      imd = std::string::npos;

    if (imd != std::string::npos) {
      auto ipm1 = min(exp.find("+"s, iop), exp.find("-"s, iop));
      auto ibi = min(exp.find("+"s, imd + 1), exp.find("-"s, imd + 1), exp.find("*"s, imd + 1), exp.find("/"s, imd + 1));

      if (ipm1 > imd)
        part = exp.substr(iop, (ibi < icl) ? ibi - iop : icl - iop);
      else {
        auto ipm2 = min(exp.find("+"s, ipm1 + 1), exp.find("-"s, ipm1 + 1));
        while (ipm2 < imd) {
          ipm1 = ipm2;
          ipm2 = min(exp.find("+"s, ipm1 + 1), exp.find("-"s, ipm1 + 1));
        }

        part = exp.substr(ipm1 + 1, (ibi < icl) ? ibi - ipm1 - 1 : icl - ipm1 - 1);
      }

      auto belem = is_elementary(part, unaries, binaries);
      if (not belem.empty()) {
        unique_emplace(expressions, belem, exp, part);
        elem = is_elementary(exp, unaries, binaries);
        return true;
      }
    }
    else {
      auto ipm1 = min(exp.find("+"s, iop), exp.find("-"s, iop));
      if (ipm1 > icl)
        imd = std::string::npos;

      if (ipm1 != std::string::npos) {
        auto ipm2 = min(exp.find("+"s, ipm1 + 1), exp.find("-"s, ipm1 + 1));
        part = exp.substr(iop, (ipm2 < icl) ? ipm2 - iop : icl - iop);

        auto belem = is_elementary(part, unaries, binaries);
        if (not belem.empty()) {
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
    if (exp == ""s)
      continue;

    // constants parsing in expression need some hacks, do that first
    // done by transforming them a technically valid (but looking so nonsensical that hopefully no none uses it) variable name
    auto iconstant = exp.find(unaries[0]);
    while (iconstant != std::string::npos) {
      auto iclose = exp.find(")"s, iconstant + unaries[0].length());
      auto real_value_str = exp.substr(iconstant + unaries[0].length(), iclose - iconstant - unaries[0].length());
      auto value_str = to_constant_str(real_value_str);
      replace(exp, real_value_str, value_str, iconstant);
      iconstant = exp.find(unaries[0], exp.find(value_str, iconstant));
    }

    // get all the 'raw' attribute names and dump them into a list, to be added to the Collection later
    // by getting rid of all the unaries, and grabbing valid variable names
    auto cexp = exp;
    for (const auto &unary : unaries)
      strip(cexp, unary);

    auto branches = split_if(cexp, [] (unsigned char c) { return std::isalnum(c) or c == '_'; });
    for (auto &branch : branches) {
      if (std::find_if(std::begin(expressions), std::end(expressions), [&branch] (auto &var) { return var[0] == branch; }) == std::end(expressions))
        expressions.emplace_back(std::vector<std::string>{branch, ""s});
    }

    auto elem = is_elementary(exp, unaries, binaries);
    std::string part = ""s;
    while (elem.empty()) {
      // a flag to indicate that an elementary expression is found and dumped into expression list
      // due to the inplace editing, we want to be analyzing it in a fresh iteration
      bool dumped = false;

      // find innermost opening and closing brackets
      auto nop = count_substring(exp, "("s), ncl = count_substring(exp, ")"s);
      if (nop != ncl) {
        std::cerr << "Expression "s << exp << " is invalid due to non-matching parenthesis. skipping..."s << std::endl;
        return;
      }

      if (nop > 0) {
        auto iop1 = exp.find("("s), icl = exp.find(")"s);
        auto iop2 = exp.find("("s, iop1 + 1);
        while (iop2 < icl) {
          iop1 = iop2;
          iop2 = exp.find("("s, iop1 + 1);
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
        if (not decompose_binaries(expressions, exp, part, elem, iop1 + 1, icl)) {
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
/// segfaults at dataset.analyze() if the source promotion is done within make_collection for yet to be understood reasons
auto variables_and_binning(const std::vector<std::string> &variables, const std::string &weight)
{
  using namespace Framework;

  static const std::vector<int> fsplits = {12, 5, 3, 2};
  auto fsplit = (variables.size() > 2) ? fsplits.back() : fsplits[variables.size() - 1];

  std::vector<std::vector<std::string>> expressions;
  std::vector<std::vector<double>> cedges, fedges;

  for (const auto &var : variables) {
    auto veb = split(var, ":"s);
    if (veb.size() != 2) {
      std::cerr << "Variable or expression "s << var << " is invalid"s << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{}, ""s);
    }

    auto ve = split(veb[0], "="s);
    if (ve.size() != 1 and ve.size() != 2) {
      std::cerr << "Variable or expression "s << veb[0] << " is invalid"s << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{}, ""s);
    }
    strip(ve[0]);
    if (not valid_name(ve[0])) {
      std::cerr << "Invalid variable name "s << ve[0] << ". Aborting. Current version considers only names containing " 
        "alphanumeric characters or underscores"s << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{}, ""s);
    }

    if (ve.size() == 2)
      strip(ve[1]);
    else
      ve.emplace_back(""s);
    expressions.emplace_back(std::move(ve));

    auto binstyle = split(veb[1], ";"s);
    int nbinf = (binstyle.size() == 2) ? std::stoi(binstyle[0]) : 0;
    const auto [minf, maxf] = [&binstyle] () {
      const auto minmax = (binstyle.size() == 2) ? split(binstyle[1], ","s) : std::vector<std::string>{};
      return (minmax.size() == 2) ? std::make_pair(std::stod(minmax[0]), std::stod(minmax[1])) : std::make_pair(0., 0.); 
    } ();

    const auto cedge = [&binning = veb[1]] {
      std::vector<double> cedge;
      if (contain(binning, ";"s))
        return cedge;

      const auto binrng = split(binning, ","s);
      for (const auto &edge : binrng)
        cedge.emplace_back(std::stod(edge));
      return cedge;
    }();

    if (nbinf > 1 and minf < maxf) {
      cedges.emplace_back( make_interval(minf, maxf, (maxf - minf) / nbinf) );
      fedges.emplace_back( make_interval(minf, maxf, (maxf - minf) / (fsplit * nbinf)) );
    }
    else if (cedge.size() > 2 and std::is_sorted(std::begin(cedge), std::end(cedge))) {
      fedges.emplace_back( [&cedge, fsplit] {
        std::vector<double> fedge = {cedge[0]};
        for (int ic = 0; ic < cedge.size() - 1; ++ic)
          fill_interval(fedge, cedge[ic + 1], (cedge[ic + 1] - cedge[ic]) / fsplit);

        return fedge;
      }() );
      cedges.emplace_back( std::move(cedge) );
    }
    else {
      std::cerr << "Binning information for variable  "s << ve[0] << " is invalid. Note that variables with only one bin is not accepted"s << std::endl;
      return std::make_tuple(std::vector<std::vector<std::string>>{}, std::vector<std::vector<double>>{}, std::vector<std::vector<double>>{}, ""s);
    }
  }

  auto we = split(weight, "="s);
  strip(we[0]);
  if (weight != ""s) {
    if (valid_name(we[0])) {
      if (we.size() == 2)
        strip(we[1]);
      else
        we.emplace_back(""s);
      expressions.emplace_back(we);
    }
    else {
      std::cerr << "Invalid weight name "s << we[0] << ". Ignoring the weight. Current version considers only names containing " 
        "alphanumeric characters or underscores"s << std::endl;
    }
  }

  /**** start testing of expression parser ****
  for (auto &exp : expressions) {
    for (auto &var : exp)
      std::cout << var << " -- "s;
    std::cout << std::endl;
  }
  std::cout << "     ----------------     "s << std::endl;
  ****/
  parse_expressions(expressions);
  /****
  for (auto &exp : expressions) {
    for (auto &var : exp)
      std::cout << var << " -- "s;
    std::cout << std::endl;
  }
  std::cout << "     ----------------     "s << std::endl;
  ****  end testing of expression parser  ****/

  return std::make_tuple(std::move(expressions), std::move(cedges), std::move(fedges), we[0]);
}



/// prepare the collection to read the variables
auto make_collection(const std::tuple<
                     std::vector<std::vector<std::string>>,
                     std::vector<std::vector<double>>,
                     std::vector<std::vector<double>>,
                     std::string> &variables_bins,
                     const std::vector<std::string> &files, const std::string &tree)
{
  using namespace Framework;

  const auto &variables = std::get<0>(variables_bins);
  if (variables.empty())
    throw std::runtime_error( "ERROR: make_histogram::make_collection: variables list is empty. Aborting."s );

  // initialize the collection, and add the source attributes ie those read in from root files
  FwkColl coll("coll"s, variables.size());
  for (const auto &var : variables) {
    if (var.size() == 2 and var[1] == ""s and not contain(var[0], "__cv__"s))
      coll.add_attribute(var[0], var[0]);
  }

  // figure out the types of the attributes based on the first file in dataset
  Dataset<TChain> dataset("dataset"s, tree);
  dataset.set_files(files, 1);
  dataset.associate(coll);
  coll.detach();

  // constant(...) declares a new attribute in the collection, but doesn't introduce one in the variable list
  // we correct for it, and ensure all other variables correspond to attributes in the collection
  const int nconst = std::count_if( std::begin(variables), std::end(variables), [] (const auto &var) {return contain(var[0], "__cv__"s);} );
  while (coll.n_attributes() + nconst != variables.size()) {
    for (const auto &var : variables) {
      if (var.size() > 2) {
        bool isthere = true;

        for (int iarg = 2; iarg < var.size(); ++iarg)
          isthere = isthere and (coll.has_attribute(var[iarg]) or var[1] == "constant("s);

        if (isthere) {
          if (coll.n_attributes() > 0 and var[1] == "constant("s) {
            const auto attrs = coll.attributes();
            std::visit([&coll, &var, &attr = attrs[0]] (auto &&arg) {
                using arg_type = typename std::decay_t<decltype(arg)>::value_type;
                coll.transform_attribute(var[0], [value = constant_str_value(var[2])] (arg_type _) -> double {(void) _; return value;}, attr);
              }, coll(attrs[0]));
          }
          else {
            if (var.size() > 3) {
              std::visit([&coll, &var] (auto &&arg1, auto &&arg2) {
                           using arg_type1 = typename std::decay_t<decltype(arg1)>::value_type;
                           using arg_type2 = typename std::decay_t<decltype(arg2)>::value_type;

                           if (var[1] == "+"s)
                             coll.transform_attribute(var[0], binaries::add<arg_type1, arg_type2>(), var[2], var[3]);
                           else if (var[1] == "-"s)
                             coll.transform_attribute(var[0], binaries::subtract<arg_type1, arg_type2>(), var[2], var[3]);
                           else if (var[1] == "*"s)
                             coll.transform_attribute(var[0], binaries::multiply<arg_type1, arg_type2>(), var[2], var[3]);
                           else if (var[1] == "/"s)
                             coll.transform_attribute(var[0], binaries::divide<arg_type1, arg_type2>(), var[2], var[3]);
              }, coll(var[2]), coll(var[3]));
            }
            else if (var.size() > 2) {
              std::visit([&coll, &var] (auto &&arg) {
                           using arg_type = typename std::decay_t<decltype(arg)>::value_type;

                           if (var[1] == "alias("s)
                             coll.alias_attribute(var[0], var[2]);
                           else if (var[1] == "exp("s)
                             coll.transform_attribute(var[0], unaries::exp<arg_type>(), var[2]);
                           else if (var[1] == "log("s)
                             coll.transform_attribute(var[0], unaries::log<arg_type>(), var[2]);
                           else if (var[1] == "log10("s)
                             coll.transform_attribute(var[0], unaries::log10<arg_type>(), var[2]);
                           else if (var[1] == "sin("s)
                             coll.transform_attribute(var[0], unaries::sin<arg_type>(), var[2]);
                           else if (var[1] == "cos("s)
                             coll.transform_attribute(var[0], unaries::cos<arg_type>(), var[2]);
                           else if (var[1] == "tan("s)
                             coll.transform_attribute(var[0], unaries::tan<arg_type>(), var[2]);
                           else if (var[1] == "asin("s)
                             coll.transform_attribute(var[0], unaries::asin<arg_type>(), var[2]);
                           else if (var[1] == "acos("s)
                             coll.transform_attribute(var[0], unaries::acos<arg_type>(), var[2]);
                           else if (var[1] == "atan("s)
                             coll.transform_attribute(var[0], unaries::atan<arg_type>(), var[2]);
                           else if (var[1] == "sqrt("s)
                             coll.transform_attribute(var[0], unaries::sqrt<arg_type>(), var[2]);
                           else if (var[1] == "abs("s)
                             coll.transform_attribute(var[0], unaries::abs<arg_type>(), var[2]);
                           else if (var[1] == "negate("s)
                             coll.transform_attribute(var[0], unaries::negate<arg_type>(), var[2]);
                           else if (var[1] == "relu("s)
                             coll.transform_attribute(var[0], unaries::relu<arg_type>(), var[2]);
                           else if (var[1] == "step("s)
                             coll.transform_attribute(var[0], unaries::step<arg_type>(), var[2]);
                           else if (var[1] == "invert("s)
                             coll.transform_attribute(var[0], unaries::invert<arg_type>(), var[2]);
                           else if (var[1] == "not("s)
                             coll.transform_attribute(var[0], unaries::contra<arg_type>(), var[2]);
              }, coll(var[2]));
            }
          }
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
                                     FwkColl &coll,
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
    for (int ivar = 0; ivar < nvar; ++ivar) {
      ivars[ivar] = coll.inquire(variables[ivar][0]);
    }

    return ivars;
  }();

  if (std::count(std::begin(ivars), std::end(ivars), -1))
    throw std::runtime_error( "ERROR: make_histogram::count_and_bin: variable indices point to nonexistent attributes. should never happen."s );

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

    auto wgt = coll.convert<double>(iweight, 0);
    for (int ivar = 0; ivar < nvar; ++ivar)
      values[ivar] = coll.convert<double>(ivars[ivar], 0);

    auto ibin = find_bin(values, bins);
    fwgt(hists, ibin, wgt);
  };

  auto f_analyze_one = [&coll, &hists, &bins, &values, &ivars, nvar, &fone] (long long entry) {
    coll.populate(entry);

    for (int ivar = 0; ivar < nvar; ++ivar)
      values[ivar] = coll.convert<double>(ivars[ivar], 0);

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



/// make a single histogram with nvar bins
/// and each bin i has a content of nbin of ith variable
Arrayhist nbin_hist(const std::tuple<
                    std::vector<std::vector<std::string>>,
                    std::vector<std::vector<double>>,
                    std::vector<std::vector<double>>,
                    std::string> &variables_bins)
{
  const auto &bins = std::get<1>(variables_bins);
  const int nvar = bins.size();

  Arrayhist hist(nvar);
  for (int iv = 0; iv < nvar; ++iv) {
    hist(iv, 0) = bins[iv].size() - 1;
    hist(iv, 1) = 0.;
  }

  return hist;
}



void make_histogram_set(const std::vector<std::string> &files,
                        const std::string &tree,
                        const std::vector<std::string> &variables, const std::string &weight,
                        char restype, const std::string &output)
{
  auto variables_bins = variables_and_binning(variables, weight);
  Framework::Dataset<TChain> dataset("dataset"s, tree);
  dataset.set_files(files);
  auto coll = make_collection(variables_bins, files, tree);
  auto histogram = count_and_bin(dataset, coll, variables_bins, 1, 1);

  const auto nbvars = std::make_tuple(std::vector<std::vector<std::string>>{{"variables"s}},
                                      std::vector<std::vector<double>>{ make_interval(0., double(std::get<1>(variables_bins).size()), 1.) },
                                      std::vector<std::vector<double>>{}, std::string{});

  save_all_as(output, array_to_root(variables_bins, ""s, std::get<1>(variables_bins), histogram[0], restype), array_to_root(nbvars, "nbin"s, std::get<1>(nbvars), nbin_hist(variables_bins)));
}

#endif
