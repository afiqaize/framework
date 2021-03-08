// -*- C++ -*-
// author: afiq anuar
// short: the front end code for the smoothing
// compile with 'make smoother', assuming the provided Makefile is used

#include "smooth_util.h"

// http://tclap.sourceforge.net/
#include "tclap/CmdLine.h"

int main(int argc, char** argv) {
  TCLAP::CmdLine cmdbase("", ' ', "0.01");
  TCLAP::MultiArg<std::string> cmdfile("", "file", std::string("input filename to be used, assumed to end with .root. ")
                                       + "use the option multiple times to sum up multiple files.", true, "string", cmdbase);
  TCLAP::ValueArg<std::string> cmdtree("", "tree", "name of tree of nominal and systematic files, if present.", true, "", "string", cmdbase);
  TCLAP::MultiArg<std::string> cmdvar("", "variable", std::string("variable(s) to use. ")
                                      + "variable names must correspond to branch names in the tree. "
                                      + "example syntax: 'cHel : nbin;min;max' or 'mtt : edge1;edge2;...;edgeN' " 
                                      + "for fixed/variable binning respectively. "
                                      + "whenever valid, fixed binning interpretation takes precedence. "
                                      + "expressions like 'c = <expression> : binning', can also be used. "
                                      + "accepted expressions are unary operations 'op(a)' where supported unary operations are: "
                                      + "exp, log, log10, sin, cos, tan, asin, acos, atan, sqrt, abs, negate and invert "
                                      + "or binary operations a op b, where supported binary operations are: +, -, * and /. "
                                      + "variables and expressions are interpreted as continuous. "
                                      + "under- and overflows are automatically added to the first and last bins respectively. "
                                      + "use the option multiple times to make a multidimensional histogram.", true, "string", cmdbase);
  TCLAP::MultiArg<std::string> cmdweight("", "weight", std::string("branch name to be used as weight. also accepts expressions (see --variable). ")
                                         + "needs to be called 3 times if --type is weight i.e. systematic variations are in the form of event weights. "
                                         + "in this case, first weight is the nominal and second and third are the up and down variations respectively.",
                                         false, "string", cmdbase);
  TCLAP::ValuesConstraint<std::string> allmode( {"histogram", "smooth"} );
  TCLAP::ValueArg<std::string> cmdmode("", "mode", std::string("runtime mode. currently available: ")
                                       + "histogram: make histograms of the provided files (--systematic related options are ignored) \n "
                                       + "smooth: perform systematic smoothing\n"
                                       + ">1D histograms are unrolled, with the first dimension presented in bins of second variable and so on.",
                                       true, "", &allmode, cmdbase);
  TCLAP::ValueArg<std::string> cmdsyst("", "systematic", std::string("systematic uncertainty name. ")
                                       + "if --type is tree (default), files containing systematically varied events are "
                                       + "assumed to be in the same directory as nominal and "
                                       + "are named filename_systematic_direction.root, where:\n " 
                                       + "filename.root is as given to --file\n "
                                       + "systematic is as given to --systematic,\n "
                                       + "and direction is up and down.", false, "", "string", cmdbase);
  TCLAP::ValuesConstraint<std::string> alltype( {"tree", "weight"} );
  TCLAP::ValueArg<std::string> cmdtype("", "type", std::string("type of systematic variations. can be either tree (default) or weight. ")
                                       + "tree is when the systematic variations are given in separate files (see --systematic) "
                                       + "and weight is when the variations are in the form of event weights.",
                                       false, "tree", &alltype, cmdbase);
  TCLAP::ValueArg<std::string> cmdsnap("", "snapshot", "file name containing previous cross-validation results to be continued from.",
                                         false, "", "string", cmdbase);
  TCLAP::MultiArg<double> cmdbwidth("", "bandwidth", std::string("fixed bandwidth to be used, in case one wants to skip cross-validation. ")
                                    + "it is expressed as a fraction 0 < bw <= 1 and needs to be provided as many times as as --variables. "
                                    + "first bandwidth is assigned to the first variable and so on. "
                                    + "in cross-validation, a range of 0.05 - 1 is tested.", false, "double", cmdbase);
  TCLAP::ValueArg<int> cmdnpart("", "npartition", "number of partitions to be used in cross-validation.", false, 10, "integer", cmdbase);
  TCLAP::ValueArg<int> cmdnrepeat("", "nrepeatcv", "number of times cross-validation is to be repeated.", false, 10, "integer", cmdbase);
  TCLAP::SwitchArg cmdsymm("", "dont-symmetrize", "independently smooth the systematic deviations without symmetrizing.", cmdbase, true);
  TCLAP::SwitchArg cmdbinom("", "binomial", "use binomial approximation in cross-validation.", cmdbase, false);
  TCLAP::ValueArg<std::string> cmdoutput("", "output", "output file name containing the result, assumed to end with .root.",
                                         false, "output.root", "string", cmdbase);
  cmdbase.parse( argc, argv );

  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  if (cmdmode.getValue() == "histogram")
    make_histogram_set(cmdfile.getValue(), cmdtree.getValue(), cmdvar.getValue(), (cmdweight.getValue().empty()) ? "" : cmdweight.getValue()[0],
                       cmdoutput.getValue());
  else if (cmdmode.getValue() == "smooth")
    smooth_templates(cmdfile.getValue(), cmdtree.getValue(), cmdvar.getValue(), cmdweight.getValue(), cmdsyst.getValue(), cmdtype.getValue(),
                     cmdsymm.getValue(), cmdbinom.getValue(), cmdnpart.getValue(), cmdnrepeat.getValue(), cmdsnap.getValue(), cmdbwidth.getValue(),
                     cmdoutput.getValue());
  return 0;
}

