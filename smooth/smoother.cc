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
  TCLAP::ValueArg<std::string> cmdtree("", "tree", "tree name of nominal and systematic files, if present.", true, "", "string", cmdbase);
  TCLAP::MultiArg<std::string> cmdvar("", "variable", std::string("variable(s) to use. ")
                                      + "variable names must correspond to branch names in the tree. "
                                      + "they must contain only alphanumeric or underscore characters. "
                                      + "example syntax: 'cHel : nbin;min,max' or 'mtt : edge1,edge2,...,edgeN' " 
                                      + "for fixed/variable binning respectively. "
                                      + "expressions like 'c = <expression> : binning', can also be used. "
                                      + "accepted expressions are unary operations 'op(a)' or binary operations 'a op b'. "
                                      + "supported unary operations: exp, log, log10, sin, cos, tan, asin, acos, atan, sqrt, abs, negate and invert. "
                                      + "supported binary operations: +, -, * and /. "
                                      + "variables and expressions are interpreted as floating-point. "
                                      + "under- and overflows are automatically added to the first and last bins respectively. "
                                      + "use the option multiple times to make a multidimensional histogram.", true, "string", cmdbase);
  TCLAP::MultiArg<std::string> cmdweight("", "weight", std::string("branch name to be used as weight. also accepts expressions (see --variable). ")
                                         + "has to be called 3 times if --type is weight i.e. systematic variations are in the form of event weights. "
                                         + "in this case, first weight is the nominal and second and third are the up and down variations respectively.",
                                         false, "string", cmdbase);
  TCLAP::ValuesConstraint<std::string> allmode( {"histogram", "systematic", "smooth"} );
  TCLAP::ValueArg<std::string> cmdmode("", "mode", std::string("runtime mode. currently available: ")
                                       + "histogram: make histograms of the provided files (--systematic related options are ignored) "
                                       + "systematic: make nominal and systematic histograms "
                                       + "smooth: as systematic and also perform smoothing on the systematic histograms "
                                       + ">1D histograms are unrolled, with the first dimension presented in bins of second variable and so on.",
                                       true, "", &allmode, cmdbase);
  TCLAP::ValueArg<std::string> cmdsyst("", "systematic", std::string("systematic uncertainty name. ")
                                       + "if --type is tree, files containing systematically varied events are "
                                       + "assumed to be in the same directory as nominal and "
                                       + "are named filename_systematic_direction.root, where: " 
                                       + "filename.root is as given to --file, "
                                       + "systematic is as given to --systematic, "
                                       + "and direction is up and down.", false, "", "string", cmdbase);
  TCLAP::ValuesConstraint<std::string> alltype( {"tree", "weight"} );
  TCLAP::ValueArg<std::string> cmdtype("", "type", std::string("type of systematic variations. ")
                                       + "tree is when the systematic variations are given in separate files (see --systematic) "
                                       + "and weight is when the variations are in the form of event weights.",
                                       false, "tree", &alltype, cmdbase);
  TCLAP::ValueArg<std::string> cmdsnap("", "snapshot", "file containing previous cross-validation run to be continued from.",
                                       false, "", "string", cmdbase);
  TCLAP::ValuesConstraint<std::string> allres( {"unroll", "normal", "all"} );
  TCLAP::ValueArg<std::string> cmdrestype("", "result-type", std::string("type of the result i.e. output histograms. ")
                                          + "unroll means >1D histograms are laid in slices of the first variable in the first bin of second and so on. "
                                          + "in this case bin edges information is lost for >1D histogram. "
                                          + "normal means 2D and 3D histogram is written out without unrolling. "
                                          + "all is unroll + normal. "
                                          + "default result-type is unroll.",
                                          false, "unroll", &allres, cmdbase);
  TCLAP::MultiArg<double> cmdbwidth("", "bandwidth", std::string("fixed bandwidth to be used, in case one wants to skip cross-validation. ")
                                    + "it is expressed as a fraction 0 < bw <= 1 and has to be provided as many times as as --variables. "
                                    + "first bandwidth is assigned to the first variable and so on. "
                                    + "in cross-validation, a range of 0.05 - 1 is tested.", false, "double", cmdbase);
  TCLAP::ValueArg<int> cmdnpart("", "npartition", "number of partitions to be used in cross-validation.", false, 10, "integer", cmdbase);
  TCLAP::ValueArg<int> cmdnrepeat("", "nrepeatcv", "number of times cross-validation is to be repeated.", false, 10, "integer", cmdbase);
  TCLAP::SwitchArg cmdoneside("", "one-side", std::string("when the systematic uncertainty is one-sided rather than up and down. ")
                              + "if --type is tree then the expected systematic file name is filename_systematic.root "
                              + "and if --type is weight then --weight needs to be called 2 times.", cmdbase, false);
  TCLAP::SwitchArg cmdsymm("", "dont-symmetrize", "independently smooth the systematic deviations without symmetrizing.", cmdbase, true);
  TCLAP::SwitchArg cmdbinom("", "binomial", "use binomial approximation in cross-validation.", cmdbase, false);
  TCLAP::ValueArg<std::string> cmdoutput("", "output", "output file name containing the result, assumed to end with .root.",
                                         false, "output.root", "string", cmdbase);
  cmdbase.parse( argc, argv );

  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  if (cmdmode.getValue() == "histogram")
    make_histogram_set(cmdfile.getValue(), cmdtree.getValue(), cmdvar.getValue(), (cmdweight.getValue().empty()) ? "" : cmdweight.getValue()[0],
                       cmdrestype.getValue()[0], cmdoutput.getValue());
  else if (cmdmode.getValue() == "systematic" or cmdmode.getValue() == "smooth") {
    const bool dosmooth = cmdmode.getValue() == "smooth";
    smooth_templates(cmdfile.getValue(), cmdtree.getValue(), cmdvar.getValue(), cmdweight.getValue(), cmdsyst.getValue(), cmdtype.getValue(), dosmooth,
                     cmdoneside.getValue(), cmdsymm.getValue(), cmdbinom.getValue(), cmdnpart.getValue(), cmdnrepeat.getValue(), cmdrestype.getValue()[0],
                     cmdsnap.getValue(), cmdbwidth.getValue(), cmdoutput.getValue());
  }
  return 0;
}

