// -*- C++ -*-
// author: afiq anuar
// short: the front end code for the smoothing


// make smoother FWK_BASE_DIR="/nfs/dust/cms/user/afiqaize/cms/bpark_nano_200218/cmssw_1103_analysis/src/fwk"

#include "smooth_util.h"

// http://tclap.sourceforge.net/
#include "tclap/CmdLine.h"

int main(int argc, char** argv) {
  TCLAP::CmdLine cmdbase("", ' ', "0.01");
  TCLAP::MultiArg<std::string> cmdfile("", "file", std::string("input filename to be used, assumed to end with .root. ")
                                       + "use the option multiple times to sum up multiple files.", true, "string", cmdbase);
  TCLAP::ValueArg<std::string> cmdtree("", "tree", "name of tree of nominal and systematic files, if present", true, "", "string", cmdbase);
  TCLAP::ValueArg<std::string> cmdsyst("", "systematic", std::string("systematic uncertainty name. ")
                                       + "assumes that the files containing systematically varied events are of the form "
                                       + "filename_systematic_direction.root, where: " 
                                       + "filename.root is as given to --file, "
                                       + "systematic is as given to --systematic, "
                                       + "and direction is up and down (both files expected to be present).", false, "", "string", cmdbase);
  TCLAP::ValueArg<std::string> cmdmode("", "mode", std::string("mode runtime mode. currently available: ")
                                       + "histogram: make (unrolled, if >1D) histograms of the provided files. --systematic is ignored in this case "
                                       + "smooth: perform systematic smoothing",
                                       true, "", "string", cmdbase);
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
  TCLAP::ValueArg<std::string> cmdweight("", "weight", "branch name to be used as weight. also accepts expressions (see --variable)",
                                         false, "", "string", cmdbase);
  TCLAP::ValueArg<std::string> cmdoutput("", "output", "output file name containing the result, assumed to end with .root",
                                         false, "output.root", "string", cmdbase);
  TCLAP::MultiArg<double> cmdbwidth("", "bandwidth", std::string("fixed bandwidth to be used, in case one wants to skip cross-validation. ")
                                    + "it is expressed as a fraction 0 < bw <= 1. needs to be called as many times as as --variables. "
                                    + "first bandwidth is assigned to the first variable and so on. ", false, "double", cmdbase);
  TCLAP::ValueArg<int> cmdnpart("", "npartition", "number of partitions to be used in cross-validation.", false, 10, "integer", cmdbase);
  TCLAP::ValueArg<int> cmdnrepeat("", "nrepeatcv", "number of times cross-validation is to be repeated.", false, 100, "integer", cmdbase);
  TCLAP::SwitchArg cmdsymm("", "dont-symmetrize", "independently smooth the systematic deviations without symmetrizing.", cmdbase, true);
  TCLAP::SwitchArg cmdbinom("", "binomial", "use binomial approximation in cross-validation.", cmdbase, false);
  cmdbase.parse( argc, argv );

  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  auto variables_bins = variables_and_binning(cmdvar.getValue(), cmdweight.getValue());

  if (cmdmode.getValue() == "histogram") {
    auto histogram = make_histogram_set(cmdfile.getValue(), cmdtree.getValue(), variables_bins);
    save_all_as(cmdoutput.getValue(), array_to_root(variables_bins, "", std::get<1>(variables_bins), histogram[0]));
  }
  else if (cmdmode.getValue() == "smooth") {
    auto result = smooth_templates(cmdfile.getValue(), cmdtree.getValue(), variables_bins, cmdsyst.getValue(),
                                   cmdsymm.getValue(), cmdbinom.getValue(), cmdnpart.getValue(), cmdnrepeat.getValue(), cmdbwidth.getValue());
    save_all_as(cmdoutput.getValue(), result); // what about plotting? a separate script?
  }
  else
    std::cerr << "Mode " << cmdmode.getValue() << " is not yet implemented" << std::endl;

  return 0;
}

