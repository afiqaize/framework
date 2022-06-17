#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RSnapshotOptions.hxx"

#include "misc/prescaler.h"
#include "misc/function_util.h"
#include "misc/string_io.h"

// 

int main() {
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df_("Events", "/eos/cms/store/group/phys_b2g/nbinnorj/TTBParkNano/v1p0/MERGEDv2/NanoAODMerged_DataUL18A_ParkingBPH1.root");

  auto df = ROOT::RDF::RNode(df_);
  df = df.Filter("run == 315974");

  // get the list of paths, and seeds
  auto ps = Prescaler<>("/afs/cern.ch/work/a/afiqaize/xxx/framework/bpark/prescales_parkingbph.json", {"HLT_Mu10p5_IP3p5_part0", "HLT_Mu8p5_IP3p5_part0"});
  std::vector<std::string> paths = ps.hlt_paths();
  std::vector<std::string> seeds = ps.l1_seeds();

  // I screwed up the json version clipping somehow, so...
  for (auto &p : paths)
    replace(p, "_v", "");

  // prepare the functions that does the booleans -> bitset conversions, and find out what type the bitset is
  auto f_hlt_bits = to_bitset<Bool_t, 2>(); // FIXME size isn't a compile-time number, need to be pre-evaluated
  auto f_l1_bits = to_bitset<Bool_t, 6>(); // FIXME size isn't a compile-time number, need to be pre-evaluated
  using bit_type = typename function_traits<decltype(f_hlt_bits)>::result_type; // assuming the above has the same type (by default std::bitset<128>)

  // and make a multithreading-friendly version of a function that calculates the prescale weights
  auto f_ps_weight = [&prescaler = ps] (uint run, uint lumi, const bit_type &hlt, const bit_type &l1) {
    thread_local auto ps = prescaler; // make a thread local copy so that the weights retrieval can be thread-safe
    return ps.weight(run, lumi, hlt, l1);
  };

  // define the columns for hlt and l1 bitsets
  df = df.Define("hlt_bits", f_hlt_bits, paths);
  df = df.Define("l1_bits", f_l1_bits, seeds);

  // and lastly the column for the ps weight
  df = df.Define("ps_weight", f_ps_weight, {"run", "luminosityBlock", "hlt_bits", "l1_bits"});

  df = df.Filter("ps_weight > 0 && nMuon >= 2 && nElectron >= 1");

  df.Snapshot("events", "lolk.root", {"ps_weight", "nMuon", "nElectron"});

  auto report = df.Report();
  report->Print();

  ROOT::RDF::SaveGraph(df, "graph_skim.dot");
}
