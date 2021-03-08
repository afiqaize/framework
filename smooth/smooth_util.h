// -*- C++ -*-
// author: afiq anuar
// short: functions related to the lowess smoothing and computing its inputs

#ifndef FWK_SMOOTH_UTIL_H
#define FWK_SMOOTH_UTIL_H

#include <thread>

#include "misc/container_util.h"

#include "make_histogram.h"
#include "fit_util.h"
#include "output_util.h"

auto bandwidth_to_test(const std::vector<std::vector<double>> &edges, const std::vector<double> &fixed_bandwidth = {})
{
  // build up the list bandwidths to be tested for each partition
  std::vector<std::vector<int>> bandwidths;

  if (fixed_bandwidth.size() == edges.size()) {
    bandwidths.emplace_back(std::vector<int>(fixed_bandwidth.size(), 0));

    for (int iv = 0; iv < edges.size(); ++iv)
      bandwidths[0][iv] = std::ceil(fixed_bandwidth[iv] * (edges[iv].size() - 1));

    return bandwidths;
  }

  const std::vector<double> grid = {0.05, 0.1, 0.2, 0.5, 1.}; 
  std::vector<std::vector<double>> dummy(edges.size(), std::vector<double>(grid.size() + 1, 0.));

  const int nbw = std::pow(grid.size(), edges.size());
  bandwidths.reserve(nbw);

  for (int ibw = 0; ibw < nbw; ++ibw) {
    auto bandwidth = std::vector<int>(edges.size(), 0);

    for (int iv = 0; iv < edges.size(); ++iv)
      bandwidth[iv] = std::ceil(grid[index_1n(ibw, iv, dummy)] * (edges[iv].size() - 1));

    if (std::all_of(std::begin(bandwidth), std::end(bandwidth), [] (auto &bw) { return bw > 1; }))
      bandwidths.emplace_back(bandwidth);
  }
  bandwidths.shrink_to_fit();

  return bandwidths;
}



auto relative_deviation(const Arrayhist &nominal, const Arrayhist &varied)
{
  auto deviation = Arrayhist(nominal.size());
  for (int ibin = 0; ibin < nominal.size(); ++ibin) {
    if (nominal(ibin, 0) > 0. and varied(ibin, 0) > 0.) {
      deviation(ibin, 0) = (varied(ibin, 0) / nominal(ibin, 0)) - 1.;

      deviation(ibin, 1) = varied(ibin, 1) / (varied(ibin, 0) * varied(ibin, 0));
      deviation(ibin, 1) += nominal(ibin, 1) / (nominal(ibin, 0) * nominal(ibin, 0));
      deviation(ibin, 1) *= (deviation(ibin, 0) + 1.) * (deviation(ibin, 0) + 1.);
    }
    else {
      deviation(ibin, 0) = 0.;
      deviation(ibin, 0) = 0.;
    }
  }

  return deviation;
}



auto symmetrize_variation(const Arrayhist &vary_h, const Arrayhist &vary_l)
{
  auto symmetric = Arrayhist(vary_h.size());
  for (int ibin = 0; ibin < vary_h.size(); ++ibin) {
    symmetric(ibin, 0) = 0.5 * (vary_h(ibin, 0) - vary_l(ibin, 0));
    symmetric(ibin, 1) = 0.25 * (vary_h(ibin, 1) + vary_l(ibin, 1));
  }

  return symmetric;
}



/// reverse operation of relative_deviation
/// for given nominal template and a set of relative deviation, constructs the varied histogram
auto apply_deviation(const Arrayhist &nominal, const Arrayhist &deviation, double scale = 1.)
{
  auto varied = Arrayhist(nominal.size());
  for (int ibin = 0; ibin < nominal.size(); ++ibin) {
    const double scd = (scale * deviation(ibin, 0));
    const double scd2 = scd * scd;
    varied(ibin, 0) = (scd + 1.) * nominal(ibin, 0);
    varied(ibin, 1) = (nominal(ibin, 1) / nominal(ibin, 0) / nominal(ibin, 0)) + (scale * scale * deviation(ibin, 1) / scd2); // assume uncorr, but hmm...
    varied(ibin, 1) *= scd2 * nominal(ibin, 0) * nominal(ibin, 0);
  }

  return varied;
}



/// fit a plane to a list of points
/// ensure that the workspace is at least as large as the entire template
std::vector<std::array<double, 2>> fit_plane(const std::vector<int> &bins, const Arrayhist &reldev,
                                             const std::vector<std::vector<double>> &edges, Workspace &wsp)
{
  const int nvar = edges.size();
  std::vector<int> idx1(nvar, -1);
  for (int iv = 0; iv < nvar; ++iv)
    idx1[iv] = index_1n(bins[0], iv, edges);
  const double rdev1 = reldev(bins[0], 0);

  // count the number of valid points to be included in the fit
  // valid points are either the first/center point and all other points whose value aren't 0
  std::vector<int> ivalid = {};
  for (int ibin = 0; ibin < bins.size(); ++ibin)
    if (reldev(bins[ibin], 0) != 0.)
      ivalid.emplace_back(bins[ibin]);

  std::vector<std::array<double, 2>> result(bins.size(), {rdev1, reldev(bins[0], 1)});

  // too few valid points, abort
  // first point i.e. plane center kept as-is
  if (ivalid.size() < nvar) {
    result.resize(1);
    return result;
  }

  // build up the needed inputs for the fit
  // values is the relative deviation, weights the inverse of its variance
  // response is something like response matrix, relating values to the bin indices
  // covariance is a matrix the gsl routine needs to store some result
  // values and idxn are 'centralized' i.e. expressed relative to the first point
  // this allows reducing the number of parameters to be fitted by 1
  double chi2 = 0.; 
  Vector values(ivalid.size()), weights(ivalid.size()), params(nvar), idxn(nvar);
  Matrix response(ivalid.size(), nvar), covariance(nvar, nvar);

  for (int ival = 0; ival < ivalid.size(); ++ival) {
    values.set(ival, reldev(ivalid[ival], 0) - rdev1);
    weights.set(ival, 1. / reldev(ivalid[ival], 1));

    for (int iv = 0; iv < nvar; ++iv)
      response.set(ival, iv, index_1n(ivalid[ival], iv, edges) - idx1[iv]);
  }

  // the fit itself - returns an int, status code of some sort?
  // chi2 value probably won't matter...
  gsl_multifit_wlinear(response, weights, values, params, covariance, &chi2, wsp);

  // update result with smooth estimates from the fit and undo the centralization
  for (int ibin = 0; ibin < bins.size(); ++ibin) {
    for (int iv = 0; iv < nvar; ++iv)
      idxn.set(iv, index_1n(bins[ibin], iv, edges) - idx1[iv]);

    // also returns an int, maybe some status code
    gsl_multifit_linear_est(idxn, params, covariance, &result[ibin][0], &result[ibin][1]);
    result[ibin][0] += rdev1;
    result[ibin][1] *= result[ibin][1];
  }

  ivalid.clear();
  return result;
}



/// compute weight of point2 relative to point1, within a given binning defined by edges
/// the weight is given by the tricubic function of the distance between the points
/// where the distance are normalized such that points at opposing corners have distance 1
double weight(const std::vector<int> &point1, const std::vector<int> &point2, const std::vector<std::vector<double>> &edges)
{
  if (point1 == point2)
    return 1.;

  static const auto maxsqd = [&edges, npnt = point1.size()] () {
    const auto origin = std::vector<int>(npnt, 0.);
    const auto corner = [&edges] () {
      std::vector<int> corner(edges.size(), -1);
      for (int iv = 0; iv < edges.size(); ++iv)
        corner[iv] = edges[iv].size() - 1;

      return corner;
    }();

    return sqdistance(origin, corner);
  }();

  return std::pow(1. - std::pow(sqdistance(point1, point2) / maxsqd, 1.5), 3.);
}



auto single_pass_smooth(const Arrayhist &rdev_tr_h, const Arrayhist &rdev_tr_l,
                        const Arrayhist &eval_n, const Arrayhist &eval_h, const Arrayhist &eval_l,
                        const std::vector<int> &bandwidth_h, const std::vector<int> &bandwidth_l, bool symmetrize,
                        const std::vector<std::vector<double>> &fine_edges, const std::vector<std::vector<double>> &coarse_edges,
                        double &chi2_h, double &chi2_l)
{
  const int nvar = coarse_edges.size(), nbin = eval_n.size();
  thread_local Workspace wsp(nbin, nvar);
  /*
  for (int ib = 0; ib < eval_n.size(); ++ib)
    std::cout << ib << " " << eval_n(ib, 0) << " " << eval_n(ib, 1) << " " << eval_h(ib, 0) << " " << eval_h(ib, 1) << " " << eval_l(ib, 0) << " " << eval_l(ib, 1) << "\n";
  */
  auto f_best_scale_num_den = [] (const Arrayhist &nominal, const Arrayhist &varied, const Arrayhist &rdev,
                                  const std::vector<std::vector<double>> &fine, const std::vector<std::vector<double>> &coarse) {
    auto numden = std::make_pair(0., 0.);
    auto &num = numden.first;
    auto &den = numden.second;

    const auto csmooth = rebin(apply_deviation(nominal, rdev), fine, coarse);
    const auto cnominal = rebin(nominal, fine, coarse);
    const auto cvaried = rebin(varied, fine, coarse);

    for (int ibin = 0; ibin < cnominal.size(); ++ibin) {
      double smnv = csmooth(ibin, 0) - cnominal(ibin, 0);
      double vmnv = cvaried(ibin, 0) - cnominal(ibin, 0);
      double vpnu = cvaried(ibin, 1) + cnominal(ibin, 1);

      if (vpnu == 0.)
        continue;

      num += smnv * vmnv / vpnu;
      den += smnv * smnv / vpnu;
    }

    return numden;
  };

  auto f_chi2 = [] (const Arrayhist &nominal, const Arrayhist &varied, const Arrayhist &rdev,
                    const std::vector<std::vector<double>> &fine, const std::vector<std::vector<double>> &coarse, double scale) {
    double sum = 0.;
    const auto csmooth = rebin(apply_deviation(nominal, rdev, scale), fine, coarse);
    const auto cnominal = rebin(nominal, fine, coarse);
    const auto cvaried = rebin(varied, fine, coarse);

    for (int ibin = 0; ibin < cnominal.size(); ++ibin) {
      double smvv = csmooth(ibin, 0) - cvaried(ibin, 0);
      double vpnu = cvaried(ibin, 1) + cnominal(ibin, 1);

      if (vpnu == 0.)
        continue;

      sum += smvv * smvv / vpnu;
    }

    return sum;
  };

  const bool samebw = (bandwidth_h == bandwidth_l);
  std::vector<int> icenter(nvar, -1), idxn(nvar, -1);

  std::vector<int> plane_h = {}, plane_l = {};
  std::vector<std::array<double, 2>> smooth_h = {}, smooth_l = {};
  plane_h.reserve(nbin);
  smooth_h.reserve(nbin);

  if (!samebw) {
    plane_l.reserve(nbin);
    smooth_l.reserve(nbin);
  }

  auto rdev_sm_h = Arrayhist(nbin), rdev_sm_l = Arrayhist(nbin);

  auto sums_h = std::vector<double>(nbin, 0.), sums_l = std::vector<double>(nbin, 0.);
  auto sum2s_h = sums_h, sum2s_l = sums_l;

  auto f_update_rdev = [] (const std::vector<int> &plane, const std::vector<std::array<double, 2>> &smooth,
                           Arrayhist &rdev, std::vector<double> &sums, std::vector<double> &sum2s,
                           const std::vector<std::vector<double>> &edges) {
    sums[plane[0]] += 1.;
    sum2s[plane[0]] += 1.;

    rdev(plane[0], 0) += smooth[0][0];
    rdev(plane[0], 1) += smooth[0][1];

    std::vector<int> idxn1(edges.size(), -1), idxn2(edges.size(), -1);
    for (int iv = 0; iv < edges.size(); ++iv)
      idxn1[iv] = index_1n(plane[0], iv, edges);

    for (int ipl = 1; ipl < plane.size(); ++ipl) {
      for (int iv = 0; iv < edges.size(); ++iv)
        idxn2[iv] = index_1n(plane[ipl], iv, edges);

      double wgt = weight(idxn1, idxn2, edges);
      sums[plane[ipl]] += wgt;
      sum2s[plane[ipl]] += wgt * wgt;

      rdev(plane[ipl], 0) += wgt * smooth[ipl][0];
      rdev(plane[ipl], 1) += wgt * wgt * smooth[ipl][1];
    }
  };

  for (int ibin = 0; ibin < nbin; ++ibin) {
    // get all the points within bandwidth centered around current bin
    for (int iv = 0; iv < nvar; ++iv)
      icenter[iv] = index_1n(ibin, iv, fine_edges);

    plane_h.emplace_back(ibin);

    if (!samebw)
      plane_l.emplace_back(ibin);

    for (int ib = 0; ib < nbin; ++ib) {
      if (ib == ibin)
        continue;

      bool withinbw_h = true, withinbw_l = true;
      for (int iv = 0; iv < nvar; ++iv) {
        idxn[iv] = index_1n(ib, iv, fine_edges);

        withinbw_h = withinbw_h and std::abs(icenter[iv] - idxn[iv]) <= bandwidth_h[iv];
        withinbw_l = !samebw and withinbw_l and std::abs(icenter[iv] - idxn[iv]) <= bandwidth_l[iv];
      }
      if (samebw)
        withinbw_l = withinbw_h;

      if (withinbw_h)
        plane_h.emplace_back(ib);

      if (!samebw and withinbw_l)
        plane_l.emplace_back(ib);
    }

    // dump into list of smooth planes
    smooth_h = fit_plane(plane_h, rdev_tr_h, fine_edges, wsp);

    // build the smooth relative deviation as the weighted average of the smooth planes
    f_update_rdev(plane_h, smooth_h, rdev_sm_h, sums_h, sum2s_h, fine_edges);

    if (!(samebw and symmetrize)) {
      smooth_l = fit_plane(plane_l, rdev_tr_l, fine_edges, wsp);
      f_update_rdev(plane_l, smooth_l, rdev_sm_l, sums_l, sum2s_l, fine_edges);
    }
    else
      f_update_rdev(plane_h, smooth_h, rdev_sm_l, sums_l, sum2s_l, fine_edges);

    plane_h.clear();
    plane_l.clear();

    smooth_h.clear();
    smooth_l.clear();
  }

  // we were just summing before, so now divide
  for (int ibin = 0; ibin < nbin; ++ibin) {
    if (sums_h[ibin] != 0.) {
      rdev_sm_h(ibin, 0) /= sums_h[ibin];
      rdev_sm_h(ibin, 1) /= sum2s_h[ibin];
    }
    else {
      rdev_sm_h(ibin, 0) = -1.; // such that apply_deviation in this bin gives 0
      rdev_sm_h(ibin, 1) = 0.;
    }

    if (sums_l[ibin] != 0.) {
      rdev_sm_l(ibin, 0) /= sums_l[ibin];
      rdev_sm_l(ibin, 1) /= sum2s_l[ibin];
    }
    else {
      rdev_sm_l(ibin, 0) = -1.;
      rdev_sm_l(ibin, 1) = 0.;
    }
  }

  // determine chi2 to the original, using the coarse binning, and dump into chi2 list
  // function is just analytical expression off AN-18-077
  const auto [snum_h, sden_h] = f_best_scale_num_den(eval_n, eval_h, rdev_sm_h, fine_edges, coarse_edges);
  const auto [snum_l, sden_l] = f_best_scale_num_den(eval_n, eval_l, rdev_sm_l, fine_edges, coarse_edges);

  double scale_h = snum_h / sden_h;
  double scale_l = snum_l / sden_l;

  const double variance_h = 1. / std::sqrt(sden_h);
  const double variance_l = 1. / std::sqrt(sden_l);

  //std::cout << "check symmetry scale " << scale_h << " " << scale_l << "\n";
  //std::cout << "check symmetry variance hi - lo " << variance_h - variance_l << "\n\n";

  if (symmetrize and samebw and (scale_h + scale_l) * (scale_h + scale_l) < variance_h + variance_l) {
    scale_h = 0.5 * (scale_h - scale_l);
    scale_l = -scale_h;
  }

  chi2_h += f_chi2(eval_n, eval_h, rdev_sm_h, fine_edges, coarse_edges, scale_h);
  chi2_l += f_chi2(eval_n, eval_l, rdev_sm_l, fine_edges, coarse_edges, scale_l);
  /*
  std::cout << "num " << snum_h << " " << snum_l << "\n";
  std::cout << "den " << sden_h << " " << sden_l << "\n";
  std::cout << "chi2 sum " << chi2_h << " " << chi2_l << "\n";
  std::cout << "scale " << scale_h << " " << scale_l << "\n";
  std::cout << "variance " << variance_h << " " << variance_l << "\n" << std::endl;
  */
  return std::make_pair(std::move(rdev_sm_h.scale(scale_h)), std::move(rdev_sm_l.scale(scale_l)));
}



/// perform the actual smoothing
/// also cross-validation over pre-determined set of bandwidths if there are more than 1 partitions
/// bandwidths is when one wants to run with a fixed bandwidth choice
/// return value is a bunch of histograms to be saved? maybe plot them?
auto smooth_templates_impl(Framework::Dataset<TChain> &data_n,
                           Framework::Dataset<TChain> &data_h,
                           Framework::Dataset<TChain> &data_l,
                           Framework::Collection<float, double> &coll,
                           const std::tuple<
                           std::vector<std::vector<std::string>>,
                           std::vector<std::vector<double>>,
                           std::vector<std::vector<double>>,
                           std::string> &variables_bins,
                           const std::string &systematic,
                           bool symmetrize,
                           bool binomial,
                           int npartition,
                           int nrepeatcv,
                           const std::vector<double> &fixed_bandwidth = {})
{
  const int npart = npartition * nrepeatcv;
  const bool runcv = npartition > 1 and npart > 1;
  std::vector<std::unique_ptr<TH1>> result;
  if (!fixed_bandwidth.empty() and runcv) {
    std::cerr << "Fixed bandwidth mode can only be ran with 1 partitions/bootstrapping i.e. without cross validation" << std::endl;
    return result;
  }

  const auto &coarse_edges = std::get<1>(variables_bins);
  const auto &fine_edges = std::get<2>(variables_bins);
  const int nvar = coarse_edges.size(), nfbin = count_nbin(fine_edges);

  if (!fixed_bandwidth.empty() and fixed_bandwidth.size() != nvar) {
    std::cerr << "List of bandwidth must be as long as the list of variables in the fixed bandwidth mode" << std::endl;
    return result;
  }

  // prepare helpers
  auto f_make_rdev = [] (const Arrayhist &nominal, const Arrayhist &vary_h, const Arrayhist &vary_l, bool symmetrize) {
    auto rdevs = std::make_tuple(relative_deviation(nominal, vary_h), relative_deviation(nominal, vary_l));
    if (symmetrize) {
      std::get<0>(rdevs) = symmetrize_variation(std::get<0>(rdevs), std::get<1>(rdevs));

      for (int ibin = 0; ibin < nominal.size(); ++ibin) {
        std::get<1>(rdevs)(ibin, 0) = -std::get<0>(rdevs)(ibin, 0);
        std::get<1>(rdevs)(ibin, 1) = std::get<0>(rdevs)(ibin, 1);
      }
    }

    return rdevs;
  };

  std::cout << "Mode smooth: making nominal histogram..." << std::endl;
  const auto vcn = count_and_bin(data_n, coll, variables_bins, 1, 1, true, true);
  const auto vfn = count_and_bin(data_n, coll, variables_bins, 1, 2, true, false);

  std::cout << "Mode smooth: making systematic " << systematic << " histograms..." << std::endl;
  const auto vch = count_and_bin(data_h, coll, variables_bins, 1, 1, true, true);
  const auto vfh = count_and_bin(data_h, coll, variables_bins, 1, 2, true, false);

  const auto vcl = count_and_bin(data_l, coll, variables_bins, 1, 1, true, true);
  const auto vfl = count_and_bin(data_l, coll, variables_bins, 1, 2, true, false);

  const auto &coarse_n = vcn[0], &coarse_h = vch[0], &coarse_l = vcl[0];
  const auto &fine_n = vfn[0], &fine_h = vfh[0], &fine_l = vfl[0];
  /*
  for (int ib = 0; ib < fine_n.size(); ++ib)
    std::cout << ib << " " << fine_n(ib, 0) << " " << fine_n(ib, 1) << " " << fine_h(ib, 0) << " " << fine_h(ib, 1) << " " << fine_l(ib, 0) << " " << fine_l(ib, 1) << "\n";
  throw std::runtime_error( "I quit my job, bye..." );
  */
  std::cout << "Translating bandwidths..." << std::endl;
  const auto bandwidths = bandwidth_to_test(fine_edges, fixed_bandwidth);
  std::vector<double> chi2s_h(bandwidths.size(), 0.), chi2s_l(bandwidths.size(), 0.);
  std::vector<std::thread> threads(bandwidths.size());

  // a couple stuff for binomial mode
  // ptr so as to not consume resources if not needed
  std::unique_ptr<std::mt19937_64> rng;
  std::unique_ptr<std::vector<Arrayhist>> unweighted_n, unweighted_h, unweighted_l;
  std::unique_ptr<Arrayhist> average_n, average_h, average_l;
  std::unique_ptr<std::binomial_distribution<int>> dist;
  std::unique_ptr<std::binomial_distribution<int>::param_type[]> params_n, params_h, params_l;
  const double invnpart = 1. / npartition;

  if (runcv) {
    const int report_every = [] (int npart) {
      npart = std::pow(10., std::round(std::log10(npart)) - 1);
      return npart > 1 ? npart : 1;
    }(npart);

    if (binomial) {
      rng = std::make_unique<std::mt19937_64>(random_generator<>());

      unweighted_n = std::make_unique<std::vector<Arrayhist>>(count_and_bin(data_n, coll, variables_bins, 1, 2, false, false));
      unweighted_h = std::make_unique<std::vector<Arrayhist>>(count_and_bin(data_h, coll, variables_bins, 1, 2, false, false));
      unweighted_l = std::make_unique<std::vector<Arrayhist>>(count_and_bin(data_l, coll, variables_bins, 1, 2, false, false));

      average_n = std::make_unique<Arrayhist>(divide(fine_n, (*unweighted_n)[0]));
      average_h = std::make_unique<Arrayhist>(divide(fine_h, (*unweighted_h)[0]));
      average_l = std::make_unique<Arrayhist>(divide(fine_l, (*unweighted_l)[0]));

      dist = std::make_unique<std::binomial_distribution<int>>();

      params_n = std::make_unique<std::binomial_distribution<int>::param_type[]>(nfbin);
      params_h = std::make_unique<std::binomial_distribution<int>::param_type[]>(nfbin);
      params_l = std::make_unique<std::binomial_distribution<int>::param_type[]>(nfbin);
      for (int ibin = 0; ibin < nfbin; ++ibin) {
        params_n[ibin] = std::binomial_distribution<int>::param_type((*unweighted_n)[0](ibin, 0), invnpart);
        params_h[ibin] = std::binomial_distribution<int>::param_type((*unweighted_h)[0](ibin, 0), invnpart);
        params_l[ibin] = std::binomial_distribution<int>::param_type((*unweighted_l)[0](ibin, 0), invnpart);
      }
    }

    std::vector<Arrayhist> trev_n, trev_h, trev_l;
    trev_n.reserve(npartition);
    trev_h.reserve(npartition);
    trev_l.reserve(npartition);
    for (int ipart = 0; ipart < npartition; ++ipart) {
      trev_n.emplace_back(Arrayhist(nfbin));
      trev_h.emplace_back(Arrayhist(nfbin));
      trev_l.emplace_back(Arrayhist(nfbin));
    }

    for (int ipart = 0; ipart < npart; ++ipart) {
      if (ipart % report_every == 0)
        std::cout << "Starting cross validation run on partition " << ipart + 1 << "/" << npart <<
          " with " << bandwidths.size() << " bandwidth choices..." << std::endl;
      const int iusepart = ipart % npartition;

      if (binomial) {
        for (int ibin = 0; ibin < nfbin; ++ibin) {
          dist->param(params_n[ibin]);
          const auto ibn = (*dist)(*rng);
          trev_n[iusepart](ibin, 0) = (*average_n)(ibin, 0) * ibn;
          trev_n[iusepart](ibin, 1) = (*average_n)(ibin, 1) * ibn * ibn;

          dist->param(params_h[ibin]);
          const auto ibh = (*dist)(*rng);
          trev_h[iusepart](ibin, 0) = (*average_h)(ibin, 0) * ibh;
          trev_h[iusepart](ibin, 1) = (*average_h)(ibin, 1) * ibh * ibh;

          dist->param(params_l[ibin]);
          const auto ibl = (*dist)(*rng);
          trev_l[iusepart](ibin, 0) = (*average_l)(ibin, 0) * ibl;
          trev_l[iusepart](ibin, 1) = (*average_l)(ibin, 1) * ibl * ibl;
        }
      }
      else if (iusepart == 0) {
        trev_n = count_and_bin(data_n, coll, variables_bins, npartition, 2, true, false);
        trev_h = count_and_bin(data_h, coll, variables_bins, npartition, 2, true, false);
        trev_l = count_and_bin(data_l, coll, variables_bins, npartition, 2, true, false);
      }

      const auto &eval_n = trev_n[iusepart], &eval_h = trev_h[iusepart], &eval_l = trev_l[iusepart];
      const auto train_n = sum(fine_n, eval_n, 1., -1.), train_h = sum(fine_h, eval_h, 1., -1.), train_l = sum(fine_l, eval_l, 1., -1.);
      const auto [rdev_tr_h, rdev_tr_l] = f_make_rdev(train_n, train_h, train_l, symmetrize);

      for (int ibw = 0; ibw < bandwidths.size(); ++ibw) {
        threads[ibw] = std::thread(single_pass_smooth,
                                   std::cref(rdev_tr_h), std::cref(rdev_tr_l),
                                   std::cref(eval_n), std::cref(eval_h), std::cref(eval_l),
                                   std::cref(bandwidths[ibw]), std::cref(bandwidths[ibw]), symmetrize,
                                   std::cref(fine_edges), std::cref(coarse_edges), std::ref(chi2s_h[ibw]), std::ref(chi2s_l[ibw]));
      }

      for (auto &thread : threads)
        thread.join();
    }
  }

  // pick either minimal sum of or individual chi2 based on symmetric fit or na
  int imin_h = -1, imin_l = -1;
  if (symmetrize) {
    double minsum = std::numeric_limits<double>::max();
    for (int ibw = 0; ibw < bandwidths.size(); ++ibw) {
      if (chi2s_h[ibw] + chi2s_l[ibw] < minsum) {
        minsum = chi2s_h[ibw] + chi2s_l[ibw];
        imin_h = ibw;
      }
    }
    imin_l = imin_h;
  }
  else {
    imin_h = std::distance(std::begin(chi2s_h), std::min_element(std::begin(chi2s_h), std::end(chi2s_h)));
    imin_l = std::distance(std::begin(chi2s_l), std::min_element(std::begin(chi2s_l), std::end(chi2s_l)));
  }

  if (runcv)
    std::cout << "\nCross validation complete. Average chi2 of the best bandwidth is " << chi2s_h[imin_h] / npart << 
      " (up) and " << chi2s_l[imin_l] / npart << " (down). Starting full smoothing with chosen bandwidths..." << std::endl;

  std::cout << "Chosen bandwidths are (in the same order as given variables, rounded to nearest discrete representation):" << "\n";
  for (int iv = 0; iv < nvar; ++iv)
    std::cout << double(bandwidths[imin_h][iv]) / (fine_edges[iv].size() - 1) << " ";
  std::cout << "(up)\n";
  for (int iv = 0; iv < nvar; ++iv)
    std::cout << double(bandwidths[imin_l][iv]) / (fine_edges[iv].size() - 1) << " ";
  std::cout << "(down)\n" << std::endl;
  // FIXME any plotting of chi2 for each bandwidth etc

  auto [rdev_h, rdev_l] = f_make_rdev(fine_n, fine_h, fine_l, symmetrize);
  double chi2_h = 0., chi2_l = 0.;
  const auto [smooth_rdev_h, smooth_rdev_l] = single_pass_smooth(rdev_h, rdev_l, fine_n, fine_h, fine_l, bandwidths[imin_h], bandwidths[imin_l],
                                                                 symmetrize, fine_edges, coarse_edges, chi2_h, chi2_l);
  std::cout << "Full smoothing complete. Chi2 values are " << chi2_h << 
    " (up) and " << chi2_l << " (down). Preparing output to be saved..." << std::endl;

  auto f_concatenate = [] (std::vector<std::unique_ptr<TH1>> &to, std::vector<std::unique_ptr<TH1>> &&from) {
    to.insert(std::end(to), std::make_move_iterator(std::begin(from)), std::make_move_iterator(std::end(from)));
  };

  f_concatenate(result, array_to_root(variables_bins, "nominal_source_template", coarse_edges, coarse_n));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_source_template", coarse_edges, coarse_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_source_template", coarse_edges, coarse_l));

  const auto coarse_smooth_h = rebin(apply_deviation(fine_n, smooth_rdev_h), fine_edges, coarse_edges);
  const auto coarse_smooth_l = rebin(apply_deviation(fine_n, smooth_rdev_l), fine_edges, coarse_edges);
  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_smooth_template", coarse_edges, coarse_smooth_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_smooth_template", coarse_edges, coarse_smooth_l));

  if (symmetrize)
    std::tie(rdev_h, rdev_l) = f_make_rdev(fine_n, fine_h, fine_l, false);
  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_fine_source_deviation", fine_edges, rdev_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_fine_source_deviation", fine_edges, rdev_l));

  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_fine_smooth_deviation", fine_edges, smooth_rdev_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_fine_smooth_deviation", fine_edges, smooth_rdev_l));

  const auto [crdev_h, crdev_l] = f_make_rdev(coarse_n, coarse_h, coarse_l, false);
  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_source_deviation", coarse_edges, crdev_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_source_deviation", coarse_edges, crdev_l));

  const auto smooth_crdev_h = divide(sum(coarse_smooth_h, coarse_n, 1., -1., false, true), coarse_n, false, true);
  const auto smooth_crdev_l = divide(sum(coarse_smooth_l, coarse_n, 1., -1., false, true), coarse_n, false, true);
  f_concatenate(result, array_to_root(variables_bins, systematic + "_up_smooth_deviation", coarse_edges, smooth_crdev_h));
  f_concatenate(result, array_to_root(variables_bins, systematic + "_down_smooth_deviation", coarse_edges, smooth_crdev_l));

  return result;
}



auto smooth_templates(const std::vector<std::string> &files,
                      const std::string &tree,
                      const std::tuple<
                      std::vector<std::vector<std::string>>,
                      std::vector<std::vector<double>>,
                      std::vector<std::vector<double>>,
                      std::string> &variables_bins,
                      const std::string &systematic,
                      bool symmetrize,
                      bool binomial,
                      int npartition,
                      int nrepeatcv,
                      const std::vector<double> &fixed_bandwidth = {})
{
  Framework::Dataset<TChain> data_n("data_n", tree);
  data_n.set_files(files);

  if (!fixed_bandwidth.empty()) {
    std::cout << "Fixed bandwidth mode detected. Ignoring --npartition and --nrepeatcv options..." << std::endl;
    npartition = 1;
    nrepeatcv = 1;
  }

  auto fvary_h = files;
  for (auto &file : fvary_h)
    replace(file, ".root", std::string("_") + systematic + "_up.root");
  Framework::Dataset<TChain> data_h("data_h", tree);
  data_h.set_files(fvary_h);

  auto fvary_l = fvary_h;
  for (auto &file : fvary_l)
    replace(file, "_up.root", "_down.root");
  Framework::Dataset<TChain> data_l("data_l", tree);
  data_l.set_files(fvary_l);

  auto coll = make_collection(data_n, variables_bins);

  // FIXME implement equiprobable binning and bin map editing here?

  return smooth_templates_impl(data_n, data_h, data_l, coll, variables_bins, systematic, symmetrize, binomial, npartition, nrepeatcv, fixed_bandwidth);
}

// FIXME to do:
// better diagnostics reporting
// more complete multithreading support (currently only smoothing, histogramming is single-thread. requires thread-safe rng)
// smoothing in equiprobable binning in addition to equidistant

#endif