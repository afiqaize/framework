// -*- C++ -*-
// author: afiq anuar
// short: functions pertaining to manipulating histograms represented as arrays
// note: the histogram representation assumes no under- and overflows 
// note: i.e. they're added to first and last bins in the relevant axis respectively

#ifndef FWK_ARRAY_HISTOGRAM_H
#define FWK_ARRAY_HISTOGRAM_H

#include <memory>
#include <vector>
#include <limits>

class Arrayhist {
public:
  Arrayhist() = delete;
  Arrayhist(int nbin_) : nbin(nbin_), hist(std::make_unique<double[]>(2 * nbin)) {}
  /*/ copy ctor, commented because not needed, but kept for reference
  Arrayhist(const Arrayhist &ah) : nbin(ah.nbin) {
    for (int ibin = 0; ibin < nbin; ++ibin) {
      hist[ibin] = ah(ibin, 0);
      hist[ibin + nbin] = ah(ibin, 1);
    }
  }
  */
  Arrayhist(Arrayhist &&ah) : nbin(ah.nbin) { std::swap(hist, ah.hist); }
  /*/ copy assign, commented because not needed, but kept for reference
  Arrayhist& operator=(const Arrayhist& ah) {
    nbin = ah.nbin;
    for (int ibin = 0; ibin < nbin; ++ibin) {
      hist[ibin] = ah(ibin, 0);
      hist[ibin + nbin] = ah(ibin, 1);
    }

    return *this;
  }
  */
  Arrayhist& operator=(Arrayhist&& ah) {
    nbin = ah.nbin;
    std::swap(hist, ah.hist);
    return *this;
  }

  double& operator()(int ibin, int variance) { return hist[ibin + (variance * nbin)]; }
  const double& operator()(int ibin, int variance) const { return hist[ibin + (variance * nbin)]; }
  int size() const { return nbin; }

  double sumw() const {
    double sum = 0.;
    for (int ibin = 0; ibin < nbin; ++ibin)
      sum += hist[ibin];
    return sum;
  }

  Arrayhist& scale(double sc) {
    for (int ibin = 0; ibin < nbin; ++ibin) {
      hist[ibin] *= sc;
      hist[ibin + nbin] *= sc * sc;
    }

    return *this;
  }

  Arrayhist& shift(double sh) {
    for (int ibin = 0; ibin < nbin; ++ibin)
      hist[ibin] += sh;

    return *this;
  }

private:
  int nbin;
  std::unique_ptr<double[]> hist;
};



/// translate nD bin index to unrolled 1D
int index_n1(const std::vector<int> &idxn, const std::vector<std::vector<double>> &edges)
{
  int idx1 = idxn[0];
  for (int ii = 1; ii < idxn.size(); ++ii) {
    int multiplier = 1;
    for (int jj = ii - 1; jj > -1; --jj)
      multiplier *= edges[jj].size() - 1;

    idx1 += idxn[ii] * multiplier;
  }

  return idx1;
}



/// and the inverse operation
int index_1n(int idx1, int dim, const std::vector<std::vector<double>> &edges)
{
  int multiplier = 1;
  for (int dd = dim - 1; dd > -1; --dd)
    multiplier *= edges[dd].size() - 1;

  return (idx1 / multiplier) % (edges[dim].size() - 1);
}



std::vector<double> center(int idx1, const std::vector<std::vector<double>> &edges)
{
  std::vector<double> result(edges.size(), std::numeric_limits<double>::max());
  for (int iv = 0; iv < edges.size(); ++iv) {
    int idxn = index_1n(idx1, iv, edges);
    result[iv] = 0.5 * (edges[iv][idxn] + edges[iv][idxn + 1]);
  }

  return result;
}



int count_nbin(const std::vector<std::vector<double>> &edges)
{
  auto fnbin = [] (int nbin, const auto &edge) {
    return nbin * (edge.size() - 1);
  };

  return std::accumulate(std::begin(edges), std::end(edges), 1, fnbin);
}



int find_bin(const std::vector<double> &values, const std::vector<std::vector<double>> &edges)
{
  std::vector<int> idxn(edges.size(), -999);
  for (int ivar = 0; ivar < edges.size(); ++ivar) {
    for (int ibin = 1; ibin < edges[ivar].size(); ++ibin) {
      if (values[ivar] < edges[ivar][ibin]) {
        idxn[ivar] = ibin - 1;
        break;
      }
    }

    if (idxn[ivar] == -999)
      idxn[ivar] = edges[ivar].size() - 2;
  }

  return index_n1(idxn, edges);
}



template <typename Number>
double sqdistance(const std::vector<Number> &point1, const std::vector<Number> &point2)
{
  double sqd = 0.;
  for (int ip = 0; ip < point1.size(); ++ip)
    sqd += (point1[ip] - point2[ip]) * (point1[ip] - point2[ip]);

  return sqd;
}



/// rebin a histogram with a fine binning (provided by its edges) to a coarse one
/// it is assumed without checks that the fine binning is a regular n-split of the coarse
auto rebin(const Arrayhist &histogram, const std::vector<std::vector<double>> &fine, const std::vector<std::vector<double>> &coarse)
{
  auto result = Arrayhist(count_nbin(coarse));
  for (int idxf = 0; idxf < histogram.size(); ++idxf) {
    auto idxc = find_bin(center(idxf, fine), coarse);

    result(idxc, 0) += histogram(idxf, 0);
    result(idxc, 1) += histogram(idxf, 1);
  }

  return result;
}




/// the booleans are for exact addition i.e. without adding up the errors too
auto sum(const Arrayhist &h1, const Arrayhist &h2, double c1 = 1., double c2 = 1., bool exact1 = false, bool exact2 = false)
{
  auto result = Arrayhist(h1.size());

  for (int ibin = 0; ibin < h1.size(); ++ibin) {
    result(ibin, 0) = (c1 * h1(ibin, 0)) + (c2 * h2(ibin, 0));
    result(ibin, 1) = (exact1) ? 0. : (c1 * c1 * h1(ibin, 1));
    result(ibin, 1) += (exact2) ? 0. : (c2 * c2 * h2(ibin, 1));
  }

  return result;
}



auto multiply(const Arrayhist &h1, const Arrayhist &h2, bool exact1 = false, bool exact2 = false)
{
  auto result = Arrayhist(h1.size());

  for (int ibin = 0; ibin < h1.size(); ++ibin) {
    result(ibin, 0) = h1(ibin, 0) * h2(ibin, 0);
    result(ibin, 1) = (exact1) ? 0. : h1(ibin, 1) * h2(ibin, 0) * h2(ibin, 0);
    result(ibin, 1) = (exact2 and h1(ibin, 0) != 0. and h2(ibin, 0) != 0.) ? result(ibin, 1) : 
      result(ibin, 0) * result(ibin, 0) * ((h1(ibin, 1) / h1(ibin, 0) / h1(ibin, 0)) + (h2(ibin, 1) / h2(ibin, 0) / h2(ibin, 0)));
  }

  return result;
}



auto divide(const Arrayhist &h1, const Arrayhist &h2, bool exact1 = false, bool exact2 = false)
{
  auto result = Arrayhist(h1.size());

  for (int ibin = 0; ibin < h1.size(); ++ibin) {
    result(ibin, 0) = (h2(ibin, 0) != 0.) ? h1(ibin, 0) / h2(ibin, 0) : 0.;
    result(ibin, 1) = (exact1 and h2(ibin, 0) != 0.) ? 0. : h1(ibin, 1) / h2(ibin, 0) / h2(ibin, 0);
    result(ibin, 1) = (exact2 and h1(ibin, 0) != 0. and h2(ibin, 0) != 0.) ? result(ibin, 1) : 
      result(ibin, 0) * result(ibin, 0) * ((h1(ibin, 1) / h1(ibin, 0) / h1(ibin, 0)) + (h2(ibin, 1) / h2(ibin, 0) / h2(ibin, 0)));
  }

  return result;
}



/// some converter methods for snapshotting
auto bandwidth_to_hist(const std::vector<std::vector<int>> &bandwidths, const std::vector<std::vector<double>> &edges, int idim)
{
  Arrayhist hist(bandwidths.size());
  for (int ibw = 0; ibw < bandwidths.size(); ++ibw) {
    hist(ibw, 0) = double(bandwidths[ibw][idim]) / (edges[idim].size() - 1);
    hist(ibw, 1) = 0.;
  }

  return hist;
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



auto vector_to_hist(const std::vector<double> &vec)
{
  Arrayhist hist(vec.size());
  for (int ibin = 0; ibin < vec.size(); ++ibin) {
    hist(ibin, 0) = vec[ibin];
    hist(ibin, 1) = 0.;
  }

  return hist;
}

#endif
