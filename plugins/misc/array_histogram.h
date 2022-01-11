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
  /// constructor
  Arrayhist() : nbin(0) {}

  Arrayhist(int nbin_) : nbin(nbin_) {
    if (nbin > 0)
      hist = std::make_unique<double[]>(2 * nbin);
    else
      nbin = 0;
  }

  Arrayhist(const Arrayhist &ah) : nbin(ah.nbin) {
    for (int ibin = 0; ibin < nbin; ++ibin) {
      hist[ibin] = ah(ibin, 0);
      hist[ibin + nbin] = ah(ibin, 1);
    }
  }

  Arrayhist(Arrayhist &&ah) : nbin(ah.nbin) { std::swap(hist, ah.hist); }

  /// assignment
  Arrayhist& operator=(const Arrayhist &ah) {
    nbin = ah.nbin;
    for (int ibin = 0; ibin < nbin; ++ibin) {
      hist[ibin] = ah(ibin, 0);
      hist[ibin + nbin] = ah(ibin, 1);
    }

    return *this;
  }

  Arrayhist& operator=(Arrayhist &&ah) = default;

  double& operator()(int ibin, bool variance) { return hist[ibin + (variance * nbin)]; }
  const double& operator()(int ibin, bool variance) const { return hist[ibin + (variance * nbin)]; }
  int size() const { return nbin; }

  void initialize(int nbin_) {
    if (nbin != 0)
      return;

    nbin = nbin_;
    hist = std::make_unique<double[]>(2 * nbin);
  }

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



/// as above, but over all dimensions in a go
std::vector<int> index_1n(int idx1, const std::vector<std::vector<double>> &edges)
{
  std::vector<int> idxn(edges.size(), -1);
  for (int iv = 0; iv < edges.size(); ++iv)
    idxn[iv] = index_1n(idx1, iv, edges);

  return idxn;
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



std::vector<std::vector<double>> centers_of(const std::vector<std::vector<double>> &edges)
{
  std::vector<std::vector<double>> centers(edges.size());
  for (int iv = 0; iv < edges.size(); ++iv) {
    std::vector<double> center(edges[iv].size() - 1, std::numeric_limits<double>::max());
    for (int ib = 0; ib < edges[iv].size() - 1; ++ib)
      center[ib] = 0.5 * (edges[iv][ib] + edges[iv][ib + 1]);

    centers.emplace_back(center);
  }

  return centers;
}



std::vector<std::vector<double>> edges_of(const std::vector<std::vector<double>> &centers)
{
  std::vector<std::vector<double>> edges(centers.size());
  for (int iv = 0; iv < centers.size(); ++iv) {
    std::vector<double> edge(centers[iv].size() + 1, std::numeric_limits<double>::max());
    edge[0] = centers[iv][0] - (0.5 * (centers[iv][1] - centers[iv][0]));

    for (int ib = 1; ib < centers[iv].size(); ++ib)
      edge[ib] = 0.5 * (centers[iv][ib] + centers[iv][ib - 1]);

    edge.back() = centers[iv].back() + (0.5 * (centers[iv].back() - centers[iv][centers[iv].size() - 2]));
    edges.emplace_back(edge);
  }

  return edges;
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
    result(ibin, 1) += (exact2) ? 0. : h2(ibin, 1) * h1(ibin, 0) * h1(ibin, 0);
  }

  return result;
}



auto divide(const Arrayhist &h1, const Arrayhist &h2, bool exact1 = false, bool exact2 = false)
{
  auto result = Arrayhist(h1.size());

  for (int ibin = 0; ibin < h1.size(); ++ibin) {
    result(ibin, 0) = (h2(ibin, 0) != 0.) ? h1(ibin, 0) / h2(ibin, 0) : 0.;
    result(ibin, 1) = (exact1 or h2(ibin, 0) == 0.) ? 0. : h1(ibin, 1) / h2(ibin, 0) / h2(ibin, 0);
    result(ibin, 1) += (exact2 or h2(ibin, 0) == 0.) ? 0. : h2(ibin, 1) * result(ibin, 0) * result(ibin, 0) / h2(ibin, 0) / h2(ibin, 0);
  }

  return result;
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
