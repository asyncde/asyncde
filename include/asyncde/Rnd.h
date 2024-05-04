/**
 This file is a part of the AsyncDE library.

 If you are using AsyncDE as part of your research, teaching,
 or other activities, we would be grateful if you could cite our work:
 Zhabitskaya, E., Zhabitsky, M. (2013).
 Asynchronous Differential Evolution with Restart.
 In: Dimov, I., Faragó, I., Vulkov, L. (eds) Numerical Analysis and Its
 Applications. NAA 2012. Lecture Notes in Computer Science, vol 8236. Springer,
 Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-41515-9_64

 The AsyncDE library is free software.
 You can redistribute it and/or modify it under the terms
 of the GNU Lesser General Public License as published
 by the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 see https://www.gnu.org/licenses/.
*/

#ifndef ASYNCDE_RND_H
#define ASYNCDE_RND_H

#include <climits>
#include <random>
#include <vector>

namespace asyncde {

class Rnd {
  //  protected:
public:
  std::mt19937 engine;
  std::uniform_real_distribution<double> unidouble;
  std::uniform_int_distribution<unsigned int> uniuint;
  std::binomial_distribution<unsigned int> binomialuint;

public:
  Rnd(unsigned long int seed) {
    engine.seed(seed);
    unidouble = std::uniform_real_distribution<double>{0.0, 1.0};
    uniuint = std::uniform_int_distribution<unsigned int>{0, UINT_MAX};
  }

  double next(double xmin = 0.0, double xmax = 1.0) {
    return xmin + unidouble(engine) * (xmax - xmin);
  }

  // returns random long int in [0; umax - 1]
  unsigned int next_uniuint(unsigned int umax) {
    return uniuint(engine) % umax;
  }

  // returns random int in [0; upbound]
  unsigned int next_binomialuint(unsigned int upbound, double psucc) {
    binomialuint.param(
        std::binomial_distribution<unsigned int>::param_type(upbound, psucc));
    return binomialuint(engine);
  }

  /*
   * Box-Muller algorithm to generate random number
   * according to truncated [xmin, xmax] Normal(mu, sigma) distribution
   */
  double BoxMuller(double mu, double sigma, double xmin, double xmax);

  // pseudorandom number according to Cauchy distribution
  double randCauchy(const double sigma);

  // pseudorandom number according to truncated Cauchy distribution
  double randCauchyTruncated(double mu, double sigma, double xmin, double xmax);

  // generate ngen distinct pseudorandom numbers in the range [0; nmax)
  int randdistinct(unsigned int ngen, unsigned int nmax,
                   std::vector<unsigned int> &nrand,
                   std::vector<unsigned int> &nrand_sorted);
};

} // namespace asyncde

#endif
