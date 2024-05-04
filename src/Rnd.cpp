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

#include <algorithm>
#include <math.h>
#include <stdlib.h>

#include "asyncde/Rnd.h"

double asyncde::Rnd::BoxMuller(double mu, double sigma, double xmin,
                               double xmax) {
  double normal_rnd;
  double u1, u2;

  do {
    u1 = unidouble(engine);
    u2 = unidouble(engine);
    normal_rnd = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2) * sigma + mu;
  } while (normal_rnd < xmin || normal_rnd > xmax);

  return normal_rnd;
}

double asyncde::Rnd::randCauchy(const double sigma) {
  double urand;
  do
    urand = unidouble(engine);
  while (urand < std::numeric_limits<double>::epsilon() ||
         urand > 1.0 - std::numeric_limits<double>::epsilon());

  return sigma * tan(M_PI * (urand - 0.5));
}

double asyncde::Rnd::randCauchyTruncated(double mu, double sigma, double xmin,
                                         double xmax) {
  double umin = atan((xmin - mu) / sigma) / M_PI + 0.5;
  double umax = atan((xmax - mu) / sigma) / M_PI + 0.5;
  double urand = umin + unidouble(engine) * (umax - umin);

  return mu + sigma * tan(M_PI * (urand - 0.5));
}

int asyncde::Rnd::randdistinct(unsigned int ngen, unsigned int nmax,
                               std::vector<unsigned int> &nrand,
                               std::vector<unsigned int> &nrand_sorted) {
  if (nmax < ngen)
    return -1;

  nrand.resize(ngen);
  nrand_sorted.resize(ngen);
  for (unsigned int irand = 0; irand < ngen; irand++) {
    unsigned int urand = uniuint(engine) % (nmax - irand);
    unsigned int j = 0;
    for (; j < irand; j++)
      if (urand >= nrand_sorted[j])
        urand++;
      else
        break;
    nrand[irand] = urand;
    std::move_backward(&nrand_sorted[j], &nrand_sorted[irand],
                       &nrand_sorted[irand + 1]);
    nrand_sorted[j] = urand;
  }

  return 0;
}
