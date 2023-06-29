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

#ifndef ASYNCDE_COMPOSITE_ROSENBROCK_H
#define ASYNCDE_COMPOSITE_ROSENBROCK_H

#include "asyncde/Function.h"

#include <stdio.h>

namespace asyncde {

/*
 Function of multidimensional real-valued arguments
 y[0] - objective function valued
 y[1]..y[y.size()-1] - summands in a case of sum of squares
*/

//#define TEST_FIXED_VARS

#define SQR(x) ((x) * (x))

class composite_rosenbrock : public asyncde::Functor1D {
protected:
  int cluster_size = 2;

public:
  double fmin;

  composite_rosenbrock(double _fmin, int _cluster_size = 2)
      : cluster_size(_cluster_size), fmin(_fmin) {}

  virtual void operator()(const std::vector<double> &cX,
                          std::vector<double> &y) const override {
    y[0] = fmin;

#ifdef TEST_FIXED_VARS
    int nfreevars = cX.size() - 1;
#else
    int nfreevars = cX.size();
#endif
    int nclusters = nfreevars / cluster_size;

    if (1 == y.size()) {
      int i, imin, icluster;

      for (icluster = 0; icluster < nclusters; icluster++) {
        imin = icluster * cluster_size;
        for (i = imin; i < imin + cluster_size - 1; i++)
          y[0] += 100.0 * SQR(cX[i + 1] - SQR(cX[i])) + SQR(cX[i] - 1.0);
      }

      imin = nclusters * cluster_size;
      for (i = imin; i < nfreevars; i++)
        y[0] += SQR(cX[i] - 1.5);
    } else {
      if (y.size() != (unsigned int)(nclusters + nfreevars -
                                     nclusters * cluster_size + 1)) {
        fprintf(stderr,
                "composite_rosenbrock::operator() y.size()==%lu doesn't match "
                "%i (cluster_size=%i, nclusters=%i, nfreevars=%i)\n",
                y.size(), nclusters + nfreevars - nclusters * cluster_size + 1,
                cluster_size, nclusters, nfreevars);
        y[0] = fmin - 1.0;
        return;
      }

      int i, imin, icluster;
      for (icluster = 0; icluster < nclusters; icluster++) {
        y[icluster + 1] = 0.0;
        imin = icluster * cluster_size;
        for (i = imin; i < imin + cluster_size - 1; i++)
          y[icluster + 1] +=
              100.0 * SQR(cX[i + 1] - SQR(cX[i])) + SQR(cX[i] - 1.0);
        y[0] += y[icluster + 1];
      }

      imin = nclusters * cluster_size;
      for (i = imin; i < nfreevars; i++) {
        y[nclusters + i - imin] += SQR(cX[i] - 1.5);
        y[0] += y[nclusters + i - imin];
      }
    }
  }
};

} // namespace asyncde

#endif
