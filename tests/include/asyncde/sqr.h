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

#ifndef ASYNCDE_SQR_H
#define ASYNCDE_SQR_H

#include "asyncde/Function.h"

#include <stdio.h>

namespace asyncde {

/*
 Function of multidimensional real-valued arguments
 y[0] - objective function valued
 y[1]..y[y.size()-1] - summands in a case of sum of squares
*/

class sqr : public asyncde::Functor1D {
public:
  double fmin;

  sqr(double _fmin) : fmin(_fmin) {}

  virtual void operator()(const std::vector<double> &cX,
                          std::vector<double> &y) const override {
    y[0] = fmin;
    double d;

    if (y.size() > 1) {
      if (y.size() != cX.size() + 1) {
        fprintf(
            stderr,
            "sqr::operator() y.size()==%lu doesn't match (x.size()=%lu)+1\n",
            y.size(), cX.size());
        y[0] = fmin - 1.0;
        return;
      }

      std::vector<double>::iterator ity = y.begin() + 1;
      for (std::vector<double>::const_iterator itx = cX.begin();
           itx != cX.end(); itx++, ity++) {
        d = *itx - 1.5;
        *ity = d * d;
        y[0] += *ity;
      }
    } else
      for (double xcoord : cX) {
        d = xcoord - 1.5;
        y[0] += d * d;
      }
  }
};

} // namespace asyncde

#endif
