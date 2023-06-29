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

/**
 This file collects simple math utils for the asyncde software
 */

#ifndef ASYNCDE_MATHEXTRA_H
#define ASYNCDE_MATHEXTRA_H

#include <cmath>
#include <limits>

namespace asyncde {

static inline int fuzzy_cmp(double v1, double v2, double epsilon) {
  double absv1, absv2;

  if (std::isnan(v1)) {
    // treat nan as a +inf
    v1 = std::numeric_limits<double>::infinity();
    absv1 = v1;
  } else
    absv1 = std::fabs(v1);

  if (std::isnan(v2)) {
    // treat nan as a +inf
    v2 = std::numeric_limits<double>::infinity();
    absv2 = v2;
  } else
    absv2 = std::fabs(v2);

  double absmax = (absv1 > absv2) ? absv1 : absv2;
  int cmp;
  double absdelta;
  if (v1 > v2) {
    absdelta = v1 - v2;
    cmp = 1;
  } else {
    absdelta = v2 - v1;
    cmp = -1;
  }

  return (std::isnan(absdelta) ||
          (std::isfinite(absdelta) && (absdelta <= absmax * epsilon)))
             ? 0
             : cmp;
}

#if 0
static inline
int fuzzy_cmp_safe(double v1, double v2, double epsilon)
{
  if (std::isfinite(v1))
    {
      if (std::isfinite(v2))
        return fuzzy_cmp(v1, v2, epsilon);
      else
        return (std::isinf(v2) && std::signbit(v2)) ? 1 : -1; // v2 == nan is same as v2 == +inf
    }
  else
    {
      if (std::isinf(v1) && std::signbit(v1))
        return (std::isinf(v2) && std::signbit(v2)) ? 0 : -1;
      else
        return (std::isfinite(v2) || (std::isinf(v2) && std::signbit(v2))) ? 1 : 0;
    }
}
#endif // 0

} // namespace asyncde

#endif
