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

#ifndef ASYNCDE_STATUS_H
#define ASYNCDE_STATUS_H

#include "asyncde/StopCriteria.h"

#include <limits>
#include <math.h>

namespace asyncde {

/// expects that IEEE 754 inf is returned for overflow due to division by zero
static inline double find_abs_ratio(double v1, double v2) {
  double absv1 = fabs(v1);
  double absv2 = fabs(v2);
  double absmax = (absv1 > absv2) ? absv1 : absv2;

  return (v1 > v2) ? (v1 - v2) / absmax : (v2 - v1) / absmax;
}

class Status : public StopCriteria {
protected:
  void Init() {
    vtr = std::numeric_limits<double>::max();
    nFE = 0;
    naccepted = 0;

    ResetAfterRestart();
  }

public:
  Status() : StopCriteria() { Init(); }

  Status(const Status &_status) : StopCriteria(_status) {}

  virtual ~Status() {}

  virtual void Reset() {
    StopCriteria::Reset();
    Init();
  }

  void ResetAfterRecovery();

  void ResetAfterRestart();

  void ResetToAddPoint() {}

  void ResetToFillInTrialPoint();

  /// return statusbits mask by confronting this with _criteria
  unsigned int FindStatus(const StopCriteria &_criteria) const;

  void Incr_nFE() { nFE++; }

  /// change status if new best is found
  void NewBestFound(int population_complete) {
    ProgressFound(population_complete);
  }

  /// change status if population has been improved
  void ProgressFound(int /* population_complete */) { naccepted++; }

  void Reset_nxdiff() { nxdiff = 0; }

  void Incr_nxdiff() { nxdiff++; }

  void Reset_nProjFeas() { nprojfeas = 0; }

  void Incr_nProjFeas() { nprojfeas++; }

  void Reset_nFailedTrialPoint() { nfailedtrialpoint = 0; }

  void Incr_nFailedTrialPoint() { nfailedtrialpoint++; }

  void Reset_nRecover() { nrecover = 0; }

  void Incr_nRecover() { nrecover++; }

  void ResetXYEpsilons() {
    xepsilon = std::numeric_limits<double>::max();
    yepsilon = 0.0;
    ydelta = 0.0;
  }

  void AddMinMaxPair(double v1, double v2, int x) {
    double new_epsilon = find_abs_ratio(v1, v2);

    if (x) {
      if (new_epsilon < xepsilon)
        xepsilon = new_epsilon;
    } else {
      if (new_epsilon > yepsilon)
        yepsilon = new_epsilon;

      double new_ydelta = fabs(v1 - v2);
      if (new_ydelta > ydelta)
        ydelta = new_ydelta;
    }
  }

  virtual int Print(FILE *stream = stdout,
                    unsigned int mask = 0xffffffff) const;
};

} // namespace asyncde

#endif
