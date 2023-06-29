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

#ifndef ASYNCDE_STOPCRITERIA_H
#define ASYNCDE_STOPCRITERIA_H

#include <limits>
#include <stdio.h>

namespace asyncde {

enum {
  CRITERION_vtr = 0x0001,
  CRITERION_nFE = 0x0002,
  CRITERION_popsizemax = 0x0004,
  CRITERION_nprojfeas = 0x0008, /* failed to generate a trial point */

  CRITERION_xepsilon = 0x0010, /* spread of at least one coordinate is below
                                |xepsilonmin * xmax_i| */
  CRITERION_nxdiff =
      0x0020, /* at least one coordinate < xepsilon in the diff vector */
  //  CRITERION_xdelta               = 0x0040, /* all spreads in coordinates are
  //  below xdelta_i (can differ for different coordinates) */

  CRITERION_yepsilon = 0x0100, /* spread in values is below |yepsilon * ymax| */
  CRITERION_ydelta = 0x0200,   /* spread in values is below ydelta */

  //  CRITERION_nFE_after_restart    = 0x01000,
  //  CRITERION_nFE_after_best       = 0x02000,
  //  CRITERION_nFE_after_progress   = 0x04000,

  CRITERION_failedtrialpoint = 0x1000, /* failed to generate a trial point */
  CRITERION_error = 0x10000,
};

class StopCriteria {
  // protected:
public:
  unsigned int statusbits; // keeps active criteria according to CRITERIA
                           // (Status uses it to keep fulfilled criteria)

  double vtr;

  long unsigned int nFE;

  double xepsilon;
  double yepsilon;
  double ydelta;

  int popsizemax;

  unsigned int nxdiff;
  unsigned int nprojfeas;
  unsigned int nrecover;
  unsigned int nfailedtrialpoint;

  int error;

protected:
  void Init() {
    statusbits = 0;
    error = 0;
    nxdiff = 0;
    nprojfeas = 0;
    nrecover = 0;
    nfailedtrialpoint = 0;
    popsizemax = -1;
  }

public:
  StopCriteria() { Init(); }

  virtual ~StopCriteria() {}

  virtual void Reset() { Init(); }

  void SetDefaultRestartCriteria() {
    SetCriterion_xepsilon(8 * std::numeric_limits<double>::epsilon());
    SetCriterion_yepsilon(8 * std::numeric_limits<double>::epsilon());
    SetCriterion_nProjFeas(100);
  }

  void SetDefaultStopCriteria() { statusbits = CRITERION_error; }

  void SetCriterion_vtr(double _vtr) {
    statusbits |= CRITERION_vtr;
    vtr = _vtr;
  }

  void SetCriterion_nFE(long unsigned int _nFE) {
    statusbits |= CRITERION_nFE;
    nFE = _nFE;
  }

  void SetCriterion_xepsilon(double _xepsilon) {
    statusbits |= CRITERION_xepsilon;
    xepsilon = _xepsilon;
  }

  void SetCriterion_nxdiff(unsigned int _nxdiff) {
    statusbits |= CRITERION_nxdiff;
    nxdiff = _nxdiff;
  }

  void SetCriterion_yepsilon(double _yepsilon) {
    statusbits |= CRITERION_yepsilon;
    yepsilon = _yepsilon;
  }

  void SetCriterion_ydelta(double _ydelta) {
    statusbits |= CRITERION_ydelta;
    ydelta = _ydelta;
  }

  void SetCriterion_nProjFeas(unsigned int _nprojfeas) {
    statusbits |= CRITERION_nprojfeas;
    nprojfeas = _nprojfeas;
  }

  void SetCriterion_PopSizeMax(unsigned int _popsizemax) {
    statusbits |= CRITERION_popsizemax;
    popsizemax = _popsizemax;
  }

  void SetCriterion_FailedTrialPoint(unsigned int _nfailedtrialpoint = 0) {
    statusbits |= CRITERION_failedtrialpoint;
    nfailedtrialpoint = _nfailedtrialpoint;
  }

  virtual int Print(FILE *stream = stdout,
                    unsigned int mask = 0xffffffff) const;
};

} // namespace asyncde

#endif
