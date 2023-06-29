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

#include "asyncde/Status.h"

#include <limits>

void asyncde::Status::ResetAfterRecovery() {
  statusbits = 0;

  xepsilon = std::numeric_limits<double>::max();
  yepsilon = std::numeric_limits<double>::max();
  ydelta = std::numeric_limits<double>::max();

  Reset_nxdiff();
  Reset_nProjFeas();
  Reset_nFailedTrialPoint();

  error = 0;
}

void asyncde::Status::ResetAfterRestart() {
  Reset_nRecover();

  ResetAfterRecovery();
}

void asyncde::Status::ResetToFillInTrialPoint() {
  Reset_nProjFeas();
  Reset_nFailedTrialPoint();
}

unsigned int asyncde::Status::FindStatus(const StopCriteria &_criteria) const {
  unsigned int mask = 0;

  if ((_criteria.statusbits & CRITERION_vtr) != 0)
    if (vtr < _criteria.vtr)
      mask |= CRITERION_vtr;

  if ((_criteria.statusbits & CRITERION_nFE) != 0)
    if (nFE >= _criteria.nFE)
      mask |= CRITERION_nFE;

  if ((_criteria.statusbits & CRITERION_nprojfeas) != 0)
    if (nprojfeas >= _criteria.nprojfeas)
      mask |= CRITERION_nprojfeas;

  if ((_criteria.statusbits & CRITERION_nxdiff) != 0)
    if (nxdiff >= _criteria.nxdiff)
      mask |= CRITERION_nxdiff;

  if ((_criteria.statusbits & CRITERION_popsizemax) != 0)
    if (popsizemax > _criteria.popsizemax)
      mask |= CRITERION_popsizemax;

  if ((_criteria.statusbits & CRITERION_xepsilon) != 0)
    if (xepsilon < _criteria.xepsilon)
      mask |= CRITERION_xepsilon;

  if ((_criteria.statusbits & CRITERION_yepsilon) != 0)
    if (yepsilon < _criteria.yepsilon)
      mask |= CRITERION_yepsilon;

  if ((_criteria.statusbits & CRITERION_ydelta) != 0)
    if (ydelta < _criteria.ydelta)
      mask |= CRITERION_ydelta;

  if ((_criteria.statusbits & CRITERION_failedtrialpoint) != 0)
    if (nfailedtrialpoint > _criteria.nfailedtrialpoint)
      mask |= CRITERION_failedtrialpoint;

  if ((_criteria.statusbits & CRITERION_error) != 0)
    if (error)
      mask |= CRITERION_error;

  return mask;
}

int asyncde::Status::Print(FILE *stream, unsigned int mask) const {
  if (!stream)
    return -1;

  int retvalue = fprintf(stream, "asyncde::Status: ");
  retvalue |= StopCriteria::Print(stream, mask);
  if (retvalue < 0)
    return retvalue;

  return retvalue;
}
