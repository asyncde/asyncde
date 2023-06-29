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

#ifndef ASYNCDE_ITERATORCONFIG_H
#define ASYNCDE_ITERATORCONFIG_H

#include "asyncde/StopCriteria.h"

#include <stdio.h>

namespace asyncde {

// forward declarations:
class AsyncIterator;
class Problem;
class Rnd;

class IteratorConfig {
public:
  StopCriteria criteriastop;
  StopCriteria criteriarestart;
  StopCriteria criteriarecover;
  Rnd *rnd;

  /// level of verbosity: 0 = quiet, 3 - maximal verbosity
  int verbose;

protected:
  void Init() {
    criteriastop.SetDefaultStopCriteria();
    criteriarestart.SetDefaultRestartCriteria();
    criteriarecover.statusbits = 0;
    verbose = 0;
  }

public:
  IteratorConfig(Rnd *_rnd, unsigned int /*_minpopsize*/ = 0) : rnd(_rnd) {
    Init();
  }

  virtual ~IteratorConfig() {}

  virtual void Reset(unsigned int /* _minpopsize */) { Init(); }

  virtual AsyncIterator *NewIterator(const Problem &_problem) const = 0;

  virtual int Print(FILE *stream = stdout) const;
};

} // namespace asyncde

#endif
