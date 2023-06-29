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

#ifndef ASYNCDE_PROBLEMMAPPING_H
#define ASYNCDE_PROBLEMMAPPING_H

#include "asyncde/Problem.h"

#include <stdio.h>

namespace asyncde {

/*
 Class with automatic exclusion of fixed variables; Internal variables contain
 only free variables
*/
class ProblemMapping : public Problem {
protected:
  /// used to convert from internal to external variables;
  std::vector<int> int2extmap; //[nfreevars]
  /// used to convert from external to internal variables
  std::vector<int> ext2intmap; // [nvars]

public:
  ProblemMapping() : Problem() {}

  virtual ~ProblemMapping() {}

  virtual int UpdateIntVars() override;

  virtual int ConvertInt2Ext(const std::vector<double> &xint,
                             std::vector<double> &xext) const override;

  virtual int ConvertExt2Int(const std::vector<double> &xext,
                             std::vector<double> &xint) const override;

  virtual int Print(FILE *stream = stdout) const override;
};

} // namespace asyncde

#endif
