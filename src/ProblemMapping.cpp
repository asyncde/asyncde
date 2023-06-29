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

#include "asyncde/ProblemMapping.h"

#include "asyncde/Variables.h"

int asyncde::ProblemMapping::UpdateIntVars() {
  unsigned int ivar;
  int retvalue = Problem::UpdateIntVars();

  if (0 != retvalue)
    return retvalue;

  const unsigned int nvars = extvars->NVariables();
  const unsigned int nfreevars = extvars->FindNFreeVariables();

  int2extmap.resize(nfreevars);
  ext2intmap.resize(nvars);

  unsigned int intcounter = 0;
  for (ivar = 0; ivar < nvars; ivar++)
    if (1 == extvars->IsFreeVariable(ivar)) {
      ext2intmap[ivar] = intcounter;
      int2extmap[intcounter] = ivar;
      intcounter++;
    } else
      ext2intmap[ivar] = -1;

  if (nfreevars != intcounter)
    return -4;

  // init intvars
  for (ivar = 0; ivar < nfreevars; ivar++)
    retvalue |= intvars.CopyVariable(ivar, *extvars, int2extmap[ivar]);

  return retvalue;
}

int asyncde::ProblemMapping::ConvertInt2Ext(const std::vector<double> &xint,
                                            std::vector<double> &xext) const {
  const unsigned int nfreevars = intvars.NVariables();

  if (extvars)
    xext = *extvars->X(); // in fact should copy only fixed parameters
  for (unsigned int ivar = 0; ivar < nfreevars; ivar++)
    xext[int2extmap[ivar]] = xint[ivar];

  return 0;
}

int asyncde::ProblemMapping::ConvertExt2Int(const std::vector<double> &xext,
                                            std::vector<double> &xint) const {
  if (!extvars)
    return -1;

  const unsigned int nvars = extvars->NVariables();
  int index;
  for (unsigned int ivar = 0; ivar < nvars; ivar++)
    if (0 <= (index = ext2intmap[ivar]))
      xint[index] = xext[ivar];

  return 0;
}

int asyncde::ProblemMapping::Print(FILE *stream) const {
  if (!stream)
    return -1;

  int retvalue = fprintf(stream, "asyncde::ProblemMapping: ");
  retvalue |= Problem::Print(stream);

  return retvalue;
}
