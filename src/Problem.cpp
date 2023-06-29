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

#include <stdlib.h>

#include "asyncde/Problem.h"
#include "asyncde/Rnd.h"
#include "asyncde/Variables.h"
#include "asyncde/mathextra.h"

int asyncde::Problem::UpdateIntVars() {
  intvars.Reset();

  if (!extvars)
    return -1;

  activeconstraints = 0;
  const std::vector<double> *xboxlow = extvars->LowerLimits();
  const std::vector<double> *xboxup = extvars->UpperLimits();
  unsigned int nextvars = NVariables();
  for (unsigned int iextvar = 0; iextvar < nextvars; iextvar++)
    if (extvars->IsFreeVariable(iextvar) &&
        ((*xboxlow)[iextvar] > extvars->DefaultLowerLimit() ||
         (*xboxup)[iextvar] < extvars->DefaultUpperLimit())) {
      activeconstraints = 1;
      break;
    }

  unsigned int nfreevars = extvars->FindNFreeVariables();
  if (nfreevars < 1)
    return -2;

  intvars.Reserve(nfreevars);

  return 0;
}

asyncde::Problem::~Problem() { delete extvars; }

int asyncde::Problem::SetExternalVariables(const Variables *_vars) {
  if (extvars)
    delete extvars;

  extvars = (_vars) ? new Variables(*_vars) : 0;
  return UpdateIntVars();
}

int asyncde::Problem::IsXFeasible(const std::vector<double> &xext) const {
  if (!extvars)
    return -1;

  if (!activeconstraints)
    return 1;

  const std::vector<double> *xboxlow = extvars->LowerLimits();
  const std::vector<double> *xboxup = extvars->UpperLimits();
  unsigned int ivar, nvars = NVariables();

  for (ivar = 0; ivar < nvars; ivar++)
    if (xext[ivar] < (*xboxlow)[ivar])
      return 0;

  for (ivar = 0; ivar < nvars; ivar++)
    if (xext[ivar] > (*xboxup)[ivar])
      return 0;

  return 1;
}

int asyncde::Problem::Print(FILE *stream) const {
  int retvalue;

  if (!stream)
    return -1;

  retvalue = fprintf(stream, "asyncde::Problem:\n");

  fprintf(stream, "External variables:");
  if (extvars)
    retvalue |= extvars->Print(stream);
  else
    fprintf(stream, " NULL\n");

  fprintf(stream, "Active constraints: %i\n", activeconstraints);

  fprintf(stream, "Internal variables:");
  retvalue |= intvars.Print(stream);

  fprintf(stream, "Objective functor (%p)\n", objfunctor);

  fprintf(stream, "Tolerance: %.3e\n", tolerance);
  if (nsigma > 0.0)
    fprintf(stream, "ErrorDef: %.3e\n", nsigma);

  return retvalue;
}
