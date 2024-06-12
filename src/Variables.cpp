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

#include <stdio.h>

#include "asyncde/Variables.h"
#include "asyncde/mathextra.h"

long int asyncde::Variables::id_counter = -1;

asyncde::Variables::Variables(unsigned int nvariables_to_reserve) {
  Resize(nvariables_to_reserve);
  SwitchToNewId();
}

asyncde::Variables::Variables(const Variables &_vars) {
  id = _vars.Id();
  name = *_vars.Names();
  x = *_vars.X();
  freevar = *_vars.Mask();
  minlimit = *_vars.LowerLimits();
  maxlimit = *_vars.UpperLimits();
  mininitial = *_vars.InitialLowerLimits();
  maxinitial = *_vars.InitialUpperLimits();
}

void asyncde::Variables::Reserve(unsigned int nvars_newmaxsize) {
  if (nvars_newmaxsize <= x.max_size())
    return;

  name.reserve(nvars_newmaxsize);
  x.reserve(nvars_newmaxsize);
  freevar.reserve(nvars_newmaxsize);
  minlimit.reserve(nvars_newmaxsize);
  maxlimit.reserve(nvars_newmaxsize);
  mininitial.reserve(nvars_newmaxsize);
  maxinitial.reserve(nvars_newmaxsize);

  return;
}

void asyncde::Variables::Resize(unsigned int nvars_newsize) {
  unsigned int ivar = x.size();

  name.resize(nvars_newsize);
  x.resize(nvars_newsize);
  freevar.resize(nvars_newsize);
  minlimit.resize(nvars_newsize);
  maxlimit.resize(nvars_newsize);
  mininitial.resize(nvars_newsize);
  maxinitial.resize(nvars_newsize);

  for (; ivar < nvars_newsize; ivar++)
    DefaultVariable(ivar);

  return;
}

int asyncde::Variables::AddVariable(const std::string _name, double value) {
  unsigned int index = NVariables();
  Resize(index + 1);

  name[index] = _name;
  x[index] = value;
  SwitchToNewId();

  return index;
}

int asyncde::Variables::CopyVariable(unsigned int intindex,
                                     const Variables &extvars, int extindex) {
  if (intindex >= NVariables())
    Resize(intindex + 1);

  if (extindex < 0 || extindex >= (int)extvars.NVariables())
    return -2;

  name[intindex] = extvars.Name(extindex);
  x[intindex] = (*extvars.X())[extindex];
  freevar[intindex] = extvars.IsFreeVariable(extindex);
  minlimit[intindex] = (*extvars.LowerLimits())[extindex];
  maxlimit[intindex] = (*extvars.UpperLimits())[extindex];
  mininitial[intindex] = (*extvars.InitialLowerLimits())[extindex];
  maxinitial[intindex] = (*extvars.InitialUpperLimits())[extindex];
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::DefaultVariable(unsigned int ivar) {
  if (ivar >= NVariables())
    return -1;

  name[ivar] = "";
  x[ivar] = 0.0;
  freevar[ivar] = 0;
  minlimit[ivar] = DefaultLowerLimit();
  maxlimit[ivar] = DefaultUpperLimit();
  mininitial[ivar] = DefaultLowerLimit();
  maxinitial[ivar] = DefaultUpperLimit();
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::SetName(unsigned int ivar, const std::string _name) {
  if (ivar >= NVariables())
    return -1;

  name[ivar] = _name;
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::Find(const std::string _name) const {
  int index = -1;
  for (std::vector<std::string>::const_iterator itname = name.cbegin();
       itname != name.cend(); itname++)
    if (_name == *itname) {
      index = itname - name.cbegin();
      break;
    }

  return index;
}

int asyncde::Variables::SetCoordinate(unsigned int ivar, double value) {
  if (ivar >= NVariables())
    return -1;

  x[ivar] = value;
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::SetX(const std::vector<double> &xnew) {
  if (x.size() != xnew.size())
    return -1;

  x = xnew;
  SwitchToNewId();

  unsigned int nvars = NVariables();
  for (unsigned int ivar = 0; ivar < nvars; ivar++)
    if (x[ivar] < minlimit[ivar] || x[ivar] > maxlimit[ivar])
      return -2;

  return 0;
}

unsigned int asyncde::Variables::FindNFreeVariables() const {
  unsigned int nfree = 0;

  for (int vfree : freevar)
    if (1 == vfree)
      nfree++;

  return nfree;
}

int asyncde::Variables::SetLowerLimit(unsigned int ivar, double value) {
  if (ivar >= NVariables())
    return -1;

  minlimit[ivar] = value;
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::SetUpperLimit(unsigned int ivar, double value) {
  if (ivar >= NVariables())
    return -1;

  maxlimit[ivar] = value;
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::SetLimits(const std::vector<double> &lowerlimits,
                                  const std::vector<double> &upperlimits) {
  unsigned int nvars = NVariables();
  if (lowerlimits.size() != nvars || upperlimits.size() != nvars)
    return -1;

  minlimit = lowerlimits;
  maxlimit = upperlimits;
  SwitchToNewId();

  int retvalue = 0;
  // consistency check
  for (unsigned int ivar = 0; ivar < nvars; ivar++)
    if (maxlimit[ivar] > minlimit[ivar])
      freevar[ivar] = 1;
    else {
      freevar[ivar] = 0;
      if (maxlimit[ivar] < minlimit[ivar])
        retvalue = -2;
    }

  return retvalue;
}

int asyncde::Variables::SetInitialLimits(
    const std::vector<double> &lowerlimits,
    const std::vector<double> &upperlimits) {
  unsigned int nvars = NVariables();
  if (lowerlimits.size() != nvars || upperlimits.size() != nvars)
    return -1;

  mininitial = lowerlimits;
  maxinitial = upperlimits;
  for (unsigned int ivar = 0; ivar < nvars; ivar++)
    freevar[ivar] = (maxlimit[ivar] > minlimit[ivar]) ? 1 : 0;

  SwitchToNewId();

  return IsConsistent();
}

int asyncde::Variables::SetInitialRange(unsigned int ivar, double min,
                                        double max) {
  if (ivar >= NVariables())
    return -1;

  mininitial[ivar] = min;
  maxinitial[ivar] = max;
  SwitchToNewId();

  return 0;
}

int asyncde::Variables::IsVariableConsistent(unsigned int ivar) const {
  if (ivar >= NVariables())
    return -1;

  if (x[ivar] < mininitial[ivar] || x[ivar] > maxinitial[ivar])
    return -2;

  if (minlimit[ivar] > mininitial[ivar] || maxlimit[ivar] < maxinitial[ivar])
    return -4;

  if (mininitial[ivar] >= maxinitial[ivar])
    return -8;

  return 1;
}

int asyncde::Variables::IsConsistent(int *index_failed) const {
  int retvalue;
  unsigned int nvars = NVariables();

  if (name.size() != nvars || x.size() != nvars || freevar.size() != nvars ||
      minlimit.size() != nvars || maxlimit.size() != nvars ||
      mininitial.size() != nvars || maxinitial.size() != nvars) {
    if (index_failed)
      *index_failed = nvars;
    return -1;
  }

  for (unsigned int ivar = 0; ivar < nvars; ivar++)
    if ((retvalue = IsVariableConsistent(ivar)) < 0) {
      if (index_failed)
        *index_failed = ivar;
      return retvalue;
    }

  return 1;
}

int asyncde::Variables::IsEquivalent(const Variables &_vars,
                                     double epsilon) const {
  if (1 != IsConsistent() || 1 != _vars.IsConsistent())
    return 0;

  unsigned int nvars = NVariables();
  if (nvars != _vars.NVariables())
    return 0;

  for (unsigned int ivar = 0; ivar < nvars; ivar++) {
    if (name[ivar] != (*_vars.Names())[ivar])
      return 0;
    if (freevar[ivar] != (*_vars.Mask())[ivar])
      return 0;

    if (0 != fuzzy_cmp(minlimit[ivar], (*_vars.LowerLimits())[ivar], epsilon) ||
        0 != fuzzy_cmp(maxlimit[ivar], (*_vars.UpperLimits())[ivar], epsilon) ||
        0 != fuzzy_cmp(mininitial[ivar], (*_vars.InitialLowerLimits())[ivar],
                       epsilon) ||
        0 != fuzzy_cmp(maxinitial[ivar], (*_vars.InitialUpperLimits())[ivar],
                       epsilon))
      return 0;
  }

  return 1;
}

int asyncde::Variables::Print(FILE *stream) const {
  if (!stream)
    return -1;

  int retvalue, status;
  unsigned int nvars = NVariables();

  if (0 > (retvalue = fprintf(
               stream, "asyncde::Variables(id=%li) Number of parameters: %i",
               id, nvars)))
    return retvalue;

  int index_failed;
  if (1 == (status = IsConsistent(&index_failed))) {
    if (0 > (retvalue = fprintf(stream, " Status is Ok\n")))
      return retvalue;
  } else {
    if (0 > (retvalue = fprintf(
                 stream, " Set is inconsistent (%i), first failed index %i\n",
                 status, index_failed)))
      return retvalue;
  }

  for (unsigned int ivar = 0; ivar < nvars; ivar++) {
    if (0 > (retvalue = fprintf(stream, "%i: %s x=%.7e", ivar,
                                name[ivar].c_str(), x[ivar])))
      return retvalue;
    if (freevar[ivar]) {
      if (0 > (retvalue = fprintf(
                   stream, " in [%.5e, %.5e] initial range in [%.5e, %.5e]",
                   minlimit[ivar], maxlimit[ivar], mininitial[ivar],
                   maxinitial[ivar])))
        return retvalue;
    } else {
      if (0 > (retvalue = fprintf(stream, " fixed")))
        return retvalue;
    }
    if (1 == (status = IsVariableConsistent(ivar))) {
      if (0 > (retvalue = fprintf(stream, " Ok\n")))
        return retvalue;
    } else if (0 > (retvalue = fprintf(stream, " Inconsistent (%i)\n", status)))
      return retvalue;
  }

  return 0;
}
