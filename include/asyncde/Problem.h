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

#ifndef ASYNCDE_PROBLEM_H
#define ASYNCDE_PROBLEM_H

#include <limits>
#include <stdio.h>
#include <vector>

#include "asyncde/Function.h"
#include "asyncde/Variables.h"

namespace asyncde {

/*
 Class to hold
   - vector of variables (external representation)
   - rules to convert from/to external to/from internal variables
   - constraints
 */
class Problem {
protected:
  /// variables (X) - external representation
  const Variables *extvars;

  /// variables (X) - internal representation (all should be free variables!)
  Variables intvars;

  /// objective functor
  const Functor1D *objfunctor;

  /// expected vector of objective function values
  std::vector<double> y;

  /// variable to bypass feasibility test for unconstrained problems
  int activeconstraints;

  /// Statistical scale used to calculate the error; typically 1 for
  /// chi-squared, 0.5 for likelihood maximization; negative - do not track
  /// errors
  double nsigma;

  /// precision used in floating point comparisons
  double tolerance;

public:
  Problem()
      : extvars(0), intvars(0), objfunctor(0), activeconstraints(0), nsigma(-1),
        tolerance(8 * std::numeric_limits<double>::epsilon()) {
    y.resize(1);
  }

  virtual ~Problem();

  unsigned int NVariables() const {
    return (extvars) ? extvars->NVariables() : 0;
  }

  unsigned int NFreeVariables() const { return intvars.NVariables(); }

  const Variables *ExtVariables() const { return extvars; }

  int SetExternalVariables(const Variables *_vars);

  const Variables *IntVariables() const { return &intvars; }

  virtual int UpdateIntVars();

  virtual int ConvertInt2Ext(const std::vector<double> &xint,
                             std::vector<double> &xext) const = 0;

  virtual int ConvertExt2Int(const std::vector<double> &xext,
                             std::vector<double> &xint) const = 0;

  // TODO: virtual double XPenalty(const std::vector<double> & xext) const;

  // Test whether point xext (external coordinates) is feasible
  virtual int IsXFeasible(const std::vector<double> &xext) const;

  // Interface to test whether point xint (internal coordinates) is feasible:
  // calls XFeasible()
  int IsXIntFeasible(const std::vector<double> &xint,
                     std::vector<double> &xexttmparray) const {
    if (!activeconstraints)
      return 1;
    ConvertInt2Ext(xint, xexttmparray);
    return IsXFeasible(xexttmparray);
  }

  const Functor1D *ObjectiveFunctor1D() const { return objfunctor; }

  void SetObjectiveFunctor1D(const Functor1D *_objfunctor) {
    objfunctor = _objfunctor;
  }

  const std::vector<double> *Y() const { return &y; }

  /*
   // TODO: dedicated class to split dimension of a function from y.size()
      unsigned int NDimensionY() const
        {return y.size();}
  */
  void SetYsize(unsigned int _newsize) {
    if (_newsize > 0)
      y.resize(_newsize);
  }

  double ErrorDef() const { return nsigma; }
  void SetErrorDef(double _nsigma) { nsigma = _nsigma; }

  double Tolerance() const { return tolerance; }
  void SetTolerance(double _tolerance) { tolerance = _tolerance; }

  virtual int Print(FILE *stream = stdout) const;
};

} // namespace asyncde

#endif
