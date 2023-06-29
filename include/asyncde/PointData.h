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

#ifndef ASYNCDE_POINTDATA_H
#define ASYNCDE_POINTDATA_H

#include <limits>
#include <stdio.h>
#include <vector>

#include <asyncde/Variables.h>

namespace asyncde {

// forward declarations

enum PointDataErrCodes {
  POINTDATAERR_XSIZE = -1,
  POINTDATAERR_YSIZE = -2,
  POINTDATAERR_XUNDEFINED = -4,
  POINTDATAERR_YUNDEFINED = -8,
  POINTDATAERR_VARSID = -16
};

class PointData {
protected:
  /// ID of a parent set of Variables
  long int vars_id;

  /// Link to the parent set of Variables
  const Variables *vars;

  /// Point coordinates
  std::vector<double> x;

  // Function values
  std::vector<double> y;

public:
  PointData(const Variables &_vars, const std::vector<double> &_y)
      : vars(&_vars) {
    vars_id = vars->Id();
    x = *vars->X();
    y = _y;
  }

  virtual ~PointData() {}

  /// Set this point by contents of _pointdata
  virtual int Set(const PointData &_pointdata) {
    int retvalue = 0;

    if (_pointdata.VariablesId() != VariablesId())
      return POINTDATAERR_VARSID;

    retvalue |= SetX(_pointdata.X());
    y = *_pointdata.Y();

    return retvalue;
  }

  /// Allocate a copy of this
  virtual PointData *Clone() const {
    PointData *newpointdata = new PointData(*vars, y);
    newpointdata->SetX(&x);

    return newpointdata;
  }

  /// Reset point (set its fields to default values)
  virtual void Reset() {
    SetX(0);
    SetY(0);
  }

  long int VariablesId() const { return vars_id; }
  const Variables *Link2Variables() const { return vars; }
  unsigned int NVariables() const { return x.size(); }
  void SetVariables(const Variables &_vars) {
    vars = &_vars;
    vars_id = _vars.Id();
  }

  const std::vector<double> *X() const { return &x; }
  std::vector<double> *Xmutable() { return &x; }

  int SetX(const std::vector<double> *_x) {
    if (_x && x.size() == _x->size()) {
      x = *_x;
      return 0;
    } else {
      std::fill(x.begin(), x.end(), -std::numeric_limits<double>::max());
      return POINTDATAERR_XUNDEFINED;
    }
  }

  unsigned int NDimensionY() const { return 1; }

  const std::vector<double> *Y() const { return &y; }

  std::vector<double> *Ymutable() { return &y; }

  int SetY(const std::vector<double> *_y) {
    if (_y && y.size() == _y->size()) {
      y = *_y;
      return 0;
    } else {
      std::fill(y.begin(), y.end(), std::numeric_limits<double>::max());
      return POINTDATAERR_YUNDEFINED;
    }
  }

  virtual int Print(FILE *stream = stdout) const;
};

} // namespace asyncde

#endif
