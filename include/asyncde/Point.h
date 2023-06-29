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

#ifndef ASYNCDE_POINT_H
#define ASYNCDE_POINT_H

#include <asyncde/PointData.h>
#include <asyncde/PointInfo.h>
#include <asyncde/Variables.h>

namespace asyncde {

class Point {
protected:
  /// contains arguments and function values
  PointData *data; // owner
  /// contains ids and information specific to algorithms
  PointInfo *info; // owner

private:
  Point(const Point &) {}
  Point(Point &) {}

public:
  /// Constructor: example new Point(*new PointData(..), new PointInfo());
  Point(PointData &_data, PointInfo *_info) : data(&_data), info(_info) {}

  virtual ~Point() {
    delete data;
    delete info;
  }

  const PointData *Data() const { return data; }

  PointData *DataMutable() { return data; }

  const PointInfo *Info() const { return info; }

  PointInfo *InfoMutable() { return info; }

  /// Allocate a copy of this
  virtual Point *Clone() const {
    return new Point(*this->data->Clone(), this->info->Clone());
  }

  /// Reset point (set its fields to default values)
  virtual void Reset() {
    data->Reset();
    info->Reset();
  }

  /// Set this point by contents of _point
  int Set(const Point &_point) {
    int retvalue = data->Set(*_point.Data());
    retvalue |= info->Set(*_point.Info());
    return retvalue;
  }

  /// Set arguments and update status
  int SetX(const std::vector<double> *_x) {
    int retvalue = data->SetX(_x);
    if (0 == retvalue) {
      if (info->Status() < POINTSTATUS_INITIALIZED)
        info->SetStatus(POINTSTATUS_INITIALIZED);
    } else
      info->SetStatus(POINTSTATUS_UNDEFINED);

    return retvalue;
  }

  /// Set function values and update status
  int SetY(const std::vector<double> *_y) {
    int retvalue = data->SetY(_y);
    if (0 == retvalue)
      info->SetStatus(POINTSTATUS_EVALUATED);
    else if (info->Status() >= POINTSTATUS_EVALUATED)
      info->SetStatus(POINTSTATUS_UNDEFINED);

    return retvalue;
  }

  virtual int Print(FILE *stream = stdout) const;
};

} // namespace asyncde

#endif
