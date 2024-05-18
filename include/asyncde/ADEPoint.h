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

#ifndef ASYNCDE_ADEPOINT_H
#define ASYNCDE_ADEPOINT_H

#include <asyncde/ADEPointInfo.h>
#include <asyncde/Point.h>
#include <asyncde/PointData.h>

namespace asyncde {

class ADEPoint : public Point {
protected:
  /// contains ids and information specific the ADE algorithm
  ADEPointInfo *ADEinfo; // owner

private:
  ADEPoint(ADEPoint &_point)
      : Point(*_point.DataMutable(), _point.InfoMutable()) {}

  ADEPoint(PointData &_data, PointInfo *_info) : Point(_data, _info) {}

public:
  /// Constructor: example new ADEPoint(*new PointData(..), new ADEPointInfo());
  ADEPoint(PointData &_data, ADEPointInfo *_info)
      : Point(_data, _info), ADEinfo(_info) {}

  const ADEPointInfo *ADEInfo() const { return ADEinfo; }

  ADEPointInfo *ADEInfoMutable() { return ADEinfo; }

  /// Allocate a copy of this
  virtual Point *Clone() const override {
    return new ADEPoint(*this->data->Clone(),
                        (ADEPointInfo *)this->ADEinfo->Clone());
  }

  virtual int Print(FILE *stream = stdout) const override;
};

} // namespace asyncde

#endif
