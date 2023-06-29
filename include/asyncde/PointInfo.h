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

#ifndef ASYNCDE_POINTINFO_H
#define ASYNCDE_POINTINFO_H

#include <stdio.h>

namespace asyncde {

enum PointStatus {
  POINTSTATUS_UNDEFINED = -1,
  POINTSTATUS_INITIALIZED = 0,
  POINTSTATUS_REQUESTED,
  POINTSTATUS_INFEASIBLE,
  POINTSTATUS_EVALUATED,
  POINTSTATUS_NSTATUSES
};
/*
 * Class to store point's information specific to a minimization algorithm
 */
class PointInfo {
public:
  /// Reserved identifier for allocated but not initialized points
  static const int NOT_INITIALIZED_POINT = -1;

protected:
  /// Stores point id used up to now (used to generate unique point id)
  static long int id_counter;

  /// Unique point identifier
  long int id;

  /// Point status
  PointStatus status;

public:
  PointInfo() : id(NOT_INITIALIZED_POINT), status(POINTSTATUS_UNDEFINED) {
    id_counter++;
    id = id_counter;
  }

  virtual ~PointInfo() {}

  /// Set this point by contents of _pointinfo
  virtual int Set(const PointInfo &_pointinfo) {
    id = _pointinfo.Id();
    SetStatus(_pointinfo.Status());

    return 0;
  }

  /// Allocate a copy of this
  virtual PointInfo *Clone() const {
    PointInfo *newpointinfo = new PointInfo();
    newpointinfo->SetId(Id());
    newpointinfo->SetStatus(Status());

    return newpointinfo;
  }

  /// Reset point (set its fields to default values)
  virtual void Reset() {
    id = NOT_INITIALIZED_POINT;
    SetStatus(POINTSTATUS_UNDEFINED);
  }

  PointStatus Status() const { return status; }
  void SetStatus(PointStatus _status) { status = _status; }

  /// Return the highest point ID already used
  long int IdCounter() const { return id_counter; }

  /// Reset global ID counter (can be used only after complete finish of
  /// iterations)
  static void ResetIdCounter() { id_counter = -1; }

  long int Id() const { return id; }

  /// Switch to a new unique ID
  virtual void SwitchToNewId() {
    id_counter++;
    id = id_counter;
  }

  void SetId(long int _id) { id = _id; }

  virtual int Print(FILE *stream = stdout) const;
};

} // namespace asyncde

#endif
