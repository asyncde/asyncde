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

#ifndef ASYNCDE_ASYNCITERATOR_H
#define ASYNCDE_ASYNCITERATOR_H

#include <stdio.h>
#include <vector>

#include "asyncde/Status.h"

namespace asyncde {

// forward declarations:
class IteratorConfig;
class Point;
class Problem;
class Status;

class AsyncIterator {
protected:
  /// problem definition: variables and constraints
  const Problem *problem; // not owner

  /// actual settings
  IteratorConfig *cfg; // owner

  /// current status of the iterator
  Status *status;

  /// status confronted to cfg->criteriastop
  unsigned int statusstopbits;
  /// status confronted to cfg->criteriarestart
  unsigned int statusrestartbits;
  /// increased by one if restart from initial population size has been
  /// performed
  int statusrestartcounter;
  /// status confronted to cfg->criteriarecover
  unsigned int statusrecoverbits;

  /// best point
  Point *bestintpoint;

protected:
  /// temporary point (internal variables) used by AddExtPoint() and
  /// FillInTrialExtPoint()
  Point *tmpintpoint;

  /// temporary point (internal variables) used by AddIntPoint() to
  /// keep the replaced point
  Point *tmpintoldpoint;

  /// temporary point (external variables) used by Minimize()
  Point *tmpextpoint;

public:
  AsyncIterator(const Problem &_problem, const IteratorConfig *_cfg);

  virtual ~AsyncIterator();

  /// Iterator allocates a point specific to this iterator (external variables)
  virtual Point *NewExtPoint() const;

  /// restart algorithm with a new population
  virtual int Restart(unsigned int _nparents) = 0;

  /// recovery algorithm for the population
  virtual int Recovery() = 0;

  /// optimization loop (single thread)
  int Minimize();

  /// optimization loop (multithread, objective function has to be thread-safe)
  int MinimizeMT(
      size_t ntrhreads = 0); /* through std::thread according to C++17 */
  //  int MinimizeOpenMP();

  /// Return a read-only reference to the best point
  const Point *BestIntPoint() const { return bestintpoint; }

  /// Return number of evaluated points
  unsigned long int NFE() const { return status->nFE; }

  /// Add point to Iterator; X-coordinates according to external variables
  int AddExtPoint(const Point &_point);

  /// Iterator fills in a trial point; X-coordinates according to external
  /// variables
  int FillInTrialExtPoint(Point &_point, int maxnrestarts = 1);

  /// return true if finished
  virtual int UpdateStatus() = 0;
  const Status *GetStatus() const { return status; }

  unsigned int StatusStopBits() const { return statusstopbits; }

  unsigned int StatusRestartBits() const { return statusrestartbits; }

  virtual int Print(FILE * /*stream = stdout*/) const { return 0; };

protected:
  /// Iterator allocates a point specific to this iterator (internal variables)
  virtual Point *NewIntPoint() const;

  /// Add point to Iterator; X-coordinates according to internal variables
  virtual int AddIntPoint(Point &_point) = 0;

  /// Iterator fills in a trial point; X-coordinates according to internal
  /// variables
  virtual int FillInTrialIntPoint(Point &_point, int maxnrestarts = 1) = 0;

  // project point xint to border of the feasible region along the line
  // (xintfeas, xint)
  int ProjectPoint2Border(std::vector<double> &xint,
                          const std::vector<double> &xintfeas,
                          std::vector<double> &xexttmparray) const;

  /// stochastically project point xint into feasible region along the line
  /// (xintfeas, xint), xintfeas is a known feasible point
  int ProjectPoint2Feasibility(std::vector<double> &xint,
                               const std::vector<double> &xintfeas,
                               std::vector<double> &xexttmparray) const;
};

} // namespace asyncde

#endif
