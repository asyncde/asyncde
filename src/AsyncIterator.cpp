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

#include "asyncde/AsyncIterator.h"

#include "asyncde/IteratorConfig.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Problem.h"
#include "asyncde/Rnd.h"
#include "asyncde/Status.h"
#include "asyncde/mathextra.h"

asyncde::AsyncIterator::AsyncIterator(const Problem &_problem,
                                      const IteratorConfig *)
    : problem(&_problem), cfg(nullptr), status(nullptr), statusstopbits(0),
      statusrestartbits(0), statusrestartcounter(0), statusrecoverbits(0),
      bestintpoint(nullptr), tmpintpoint(nullptr), tmpextpoint(nullptr) {
  status = new Status();
}

asyncde::AsyncIterator::~AsyncIterator() {
  delete cfg;
  delete status;
  delete bestintpoint;
  delete tmpintpoint;
  delete tmpextpoint;
}

asyncde::Point *asyncde::AsyncIterator::NewExtPoint() const {
  Point *newpoint = new Point(
      *new PointData(*problem->ExtVariables(), *problem->Y()), new PointInfo());
  newpoint->Reset();
  return newpoint;
}

asyncde::Point *asyncde::AsyncIterator::NewIntPoint() const {
  Point *newpoint = new Point(
      *new PointData(*problem->IntVariables(), *problem->Y()), new PointInfo());
  newpoint->Reset();
  return newpoint;
}

int asyncde::AsyncIterator::Minimize() {
  const Functor1D *fitnessfunctor = problem->ObjectiveFunctor1D();
  if (!fitnessfunctor) {
    statusstopbits = CRITERION_error;
    return -1;
  }
  statusstopbits = 0;

  int retvalue = 0;
  std::vector<double> y(problem->Y()->size());
  if (!tmpextpoint)
    tmpextpoint = NewExtPoint();
  for (;;) {
    FillInTrialExtPoint(*tmpextpoint);

    if ((retvalue |= StatusStopBits()))
      break;

    (*fitnessfunctor)(*tmpextpoint->Data()->X(), y);
    tmpextpoint->SetY(&y);
    AddExtPoint(*tmpextpoint);

    if ((retvalue |= StatusStopBits()))
      break;
  }

  if (cfg->verbose) {
    const Point *bestpoint = BestIntPoint();
    if (bestpoint) {
      printf("nevals = %li  bestvalue = %.15e\n", NFE(),
             (*bestpoint->Data()->Y())[0]);
      problem->ConvertInt2Ext(*bestpoint->Data()->X(),
                              *tmpextpoint->DataMutable()->Xmutable());
      for (std::vector<double>::const_iterator it =
               tmpextpoint->Data()->X()->begin();
           it != tmpextpoint->Data()->X()->end(); it++) {
        if (it != tmpextpoint->Data()->X()->begin())
          printf("  ");
        printf("%.15e", *it);
      }
      printf("\n");
    }
  }

  return retvalue;
}

int asyncde::AsyncIterator::MinimizeOpenMP() {
  const Functor1D *fitnessfunctor = problem->ObjectiveFunctor1D();
  if (!fitnessfunctor) {
    statusstopbits = CRITERION_error;
    return -1;
  }
  statusstopbits = 0;

  int retvalue = 0;
  Point *threadextpoint = 0;
  std::vector<double> *Ythread = 0;

#pragma omp parallel private(Ythread, threadextpoint) shared(retvalue)
  {
#pragma omp critical
    {
      threadextpoint = NewExtPoint();
      Ythread = new std::vector<double>(problem->Y()->size());
      (*Ythread)[0] = std::numeric_limits<double>::max();
    }

    while (retvalue == 0) {
#pragma omp critical(de_internal_op)
      { FillInTrialExtPoint(*threadextpoint); }

      if (threadextpoint->Info()->Id() != POINTSTATUS_UNDEFINED) {
        (*fitnessfunctor)(*threadextpoint->Data()->X(), *Ythread);
        threadextpoint->SetY(Ythread);
      }

      if (retvalue == 0 &&
          threadextpoint->Info()->Id() != POINTSTATUS_UNDEFINED)
#pragma omp critical(de_internal_op)
      {
        AddExtPoint(*threadextpoint);
        retvalue |= StatusStopBits();
      }
    }

#pragma omp critical
    {
      delete threadextpoint;
      delete Ythread;
    }
  }

  if (cfg->verbose) {
    const Point *bestpoint = BestIntPoint();
    if (bestpoint) {
      printf("nevals = %li  bestvalue = %.15e\n", NFE(),
             (*bestpoint->Data()->Y())[0]);
      problem->ConvertInt2Ext(*bestpoint->Data()->X(),
                              *tmpextpoint->DataMutable()->Xmutable());
      for (std::vector<double>::const_iterator it =
               tmpextpoint->Data()->X()->begin();
           it != tmpextpoint->Data()->X()->end(); it++) {
        if (it != tmpextpoint->Data()->X()->begin())
          printf("  ");
        printf("%.15e", *it);
      }
      printf("\n");
    }
  }

  return retvalue;
}

int asyncde::AsyncIterator::AddExtPoint(const Point &_point) {
  if (!tmpintpoint)
    tmpintpoint = NewIntPoint();

  tmpintpoint->InfoMutable()->Set(*_point.Info());
  tmpintpoint->DataMutable()->SetVariables(*problem->IntVariables());
  problem->ConvertExt2Int(*_point.Data()->X(),
                          *tmpintpoint->DataMutable()->Xmutable());
  tmpintpoint->SetY(_point.Data()->Y());

  status->ResetToAddPoint();

  return AddIntPoint(*tmpintpoint);
}

int asyncde::AsyncIterator::FillInTrialExtPoint(Point &_point,
                                                int maxnrestarts) {
  if (!tmpintpoint)
    tmpintpoint = NewIntPoint();

  status->ResetToFillInTrialPoint();
  int retvalue = FillInTrialIntPoint(*tmpintpoint, maxnrestarts);

  _point.InfoMutable()->Set(*tmpintpoint->Info());
  _point.DataMutable()->SetVariables(*problem->ExtVariables());
  problem->ConvertInt2Ext(*tmpintpoint->Data()->X(),
                          *_point.DataMutable()->Xmutable());
  _point.SetY(tmpintpoint->Data()->Y());

  UpdateStatus();

  return retvalue;
}

int asyncde::AsyncIterator::ProjectPoint2Border(
    std::vector<double> &xint, const std::vector<double> &xintfeas,
    std::vector<double> &xexttmparray) const {
  /*
   *  if (1 == IsXIntFeasible(xint, xexttmparray))
   * return 0;
   *
   *  if (!IsXIntFeasible(xintfeas, xexttmparray))
   *    return -1;
   */
  unsigned int ivar, nvars = problem->NFreeVariables();
  int same;
  double r = -1.0, rfeas = 0.0, rinfeas = 1.0;
  double xiold;

  xexttmparray = xint;

  do {
    r = 0.5 * (rfeas + rinfeas);

    same = 1;
    for (ivar = 0; ivar < nvars; ivar++) {
      xiold = xint[ivar];
      xint[ivar] = xintfeas[ivar] + r * (xexttmparray[ivar] - xintfeas[ivar]);
      if (same && (0 != fuzzy_cmp(xint[ivar], xiold, problem->Tolerance())))
        same = 0;
    }

    if (1 == problem->IsXIntFeasible(xint, xexttmparray))
      rfeas = r;
    else
      rinfeas = r;

  } while (!same);

  return 0;
}

int asyncde::AsyncIterator::ProjectPoint2Feasibility(
    std::vector<double> &xint, const std::vector<double> &xintfeas,
    std::vector<double> &xexttmparray) const {
  /*
   *  if (!IsXIntFeasible(xintfeas, xexttmparray))
   *    return -1;
   */
  unsigned int ivar, nvars = problem->NFreeVariables();
  double r;
  int same;

  do {
    r = cfg->rnd->next(0.5, 1.0);

    for (ivar = 0; ivar < nvars; ivar++)
      xint[ivar] = xintfeas[ivar] + r * (xint[ivar] - xintfeas[ivar]);

    same = 1;
    for (ivar = 0; ivar < nvars; ivar++)
      if (0 != fuzzy_cmp(xint[ivar], xintfeas[ivar], problem->Tolerance())) {
        same = 0;
        break;
      }

    if (!same && 1 == problem->IsXIntFeasible(xint, xexttmparray))
      return 0;
  } while (!same);

  return -1;
}
