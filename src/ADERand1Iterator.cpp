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

#include <algorithm>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "asyncde/ADERand1Config.h"
#include "asyncde/ADERand1Iterator.h"
#include "asyncde/AsyncIterator.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Problem.h"
#include "asyncde/Rnd.h"
#include "asyncde/Status.h"
#include "asyncde/StopCriteria.h"
#include "asyncde/Variables.h"
#include "asyncde/mathextra.h"

#define messages_stream stderr
#define errors_stream stderr

// #define LC_DEBUG

asyncde::ADERand1Iterator::ADERand1Iterator(const Problem &_problem,
                                            const IteratorConfig *_cfg)
    : AsyncIterator(_problem, _cfg), aderand1cfg(nullptr),
      base_population_actual_size(0), bestpointcpop(0),
      lastidbeforerestart(PointInfo::NOT_INITIALIZED_POINT) {
  const Variables *vars = problem->IntVariables();
  unsigned int nfreeparams = 0;
  unsigned int nfreevars = (vars) ? vars->NVariables() : 0;

  if (nfreevars > 0)
    for (unsigned int ivar = 0; ivar < nfreevars; ivar++)
      if (vars->IsFreeVariable(ivar) > 0)
        nfreeparams++;

  if (nfreeparams < 1 || (nfreeparams != nfreevars)) {
    fprintf(errors_stream,
            "Error in asyncde::ADERand1Iterator::ADERand1Iterator(): "
            "nfreeparams(%i) != "
            "nfreevars(%i)\n",
            nfreeparams, nfreevars);
    statusstopbits |= CRITERION_error;
  }

  const ADERand1Config *_aderand1cfg =
      dynamic_cast<const ADERand1Config *>(_cfg);
  aderand1cfg = (_aderand1cfg) ? new ADERand1Config(*_aderand1cfg)
                               : new ADERand1Config(*_cfg);
  cfg = aderand1cfg;

  // increase population size to minimal number of points required by algorithm
  unsigned int min_nparents = aderand1cfg->MinimalADEPopulationSize();
  if (aderand1cfg->MinPopSize() < min_nparents)
    aderand1cfg->SetMinPopSize(min_nparents);

  if (nfreeparams > 0) {
    parlow = *vars->InitialLowerLimits();
    parup = *vars->InitialUpperLimits();
    parlowindex.resize(nfreeparams, (int)PointInfo::NOT_INITIALIZED_POINT);
    parupindex.resize(nfreeparams, (int)PointInfo::NOT_INITIALIZED_POINT);
  }

  xexttmpvector.resize(problem->NVariables());

  status->popsizemax = aderand1cfg->MinPopSize();
  base_population.resize(status->popsizemax, 0);
  for (Point *&pointptr : base_population) {
    pointptr =
        new Point(*new PointData(*problem->IntVariables(), *problem->Y()),
                  new PointInfo());
    pointptr->Reset();
  }

  bestintpoint = new Point(
      *new PointData(*problem->IntVariables(), *problem->Y()), new PointInfo());
  bestintpoint->Reset();
}

asyncde::ADERand1Iterator::~ADERand1Iterator() {
  for (Point *pointptr : base_population)
    delete pointptr;
}

int asyncde::ADERand1Iterator::UpdateStatus() {
  statusstopbits |= status->FindStatus(cfg->criteriastop);
  statusrestartbits |= status->FindStatus(cfg->criteriarestart);
  statusrecoverbits |= status->FindStatus(cfg->criteriarecover);

  int retvalue = -1;

  if (statusstopbits == 0) {
    if (statusrecoverbits)
      retvalue = Recovery();

    if (0 != retvalue && statusrestartbits)
      retvalue = Restart(2 * aderand1cfg->MinPopSize());

    if (0 != retvalue && (0 != statusrecoverbits || 0 != statusrestartbits))
      statusstopbits = CRITERION_error; // enforce stop;
  }

  return statusstopbits != 0;
}

int asyncde::ADERand1Iterator::lcResize(unsigned int _nparents) {
  if (problem->NFreeVariables() < 1 || problem->Y()->size() < 1)
    return -2;

  unsigned int oldsize = base_population.size();
  if (_nparents < oldsize) {
    for (unsigned int ip = _nparents; ip < oldsize; ip++)
      delete base_population[ip];
    base_population.resize(_nparents, 0);
  } else {
    base_population.resize(_nparents, 0);
    for (unsigned int ip = oldsize; ip < _nparents; ip++)
      base_population[ip] = NewIntPoint();
  }

  for (Point *pointptr : base_population)
    pointptr->Reset();
  base_population_actual_size = 0;

  /*
  fprintf(stderr, "ADERand1Iterator::lcResize():\n");
  for (std::vector<Point *>::iterator it = base_population.begin();
       it != base_population.end(); it++)
    (*it)->Print();
  */

  return 0;
}

int asyncde::ADERand1Iterator::IsContainsX(const std::vector<double> &x,
                                           long int &_id,
                                           Point **pointptr) const {
  _id = PointInfo::NOT_INITIALIZED_POINT;
  if (pointptr)
    *pointptr = 0;

  const double tolerance = problem->Tolerance();
  const unsigned int nfreeparams = problem->NFreeVariables();
  for (unsigned int ip = 0; ip < base_population_actual_size; ip++) {
    int same = 1;
    const std::vector<double> *_point_x = base_population[ip]->Data()->X();
    for (unsigned int idim = 0; idim < nfreeparams; idim++)
      if (0 != fuzzy_cmp(x[idim], (*_point_x)[idim], tolerance)) {
        same = 0;
        break;
      }

    if (same) {
      _id = base_population[ip]->Info()->Id();
      if (pointptr)
        *pointptr = base_population[ip];
      return 1;
    }
  }

  return 0;
}

int asyncde::ADERand1Iterator::FillInTrialIntPoint(Point &_point,
                                                   int maxnrestarts) {
  int retvalue = 0;

  if (_point.Data()->NVariables() != problem->NFreeVariables())
    return -1;

  if (base_population_actual_size < aderand1cfg->MinPopSize()) {
    // generate a random (initial) point
    retvalue = GenerateUniformRandomPoint(_point);
  } else {
    // generate a point according to an ADE strategy
    retvalue = FillInTrialIntADEPoint(_point);
  }

  if (_point.Info()->Id() == PointInfo::NOT_INITIALIZED_POINT ||
      !(_point.Info()->Status() == POINTSTATUS_INITIALIZED ||
        _point.Info()->Status() == POINTSTATUS_REQUESTED)) {
    if (maxnrestarts > 0) {
      maxnrestarts--;

      UpdateStatus();
      retvalue = FillInTrialIntPoint(_point, maxnrestarts);
    } else {
      fprintf(errors_stream,
              "error in asyncde::ADERand1Iterator::FillInTrialPoint()\n");
      _point.Print(errors_stream);
      status->error = 1;
      statusstopbits |= CRITERION_error;
    }
  }

  return retvalue;
}

int asyncde::ADERand1Iterator::GenerateUniformRandomPoint(Point &_point) {
  const std::vector<double> *lowerbound =
      problem->IntVariables()->InitialLowerLimits();
  const std::vector<double> *upperbound =
      problem->IntVariables()->InitialUpperLimits();

  _point.Reset();

  // check whether all lower and upper bounds are specified
  if (!lowerbound || !upperbound)
    return -1;

  std::vector<double> *x = _point.DataMutable()->Xmutable();

  status->nprojfeas = 0;
  unsigned int nmax_attempts = std::numeric_limits<unsigned int>::max();
  if ((cfg->criteriastop.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriastop.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriastop.nprojfeas;
  if ((cfg->criteriarestart.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriarestart.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriarestart.nprojfeas;
  long int tag;

  // cycle to test that new point is unique in the population
  const unsigned int nfreeparams = problem->NFreeVariables();
  do {
    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
      (*x)[ivar] = cfg->rnd->next((*lowerbound)[ivar], (*upperbound)[ivar]);
    _point.InfoMutable()->SetStatus(POINTSTATUS_INITIALIZED);

    tag = PointInfo::NOT_INITIALIZED_POINT;
    IsContainsX(*x, tag);

    if (1 != problem->IsXIntFeasible(*x, xexttmpvector))
      tag = 0;
    status->Incr_nProjFeas();
  } while ((tag != PointInfo::NOT_INITIALIZED_POINT) &&
           status->nprojfeas < nmax_attempts);

  if (tag == PointInfo::NOT_INITIALIZED_POINT)
    _point.InfoMutable()->SwitchToNewId();
  else
    _point.InfoMutable()->SetStatus(POINTSTATUS_INFEASIBLE);

#ifdef LC_DEBUG
  if (cfg->verbose >= 1) {
    fprintf(messages_stream,
            "asyncde::ADERand1Iterator::GenerateUniformRandomPoint():\n");
    fprintf(messages_stream, "new point:\n");
    _point.Print(messages_stream);

    if (cfg->verbose > 2) {
      fprintf(messages_stream, "Initial lower bounds:\n");
      for (double xlow : *lowerbound)
        fprintf(messages_stream, "%.3e  ", xlow);
      fprintf(messages_stream, "\nInitial upper bounds:\n");
      for (double xup : *upperbound)
        fprintf(messages_stream, "%.3e  ", xup);
      fprintf(messages_stream, "\n");
    }
  }
#endif // LC_DEBUG

  return 0;
}

int asyncde::ADERand1Iterator::FillInTrialIntADEPoint(Point &_point) {
  _point.Reset();

  if (!base_population[aderand1cfg->MinPopSize() - 1] ||
      base_population[aderand1cfg->MinPopSize() - 1]->Info()->Id() ==
          PointInfo::NOT_INITIALIZED_POINT)
    return -1;

  std::vector<double> *x = _point.DataMutable()->Xmutable();

  status->nprojfeas = 0;
  unsigned int nmax_attempts = std::numeric_limits<unsigned int>::max();
  if ((cfg->criteriastop.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriastop.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriastop.nprojfeas;
  if ((cfg->criteriarestart.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriarestart.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriarestart.nprojfeas;
  long int tag = PointInfo::NOT_INITIALIZED_POINT;

  const unsigned int nfreevars = problem->NFreeVariables();
  do {
    // Generate a mutant vector
    cfg->rnd->randdistinct(3, base_population_actual_size, nrand, nrand_sorted);

    const std::vector<double> *x1 = base_population[nrand[0]]->Data()->X();
    const std::vector<double> *x2 = base_population[nrand[1]]->Data()->X();
    const std::vector<double> *x3 = base_population[nrand[2]]->Data()->X();
    const double F = aderand1cfg->GetF();
    int similarvar = -1;
    for (unsigned int ix = 0; ix < nfreevars; ix++) {
      (*x)[ix] = (*x1)[ix] + F * ((*x2)[ix] - (*x3)[ix]);
      if (similarvar == -1 &&
          0 == fuzzy_cmp((*x2)[ix], (*x3)[ix], cfg->criteriarestart.xepsilon))
        similarvar = ix;
    }

    if (similarvar != -1)
      status->Incr_nxdiff();

    _point.InfoMutable()->SetStatus(POINTSTATUS_INITIALIZED);

    tag = PointInfo::NOT_INITIALIZED_POINT;

    // project the mutant vector into the feasible region
    if (1 != problem->IsXIntFeasible(*x, xexttmpvector)) {
      if (0 != ProjectPoint2Feasibility(*x, *x1, xexttmpvector))
        tag = 0;
    }

    // check whether the mutant vector doesn't coincide with any vector within
    // population
    //    IsContainsX(*x, tag);
    status->Incr_nProjFeas();
  } while ((tag != PointInfo::NOT_INITIALIZED_POINT) &&
           status->nprojfeas < nmax_attempts);

  if (tag == PointInfo::NOT_INITIALIZED_POINT)
    _point.InfoMutable()->SwitchToNewId();
  else
    _point.InfoMutable()->SetStatus(POINTSTATUS_INFEASIBLE);

  return 0;
}

void asyncde::ADERand1Iterator::UpdatePopulationSpread(
    const long int newpoint_index) {
  unsigned int ip;
  double vmin, vmax, vnew;
  int population_complete =
      base_population_actual_size == aderand1cfg->MinPopSize();
  const unsigned int nfreeparams = problem->NFreeVariables();

  status->ResetXYEpsilons();

  if (bestpointcpop) {
    const std::vector<double> *newvX =
        base_population[newpoint_index]->Data()->X();

    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++) {
      int needcycle = 0;
      vnew = (*newvX)[ivar];
      vmin = std::numeric_limits<double>::max();
      vmax = -std::numeric_limits<double>::max();

      // range from the previous iteration
      if (parlowindex[ivar] > PointInfo::NOT_INITIALIZED_POINT)
        vmin = parlow[ivar];
      if (parupindex[ivar] > PointInfo::NOT_INITIALIZED_POINT)
        vmax = parup[ivar];

      if (vnew < vmin) {
        vmin = vnew;
        parlow[ivar] = vmin;
        parlowindex[ivar] = ivar;
      } else if (parlowindex[ivar] == (int)ivar)
        needcycle = 1;

      if (vnew > vmax) {
        vmax = vnew;
        parup[ivar] = vmax;
        parupindex[ivar] = ivar;
      } else if (parupindex[ivar] == (int)ivar)
        needcycle = 1;

      if (needcycle) {
        // analyse spread in the population
        vmin = std::numeric_limits<double>::max();
        parlowindex[ivar] = PointInfo::NOT_INITIALIZED_POINT;
        vmax = -std::numeric_limits<double>::max();
        parupindex[ivar] = PointInfo::NOT_INITIALIZED_POINT;
        for (ip = 0; ip < base_population_actual_size; ip++) {
          vnew = (*base_population[ip]->Data()->X())[ivar];
          if (vnew < vmin) {
            vmin = vnew;
            parlowindex[ivar] = ip;
          }
          if (vnew > vmax) {
            vmax = vnew;
            parupindex[ivar] = ip;
          }
        }
        parlow[ivar] = vmin;
        parup[ivar] = vmax;
      }

      if (population_complete)
        status->AddMinMaxPair(vmin, vmax, 1);
      else
        status->AddMinMaxPair(-std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::max(), 1);
      /*
            // debug spread in population:
            if (vmin != parlow[ivar])
              fprintf(errors_stream,
                      "error in
         asyncde::ADERand1Iterator::UpdatePopulationSpread(): " "vmin !=
         parlow[%i]\n", ivar); if (vmax != parup[ivar]) fprintf(errors_stream,
                      "error in
         asyncde::ADERand1Iterator::UpdatePopulationSpread(): " "vmax !=
         parup[%i]\n", ivar); for (ip = 0; ip < base_population_actual_size;
         ip++) { vnew = (*base_population[ip]->Data()->X())[ivar]; if (vnew <
         parlow[ivar]) fprintf( errors_stream, "error in
         asyncde::ADERand1Iterator::UpdatePopulationSpread(): " "vnew <
         parlow[%i]\n", ivar); if (vnew > parup[ivar]) fprintf( errors_stream,
                    "error in
         asyncde::ADERand1Iterator::UpdatePopulationSpread(): " "vnew >
         parup[%i]\n", ivar);
            }
      */
    }

    // update Y
    vmin = (*base_population[0]->Data()->Y())[0];
    vmax = vmin;

    for (ip = 1; ip < base_population_actual_size; ip++) {
      vnew = (*base_population[ip]->Data()->Y())[0];
      if (vnew > vmax)
        vmax = vnew;
      else if (vnew < vmin)
        vmin = vnew;
    }

    if (population_complete)
      status->AddMinMaxPair(vmin, vmax, 0);
    else
      status->AddMinMaxPair(-std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max(), 0);
  }
}

int asyncde::ADERand1Iterator::AddIntPoint(Point &_point) {
  if (_point.Info()->Status() < POINTSTATUS_EVALUATED)
    return -2;

  if ((*_point.Data()->Y())[0] < (*bestintpoint->Data()->Y())[0]) {
    bestintpoint->Set(_point);
    status->vtr = (*bestintpoint->Data()->Y())[0];
  }

  // add points generated before restart
  if (_point.Info()->Id() <= lastidbeforerestart) {
    status->Incr_nFE();
    return 0;
  }

  int selected = 0;
  int ipos = -1;
  if (base_population_actual_size < aderand1cfg->MinPopSize()) {
    selected = 1;
    ipos = base_population_actual_size;
    base_population[ipos]->Set(_point);
    base_population_actual_size++;
  } else {
    // compare to a random population member
    ipos = cfg->rnd->next_ulong(base_population_actual_size);
    if ((*_point.Data()->Y())[0] < (*base_population[ipos]->Data()->Y())[0]) {
      selected = 1;
      base_population[ipos]->Set(_point);
    }
  }

  status->Incr_nFE();

  if (selected) {
    if (!bestpointcpop || (*base_population[ipos]->Data()->Y())[0] <
                              (*bestpointcpop->Data()->Y())[0]) {
      bestpointcpop = base_population[ipos];
      status->NewBestFound(base_population_actual_size ==
                           aderand1cfg->MinPopSize());
      /*
      fprintf(stderr, "%li: Improvement -> %.15e (%.15e)\n",
              status->nFE, (*bestpointcpop->Data()->Y())[0], status->vtr);
      */
    }
    UpdatePopulationSpread(ipos);
  }

  UpdateStatus();

  return selected;
}

unsigned int
asyncde::ADERand1Iterator::PopulationSize(int *restartcounter) const {
  if (restartcounter)
    *restartcounter = statusrestartcounter;

  return aderand1cfg->MinPopSize();
}

int asyncde::ADERand1Iterator::Restart(unsigned int _nparents) {
  if (problem->NFreeVariables() < 1 || problem->Y()->size() < 1)
    return -2;

  if (bestpointcpop &&
      bestpointcpop->Info()->Id() != PointInfo::NOT_INITIALIZED_POINT)
    lastidbeforerestart = bestpointcpop->Info()->IdCounter();

  // refine _nparents
  if ((cfg->criteriastop.statusbits & CRITERION_popsizemax) ==
          CRITERION_popsizemax &&
      cfg->criteriastop.popsizemax > 0 &&
      (int)_nparents > cfg->criteriastop.popsizemax) {
    statusstopbits |= CRITERION_popsizemax;
    return -4;
  }
  if ((cfg->criteriarestart.statusbits & CRITERION_popsizemax) ==
          CRITERION_popsizemax &&
      cfg->criteriarestart.popsizemax > 0 &&
      (int)_nparents > aderand1cfg->criteriarestart.popsizemax) {
    _nparents = aderand1cfg->criteriarestart.popsizemax;
    statusrestartcounter++;
  }

  if (cfg->verbose) {
    fprintf(messages_stream,
            "restart(): mask=0x%x  nFE=%li  new_size=%i  "
            "value[id=%li]=%.15e \n",
            statusrestartbits, status->nFE, _nparents,
            bestpointcpop->Info()->Id(), (*bestpointcpop->Data()->Y())[0]);
    //    bestpointcpop->Print(messages_stream);
  }

  const Variables *intvars = problem->IntVariables();
  parlow = *intvars->InitialLowerLimits();
  parup = *intvars->InitialUpperLimits();
  std::fill(parlowindex.begin(), parlowindex.end(),
            PointInfo::NOT_INITIALIZED_POINT);
  std::fill(parupindex.begin(), parupindex.end(),
            PointInfo::NOT_INITIALIZED_POINT);

  lcResize(_nparents);

  aderand1cfg->SetMinPopSize(_nparents);
  status->popsizemax = _nparents;

  status->ResetAfterRestart();
  statusrestartbits = 0;
  statusrecoverbits = 0;

  return 0;
}

int asyncde::ADERand1Iterator::Recovery() {
  if (cfg->verbose)
    fprintf(messages_stream, "recovery(): mask=0x%x  nFE=%li\n",
            statusrecoverbits, status->nFE);

  if (status->nrecover >= cfg->criteriarecover.nrecover) {
    if (cfg->verbose)
      fprintf(messages_stream,
              "Recovery(): not implemented due to the maximal number of "
              "attempts (%i)\n",
              cfg->criteriarecover.nrecover);
    return -1;
  }

  int nrec = status->nrecover;
  nrec++;
  int retvalue = Restart(aderand1cfg->MinPopSize());
  status->nrecover = nrec;
  if (retvalue < 0)
    return retvalue;

  if (bestintpoint &&
      bestintpoint->Info()->Id() != PointInfo::NOT_INITIALIZED_POINT) {
    if (!tmpintpoint)
      tmpintpoint = NewIntPoint();
    tmpintpoint->Set(*bestintpoint);
    tmpintpoint->InfoMutable()->SwitchToNewId();
    AddIntPoint(*tmpintpoint);
  }

  return 0;
}
