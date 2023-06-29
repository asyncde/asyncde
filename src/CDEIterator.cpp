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

/**
 CDEIterator implements asynchronous variants
 of (Classic) Differential Evolution with the binomial crossover
 */

#include <algorithm>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <utility>

#include "asyncde/CDEIterator.h"

#include "asyncde/ADEArchive.h"
#include "asyncde/ADEConfig.h"
#include "asyncde/ADEPoint.h"
#include "asyncde/ADEPointInfo.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Problem.h"
#include "asyncde/Rnd.h"
#include "asyncde/Status.h"
#include "asyncde/StopCriteria.h"
#include "asyncde/Variables.h"
#include "asyncde/mathextra.h"

//#define LC_DEBUG

//#define messages_stream stdout
#define messages_stream stderr
#define errors_stream stderr

asyncde::CDEIterator::CDEIterator(const Problem &_problem,
                                  const IteratorConfig *_cfg)
    : AsyncIterator(_problem, _cfg), adecfg(nullptr), nfreeparams(0),
      base_population_actual_size(0),
      lastidbeforerestart(PointInfo::NOT_INITIALIZED_POINT), archive(nullptr),
      statsumsupdates_counter(0), statsumsupdates_counter_max(1000), Fsum(0.0),
      F2sum(0.0), FmuJADE(0.5), CRn(0), CRsum(0.0), CR2sum(0.0), CRmu(0.0) {
  const Variables *vars = problem->IntVariables();
  unsigned int nfreevars = (vars) ? vars->NVariables() : 0;

  if (nfreevars > 0)
    for (unsigned int ivar = 0; ivar < nfreevars; ivar++)
      if (vars->IsFreeVariable(ivar) > 0)
        nfreeparams++;

  if (nfreeparams < 1 || (nfreeparams != nfreevars)) {
    fprintf(errors_stream,
            "Error in asyncde::CDEIterator::CDEIterator(): nfreeparams(%i) != "
            "nfreevars(%i)\n",
            nfreeparams, nfreevars);
    statusstopbits |= CRITERION_error;
  }

  const ADEConfig *_adecfg = dynamic_cast<const ADEConfig *>(_cfg);
  adecfg = (_adecfg) ? new ADEConfig(*_adecfg)
                     : new ADEConfig(_cfg->rnd, nfreeparams);
  cfg = adecfg;

  // increase population size to minimal number of points required by algorithm
  unsigned int min_nparents = adecfg->MinimalADEPopulationSize();
  if (adecfg->MinPopSize() < min_nparents)
    adecfg->SetMinPopSize(min_nparents);

  if (adecfg->archive)
    archive = new ADEArchive(adecfg->MinPopSize(), *adecfg);

  if (nfreeparams > 0) {
    parlow = *vars->InitialLowerLimits();
    parup = *vars->InitialUpperLimits();
    const long int piddefault = PointInfo::NOT_INITIALIZED_POINT;
    parlowpid.resize(nfreeparams, piddefault);
    paruppid.resize(nfreeparams, piddefault);

    CDEIterator::FCrDefaultSettings();
  }

  xexttmpvector.resize(problem->NVariables());

  CDEIterator::lcResize(adecfg->MinPopSize());
  status->popsizemax = adecfg->MinPopSize();

  bestintpoint = new Point(*base_population[0]->Data()->Clone(),
                           base_population[0]->Info()->Clone());
}

asyncde::CDEIterator::~CDEIterator() {
  for (ADEPoint *pointptr : base_population)
    delete pointptr;

  delete archive;
}

asyncde::Point *asyncde::CDEIterator::NewExtPoint() const {
  ADEPoint *newpoint =
      new ADEPoint(*new PointData(*problem->ExtVariables(), *problem->Y()),
                   new ADEPointInfo(nfreeparams));
  newpoint->Reset();
  return newpoint;
}

asyncde::Point *asyncde::CDEIterator::NewIntPoint() const {
  ADEPoint *newpoint =
      new ADEPoint(*new PointData(*problem->IntVariables(), *problem->Y()),
                   new ADEPointInfo(nfreeparams));
  newpoint->Reset();
  return newpoint;
}

unsigned int asyncde::CDEIterator::PopulationSize(int *restartcounter) const {
  if (restartcounter)
    *restartcounter = statusrestartcounter;

  return adecfg->MinPopSize();
}

int asyncde::CDEIterator::UpdateStatus() {
  statusstopbits |= status->FindStatus(cfg->criteriastop);
  statusrestartbits |= status->FindStatus(cfg->criteriarestart);
  statusrecoverbits |= status->FindStatus(cfg->criteriarecover);

  int retvalue = -1;

  if (statusstopbits == 0) {
    if (statusrecoverbits)
      retvalue = Recovery();

    if (0 != retvalue && statusrestartbits)
      retvalue = Restart(2 * adecfg->MinPopSize());

    if (0 != retvalue && (0 != statusrecoverbits || 0 != statusrestartbits))
      statusstopbits = CRITERION_error; // enforce stop;
  }

  return statusstopbits != 0;
}

int asyncde::CDEIterator::lcResize(unsigned int _nparents) {
  if (nfreeparams < 1 || problem->Y()->size() < 1)
    return -2;

  vetovector.reserve(_nparents);

  unsigned int oldsize = base_population.size();
  base_population.resize(_nparents, 0);
  for (unsigned int ip = oldsize; ip < _nparents; ip++)
    base_population[ip] =
        new ADEPoint(*new PointData(*problem->IntVariables(), *problem->Y()),
                     new ADEPointInfo(nfreeparams));
  for (ADEPoint *pointptr : base_population)
    pointptr->Reset();
  base_population_actual_size = 0;
  std::fill(pid_key.begin(), pid_key.end(),
            (long int)PointInfo::NOT_INITIALIZED_POINT);
  pid_key.resize(_nparents, (long int)PointInfo::NOT_INITIALIZED_POINT);
  std::fill(value_key.begin(), value_key.end(),
            std::numeric_limits<double>::max());
  value_key.resize(_nparents, std::numeric_limits<double>::max());
  std::fill(pid_index_for_value_key.begin(), pid_index_for_value_key.end(),
            (int)PointInfo::NOT_INITIALIZED_POINT);
  pid_index_for_value_key.resize(_nparents,
                                 (int)PointInfo::NOT_INITIALIZED_POINT);
  std::fill(value_index_for_pid_key.begin(), value_index_for_pid_key.end(),
            (int)PointInfo::NOT_INITIALIZED_POINT);
  value_index_for_pid_key.resize(_nparents,
                                 (int)PointInfo::NOT_INITIALIZED_POINT);

  if (archive) {
    archive->Resize(_nparents);
    archive->Reset();
  }

  return 0;
}

int asyncde::CDEIterator::Restart(unsigned int _nparents) {
  if (nfreeparams < 1 || problem->Y()->size() < 1)
    return -2;

  const ADEPoint *bestpoint = BestIntPointCurrentPopulation();
  if (bestpoint && bestpoint->Info()->Id() != PointInfo::NOT_INITIALIZED_POINT)
    lastidbeforerestart = bestpoint->Info()->IdCounter();

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
      (int)_nparents > adecfg->criteriarestart.popsizemax) {
    _nparents = adecfg->criteriarestart.popsizemax;
    statusrestartcounter++;
  }

  if (cfg->verbose) {
    fprintf(messages_stream,
            "restart(): mask=0x%x  nFE=%li  new_size=%i  "
            "value[id=%li,pid=%li]=%.15e  nxdiff=%i\n",
            statusrestartbits, status->nFE, _nparents, bestpoint->Info()->Id(),
            bestpoint->ADEInfo()->ParentId(), (*bestpoint->Data()->Y())[0],
            status->nxdiff);
    if (cfg->verbose > 1)
      bestpoint->Print(messages_stream);
#ifdef LC_DEBUG
      /*
        if (cfg->verbose >= 2) {
          long int nentries;
          double mean, stddev, meanerr;
          if ((statusrestartbits & CRITERION_nFE_after_best) != 0) {
            nentries = status->nFE_after_best_history->get_n_entries();
            status->nFE_after_best_history->find_stat_properties(mean, stddev,
                                                                 &meanerr);
            fprintf(messages_stream,
                    "nFE_after_best = %li   nFE_after_best_max = %.2e  [meanerr
        "
                    "= %.2e]  sigma = %.2e  nentries = %lu\n",
                    status->nFE_after_best,
                    cfg->criteriarestart.nFE_after_best_factor * mean, meanerr,
                    stddev, nentries);
          }

          if ((statusrestartbits & CRITERION_nFE_after_progress) != 0) {
            nentries = status->nFE_after_progress_history->get_n_entries();
            status->nFE_after_progress_history->find_stat_properties(mean,
        stddev, &meanerr); fprintf(messages_stream, "nFE_after_progress = %li
        nFE_after_progress_max = %.2e  "
                    "[meanerr = %.2e]  sigma = %.2e  nentries = %lu\n",
                    status->nFE_after_progress,
                    cfg->criteriarestart.nFE_after_progress_factor * mean,
        meanerr, stddev, nentries);
          }
        }
    */
#endif // LC_DEBUG
  }

  const Variables *intvars = problem->IntVariables();
  parlow = *intvars->InitialLowerLimits();
  parup = *intvars->InitialUpperLimits();
  std::fill(parlowpid.begin(), parlowpid.end(),
            (long int)PointInfo::NOT_INITIALIZED_POINT);
  std::fill(paruppid.begin(), paruppid.end(),
            (long int)PointInfo::NOT_INITIALIZED_POINT);

  lcResize(_nparents);

  adecfg->SetMinPopSize(_nparents);
  status->popsizemax = _nparents;

  status->ResetAfterRestart();
  statusrestartbits = 0;
  statusrecoverbits = 0;

  FCrDefaultSettings();

  return 0;
}

int asyncde::CDEIterator::Recovery() {
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
  int retvalue = Restart(adecfg->MinPopSize());
  status->nrecover = nrec;
  if (retvalue < 0)
    return retvalue;

  const Point *bestpoint = BestIntPoint();
  if (bestpoint &&
      bestpoint->Info()->Id() != PointInfo::NOT_INITIALIZED_POINT) {
    if (!tmpintpoint)
      tmpintpoint = NewIntPoint();
    tmpintpoint->Set(*bestpoint);
    tmpintpoint->InfoMutable()->SwitchToNewId();
    AddIntPoint(*tmpintpoint);
  }

  if (archive)
    archive->Reset();

  return 0;
}

const asyncde::ADEPoint *
asyncde::CDEIterator::BestIntPointCurrentPopulation() const {
  if (pid_index_for_value_key.size() < 1 || pid_index_for_value_key[0] < 0 ||
      pid_index_for_value_key[0] >= (int)adecfg->MinPopSize())
    return 0;

  const ADEPoint *bestpoint = base_population[pid_index_for_value_key[0]];
  if (!bestpoint || bestpoint->Info()->Id() == PointInfo::NOT_INITIALIZED_POINT)
    return 0;

  return bestpoint;
}

int asyncde::CDEIterator::FCrDefaultSettings() {
  Fsum = 0.0;
  F2sum = 0.0;
  FmuJADE = adecfg->Fmu;

  CRn = 0;
  CRsum = 0.0;
  CR2sum = 0.0;
  CRmu = adecfg->CRmu;

  return 0;
}

void asyncde::CDEIterator::UpdatePopulationSpread(const long int iposition) {
  double vmin, vmax, vnew;
  int population_complete = base_population_actual_size == adecfg->MinPopSize();

  status->ResetXYEpsilons();

  if (base_population_actual_size > 0) {
    const long int newpoint_parentid = pid_key[iposition];
    const std::vector<double> *newvX = base_population[iposition]->Data()->X();

    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++) {
      int needcycle = 0;
      vnew = (*newvX)[ivar];
      vmin = std::numeric_limits<double>::max();
      vmax = -std::numeric_limits<double>::max();

      // range from the previous iteration
      if (parlowpid[ivar] > PointInfo::NOT_INITIALIZED_POINT)
        vmin = parlow[ivar];
      if (paruppid[ivar] > PointInfo::NOT_INITIALIZED_POINT)
        vmax = parup[ivar];

      if (vnew < vmin) {
        vmin = vnew;
        parlow[ivar] = vmin;
        parlowpid[ivar] = newpoint_parentid;
      } else if (parlowpid[ivar] == newpoint_parentid)
        needcycle = 1;

      if (vnew > vmax) {
        vmax = vnew;
        parup[ivar] = vmax;
        paruppid[ivar] = newpoint_parentid;
      } else if (paruppid[ivar] == newpoint_parentid)
        needcycle = 1;

      if (needcycle) {
        // analyse spread in the population
        vmin = std::numeric_limits<double>::max();
        parlowpid[ivar] = PointInfo::NOT_INITIALIZED_POINT;
        vmax = -std::numeric_limits<double>::max();
        paruppid[ivar] = PointInfo::NOT_INITIALIZED_POINT;
        for (unsigned int ip = 0; ip < base_population_actual_size; ip++) {
          vnew = (*base_population[ip]->Data()->X())[ivar];
          if (vnew < vmin) {
            vmin = vnew;
            parlowpid[ivar] = pid_key[ip];
          }
          if (vnew > vmax) {
            vmax = vnew;
            paruppid[ivar] = pid_key[ip];
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
         asyncde::CDEIterator::UpdatePopulationSpread(): " "vmin !=
         parlow[%i]\n", ivar); if (vmax != parup[ivar]) fprintf(errors_stream,
                      "error in
         asyncde::CDEIterator::UpdatePopulationSpread(): " "vmax !=
         parup[%i]\n", ivar); for (ip = 0; ip < base_population_actual_size;
         ip++) { vnew = (*base_population[ip]->Data()->X())[ivar]; if (vnew <
         parlow[ivar]) fprintf( errors_stream, "error in
         asyncde::CDEIterator::UpdatePopulationSpread(): " "vnew <
         parlow[%i]\n", ivar); if (vnew > parup[ivar]) fprintf( errors_stream,
                    "error in
         asyncde::CDEIterator::UpdatePopulationSpread(): " "vnew >
         parup[%i]\n", ivar);
            }

            //      fprintf(errors_stream, "(vmax-vmin)=%.3e\n", vmax - vmin);
      */
    }

    // update Y
    vmin = (*base_population[0]->Data()->Y())[0];
    vmax = (*base_population[base_population_actual_size - 1]->Data()->Y())[0];

    if (population_complete)
      status->AddMinMaxPair(vmin, vmax, 0);
    else
      status->AddMinMaxPair(-std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::max(), 0);
  }
}

int asyncde::CDEIterator::Selection(ADEPoint &_adepoint, int &ipidpos2replace) {
  ipidpos2replace = -1;

  if (base_population_actual_size < adecfg->MinPopSize())
    return 1;

  switch (_adepoint.ADEInfo()->selectiontype) {
  case ADE_SELECTION_RANDOM:
    ipidpos2replace = cfg->rnd->next_ulong(base_population_actual_size);
    _adepoint.ADEInfoMutable()->SetParentId(pid_key[ipidpos2replace]);
    break;
  case ADE_SELECTION_WORST:
    ipidpos2replace = pid_index_for_value_key[base_population_actual_size - 1];
    _adepoint.ADEInfoMutable()->SetParentId(pid_key[ipidpos2replace]);
    break;
  case ADE_SELECTION_PARENT:
  default:
    std::vector<long int>::iterator pidpos = std::lower_bound(
        pid_key.begin(), pid_key.begin() + base_population_actual_size,
        _adepoint.ADEInfo()->ParentId());
    if (pidpos < pid_key.end())
      ipidpos2replace = pidpos - pid_key.begin();
    else
      return 0; // parent not found
  }

  return (*_adepoint.Data()->Y())[0] <
         (*base_population[ipidpos2replace]->Data()->Y())[0];
}

int asyncde::CDEIterator::AddIntPoint(Point &_point) {
  if (_point.Info()->Status() < POINTSTATUS_EVALUATED)
    return -2;

  if ((*_point.Data()->Y())[0] < (*bestintpoint->Data()->Y())[0]) {
    bestintpoint->Set(_point);
    status->vtr = (*bestintpoint->Data()->Y())[0];
  }

#ifdef LC_DEBUG
  if (cfg->verbose > 1) {
    fprintf(messages_stream, "asyncde::CDEIterator::AddIntPoint():\n");
    fprintf(messages_stream, "new point:\n");
    _point.Print(messages_stream);
    fprintf(messages_stream, "initial population:\n");
    for (unsigned int ip = 0; ip < base_population_actual_size; ip++)
      base_population[ip]->Print(messages_stream);
    fprintf(messages_stream, "----------\n");
  }
#endif // LC_DEBUG

  if (_point.Info()->Id() <= lastidbeforerestart) {
    status->Incr_nFE();
    return 0;
  }

  int selected = 0;
  ADEPoint *_adepoint = dynamic_cast<ADEPoint *>(&_point);
  if (!_adepoint)
    return -4;

  int iposition = -1;
  int ipos2replace = -1;
  double newY = (*_point.Data()->Y())[0];
  status->Incr_nFE();

  if (0 == base_population_actual_size) {
    iposition = 0;
    selected = 1;

    _adepoint->ADEInfoMutable()->selectiontype = ADE_SELECTION_PARENT;
    base_population[0]->Set(*_adepoint);
    pid_key[0] = _adepoint->ADEInfo()->ParentId();
    value_key[0] = newY;
    pid_index_for_value_key[0] = 0;
    value_index_for_pid_key[0] = 0;
    base_population_actual_size++;
    UpdateStatSums(_adepoint, 0);
  } else {
    selected = Selection(*_adepoint, ipos2replace);

    if (selected) {
      if (ipos2replace < 0) {
        // find position of the new point in the current population
        std::vector<long int>::iterator newpidpos = std::lower_bound(
            pid_key.begin(), pid_key.begin() + base_population_actual_size,
            _adepoint->ADEInfo()->ParentId());
        iposition = newpidpos - pid_key.begin();
      } else {
        iposition = ipos2replace;
      }

      int ioldvaluepos =
          (ipos2replace >= 0) ? value_index_for_pid_key[ipos2replace] : -1;
      std::vector<double>::iterator newvaluepos = std::lower_bound(
          value_key.begin(), value_key.begin() + base_population_actual_size,
          newY);
      int inewvaluepos = newvaluepos - value_key.begin();

      if (iposition >= 0 &&
          base_population_actual_size < adecfg->MinPopSize()) {
        // new entry; re-order previous data
        ADEPoint *tmp_point_ptr = base_population[base_population_actual_size];
        std::move_backward(
            base_population.begin() + iposition,
            base_population.begin() + base_population_actual_size,
            base_population.begin() + base_population_actual_size + 1);
        base_population[iposition] = tmp_point_ptr;
        std::move_backward(pid_key.begin() + iposition,
                           pid_key.begin() + base_population_actual_size,
                           pid_key.begin() + base_population_actual_size + 1);
        std::move_backward(
            value_index_for_pid_key.begin() + iposition,
            value_index_for_pid_key.begin() + base_population_actual_size,
            value_index_for_pid_key.begin() + base_population_actual_size + 1);
        for (unsigned int ipos = iposition + 1;
             ipos <= base_population_actual_size; ipos++)
          pid_index_for_value_key[value_index_for_pid_key[ipos]] = ipos;
      }
      if (ioldvaluepos >= 0) {
        // replace old value
        if (inewvaluepos < ioldvaluepos) {
          std::move_backward(value_key.begin() + inewvaluepos,
                             value_key.begin() + ioldvaluepos,
                             value_key.begin() + ioldvaluepos + 1);
          std::move_backward(pid_index_for_value_key.begin() + inewvaluepos,
                             pid_index_for_value_key.begin() + ioldvaluepos,
                             pid_index_for_value_key.begin() + ioldvaluepos +
                                 1);
          for (unsigned int ipos = inewvaluepos + 1; (int)ipos <= ioldvaluepos;
               ipos++)
            value_index_for_pid_key[pid_index_for_value_key[ipos]] = ipos;
        }
      } else {
        if (inewvaluepos < (int)base_population_actual_size &&
            base_population_actual_size < adecfg->MinPopSize()) {
          // new entry
          std::move_backward(value_key.begin() + inewvaluepos,
                             value_key.begin() + base_population_actual_size,
                             value_key.begin() + base_population_actual_size +
                                 1);
          std::move_backward(pid_index_for_value_key.begin() + inewvaluepos,
                             pid_index_for_value_key.begin() +
                                 base_population_actual_size,
                             pid_index_for_value_key.begin() +
                                 base_population_actual_size + 1);
          for (unsigned int ipos = inewvaluepos + 1;
               ipos <= base_population_actual_size; ipos++)
            value_index_for_pid_key[pid_index_for_value_key[ipos]] = ipos;
        }
      }

      // update statistical and adaptive info
      UpdateStatSums(_adepoint,
                     (ipos2replace >= 0) ? base_population[ipos2replace] : 0);

      // set new values
      base_population[iposition]->Set(*_adepoint);
      pid_key[iposition] = _adepoint->ADEInfo()->ParentId();
      value_key[inewvaluepos] = newY;
      pid_index_for_value_key[inewvaluepos] = iposition;
      value_index_for_pid_key[iposition] = inewvaluepos;

      if (ipos2replace < 0)
        base_population_actual_size++;
    } // if (selected)
  }

  if (selected) {
    if (_adepoint->Info()->Id() ==
        BestIntPointCurrentPopulation()->Info()->Id())
      status->NewBestFound(base_population_actual_size == adecfg->MinPopSize());
    else
      status->ProgressFound(base_population_actual_size ==
                            adecfg->MinPopSize());

    UpdatePopulationSpread(iposition);

#ifdef LC_DEBUG
    if (cfg->verbose) {
      if (_adepoint->Info()->Id() ==
          BestIntPointCurrentPopulation()->Info()->Id()) {
        fprintf(messages_stream, "%li: New best point found:\n", status->nFE);
        _adepoint->Print(messages_stream);
      }
      if (cfg->verbose > 1) {
        fprintf(messages_stream, "New population (PID ordered):\n");
        for (unsigned int ip = 0; ip < base_population_actual_size; ip++)
          if (base_population[ip])
            base_population[ip]->Print(messages_stream);
        fprintf(messages_stream, "New population (value ordered):\n");
        for (unsigned int ip = 0; ip < base_population_actual_size; ip++) {
          fprintf(messages_stream, "pid_index_for_value_key[%i]=%i\n", ip,
                  pid_index_for_value_key[ip]);
          if (base_population[pid_index_for_value_key[ip]])
            base_population[pid_index_for_value_key[ip]]->Print(
                messages_stream);
        }
        fprintf(messages_stream, "----------\n");
      }
    }
#endif // LC_DEBUG
  }    // if (selected)

  UpdateStatus();

#ifdef LC_DEBUG
  double ybest = (*base_population[pid_index_for_value_key[0]]->Data()->Y())[0];
  fprintf(stdout, "%.5e [%.5e]:", ybest, (*BestIntPoint()->Data()->Y())[0]);
  for (unsigned int ip = 1; ip < base_population_actual_size; ip++)
    fprintf(stdout, " %.3e",
            (*base_population[pid_index_for_value_key[ip]]->Data()->Y())[0] -
                ybest);
  fprintf(stdout, "\n");
#endif // LC_DEBUG

  return selected;
}

int asyncde::CDEIterator::IsContainsX(const std::vector<double> &x,
                                      long int &_id, Point **pointptr) const {
  _id = PointInfo::NOT_INITIALIZED_POINT;
  if (pointptr)
    *pointptr = 0;

  const double tolerance = problem->Tolerance();
  for (unsigned int ip = 0; ip < base_population_actual_size; ip++) {
    /*
        if (!base_population[ip] ||
       base_population[ip]->Info()->Id() ==
       PointInfo::NOT_INITIALIZED_POINT) { fprintf(errors_stream, "error in
       asyncde::CDEIterator::containsX()\n"); break;
        }
    */
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

int asyncde::CDEIterator::FillInTrialIntPoint(Point &_point, int maxnrestarts) {
  int retvalue = 0;

  if (_point.Data()->NVariables() != nfreeparams)
    return -1;

  ADEPoint *adepoint = dynamic_cast<ADEPoint *>(&_point);
  if (!adepoint)
    return -2;

  if (base_population_actual_size < adecfg->MinPopSize()) {
    // generate a random (initial) point
    retvalue = GenerateUniformRandomPoint(*adepoint);
  } else {
    // generate a point according to an ADE strategy
    retvalue = FillInTrialIntADEPoint(*adepoint);
  }

  if (_point.Info()->Id() == PointInfo::NOT_INITIALIZED_POINT ||
      !(_point.Info()->Status() == POINTSTATUS_INITIALIZED ||
        _point.Info()->Status() == POINTSTATUS_REQUESTED)) {
    if (maxnrestarts > 0) {
      maxnrestarts--;

      fprintf(stderr,
              "recursive call of FillInTrialIntPoint(), maxnrestarts=%i\n",
              maxnrestarts);

      UpdateStatus();
      retvalue = FillInTrialIntPoint(_point, maxnrestarts);
    } else {
      fprintf(errors_stream,
              "error in asyncde::CDEIterator::FillInTrialPoint()\n");
      status->error = 1;
      statusstopbits |= CRITERION_error;
    }
  }

  return retvalue;
}

int asyncde::CDEIterator::GenerateUniformRandomPoint(ADEPoint &_point) {
  const std::vector<double> *lowerbound =
      problem->IntVariables()->InitialLowerLimits();
  const std::vector<double> *upperbound =
      problem->IntVariables()->InitialUpperLimits();

  _point.Reset();

  // check whether all lower and upper bounds are specified
  if (!lowerbound || !upperbound)
    //    return generate_MC_random_point(_point);
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
            "asyncde::CDEIterator::GenerateUniformRandomPoint():\n");
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

int asyncde::CDEIterator::AssignVector(const int vector_type,
                                       const std::vector<double> *&vector,
                                       const ADEPoint **_point) {
  int index = -1;
  vector = 0;
  const ADEPoint *newpoint = 0;
  int poolsize = adecfg->MinPopSize() - vetovector.size();
  int maxsize;

  switch (vector_type) {
  case ADE_VECTOR_BEST:
    index = 0;
    newpoint = base_population[pid_index_for_value_key[index]];
    break;
  case ADE_VECTOR_WORST:
    index = adecfg->MinPopSize() - 1;
    newpoint = base_population[pid_index_for_value_key[index]];
    break;
  case ADE_VECTOR_LINBEST:
    index = int(poolsize * sqrt(cfg->rnd->next(0.0, 1.0)));
    if (index > poolsize - 1)
      index = poolsize - 1;
    index = poolsize - 1 - index;
    break;
  case ADE_VECTOR_LINWORST:
    index = int(poolsize * sqrt(cfg->rnd->next(0.0, 1.0)));
    if (index > poolsize - 1)
      index = poolsize - 1;
    break;
  case ADE_VECTOR_PWORST:
    maxsize = int(adecfg->pworst * adecfg->MinPopSize() + 0.5);
    if (maxsize < 2)
      index = poolsize - 1;
    else
      index = (maxsize < poolsize)
                  ? poolsize - 1 - cfg->rnd->next_ulong(maxsize)
                  : cfg->rnd->next_ulong(poolsize);
    break;
  case ADE_VECTOR_CURRENTTOPBEST:
  case ADE_VECTOR_TARGETTOPBEST:
    maxsize = int(adecfg->pbest * adecfg->MinPopSize() + 0.5);
    if (maxsize < 2)
      index = 0;
    else
      index = (maxsize < poolsize) ? cfg->rnd->next_ulong(maxsize)
                                   : cfg->rnd->next_ulong(poolsize);
    break;
  case ADE_VECTOR_ARAND:
    if (archive) {
      index = cfg->rnd->next_ulong(poolsize + archive->Size());
      if (index >= poolsize)
        newpoint = (const ADEPoint *)archive->GetPoint(index - poolsize);
      break;
    }
    /* FALLTHRU */
  case ADE_VECTOR_RAND:
  case ADE_VECTOR_RANDTOBEST:
  case ADE_VECTOR_CURRENTTORAND:
  case ADE_VECTOR_TARGETTORAND:
  case ADE_VECTOR_CURRENTTOBEST:
  case ADE_VECTOR_TARGETTOBEST:
  default:
    index = cfg->rnd->next_ulong(poolsize);
  }

  if (!newpoint) {
    // refine index value
    for (int ind : vetovector)
      if (index >= ind)
        index++;

    if (index >= (int)adecfg->MinPopSize())
      fprintf(errors_stream, "Error in asyncde::CDEIterator::AssignVector()\n");

    newpoint = base_population[pid_index_for_value_key[index]];

    InsertIntoVetoVector(index);
  }

  vector = newpoint->Data()->X();
  if (_point)
    *_point = newpoint;

  return 0;
}

/*
// classical crossover O(nfreeparams)
int asyncde::CDEIterator::CrossoverMaskUniform(double CR,
                                               std::vector<char> &mask) const {
  int nmutantcoords = nfreeparams;

  if (nfreeparams < 1) {
    std::fill(mask.begin(), mask.end(), -1);
    return -2;
  }

  // enforce one coordinate from a mutant vector
  unsigned int R_internal = cfg->rnd->next_ulong(nfreeparams);

  for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
    if (cfg->rnd->next(0.0, 1.0) > CR && ivar != R_internal) {
      mask[ivar] = 0;
      nmutantcoords--;
    } else
      mask[ivar] = 1;

  return nmutantcoords;
}
*/

/*
// accelerated CrossoverMaskUniform (reduced number of next_ulong() calls)
int asyncde::CDEIterator::CrossoverMaskUniform(double CR,
                                               std::vector<char> &mask) const {
  if (nfreeparams < 1) {
    std::fill(mask.begin(), mask.end(), -1);
    return -2;
  }

  unsigned int nmutantcoords = cfg->rnd->next_binomialuint(nfreeparams - 1, CR);
  // enforce one coordinate from a mutant vector
  if (nmutantcoords < 1)
    nmutantcoords = 1;

  int value;
  unsigned int nfills;
  if (nmutantcoords < nfreeparams / 2) {
    value = 1;
    nfills = nmutantcoords;
    std::fill(mask.begin(), mask.end(), 0);
  } else {
    value = 0;
    nfills = nfreeparams - nmutantcoords;
    std::fill(mask.begin(), mask.end(), 1);
  }

  for (unsigned int irand = 0; irand < nfills; irand++) {
    unsigned int urand = cfg->rnd->next_ulong(nfreeparams - irand);
    unsigned int counter = 0;
    for (unsigned int j = 0; j < nfreeparams; j++)
      if (mask[j] != value) {
        if (counter == urand) {
          mask[j] = value;
          break;
        }
        counter++;
      }
  }

  return nmutantcoords;
}
*/

// accelerated CrossoverMaskUniform (requires std::vector<int> tmpindices)
int asyncde::CDEIterator::CrossoverMaskUniform(const double CR,
                                               std::vector<char> &mask) {
  if (nfreeparams < 1) {
    std::fill(mask.begin(), mask.end(), -1);
    return -2;
  }

  unsigned int nmutantcoords = cfg->rnd->next_binomialuint(nfreeparams - 1, CR);
  // enforce one coordinate from a mutant vector
  if (nmutantcoords < 1)
    nmutantcoords = 1;

  int value;
  unsigned int nfills;
  if (nmutantcoords < nfreeparams / 2) {
    value = 1;
    nfills = nmutantcoords;
    std::fill(mask.begin(), mask.end(), 0);
  } else {
    value = 0;
    nfills = nfreeparams - nmutantcoords;
    std::fill(mask.begin(), mask.end(), 1);
  }

  tmpindices.resize(nfreeparams);
  std::iota(tmpindices.begin(), tmpindices.end(), 0);

  for (unsigned int irand = 0; irand < nfills; irand++) {
    unsigned int urand = cfg->rnd->next_ulong(nfreeparams - irand);
    mask[tmpindices[urand]] = value;
    tmpindices[urand] = tmpindices[nfreeparams - irand - 1];
  }

  return nmutantcoords;
}

int asyncde::CDEIterator::CrossoverMask(const ADEPoint &target_point,
                                        ADEPoint &_point) {
  double CRlocal;

  switch (adecfg->CRupdatetype) {
  case ADE_CROSSOVER_UPDATE_jDE:
    CRlocal = target_point.ADEInfo()->CR;
    if (CRlocal < adecfg->CRmin || CRlocal > adecfg->CRmax ||
        cfg->rnd->next(0.0, 1.0) > adecfg->tauCR)
      CRlocal =
          (adecfg->CRmax - adecfg->CRmin > std::numeric_limits<double>::min())
              ? cfg->rnd->next(adecfg->CRmin, adecfg->CRmax)
              : adecfg->CRmin;
    break;
  case ADE_CROSSOVER_UPDATE_JADE:
    CRlocal = (adecfg->CRmax - adecfg->CRmin > 0.1)
                  ? cfg->rnd->BoxMuller(CRmu, adecfg->CRsigma, adecfg->CRmin,
                                        adecfg->CRmax)
                  : adecfg->CRmin;
    break;
  case ADE_CROSSOVER_UPDATE_CRCauchy:
  default:
    CRlocal = (adecfg->CRmax > adecfg->CRmin)
                  ? cfg->rnd->randCauchyTruncated(adecfg->CRmu, adecfg->CRsigma,
                                                  adecfg->CRmin, adecfg->CRmax)
                  : adecfg->CRmu;
  }

  _point.ADEInfoMutable()->CR = CRlocal;
  int nmutantcoords = CrossoverMaskUniform(_point.ADEInfoMutable()->CR,
                                           _point.ADEInfoMutable()->mask);

  if (nmutantcoords < 1) {
    statusstopbits |= CRITERION_error;
    return -2;
  }

#ifdef LC_DEBUG
  if (cfg->verbose >= 3) {
    fprintf(messages_stream, "F=%.5e  CR=%.5e\n", _point.ADEInfo()->F,
            _point.ADEInfoMutable()->CR);
    fprintf(messages_stream, "mutant vector:");
    for (double xcoord : *_point.Data()->X())
      fprintf(messages_stream, " %.3e", xcoord);
    fprintf(messages_stream, "\n");
  }
#endif // LC_DEBUG

  return nmutantcoords;
}

int asyncde::CDEIterator::CrossoverApplyMask(
    const std::vector<double> &target_vector,
    std::vector<double> &mutant_vector, const std::vector<char> &mask,
    int &different) const {
  int nmutantcoords = nfreeparams;

  for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
    if (0 == mask[ivar]) {
      mutant_vector[ivar] = target_vector[ivar];
      nmutantcoords--;
    } else if (!different) {
      if (0 != fuzzy_cmp(mutant_vector[ivar], target_vector[ivar],
                         problem->Tolerance()))
        different = 1;
    }

  return nmutantcoords;
}

int asyncde::CDEIterator::CrossoverIntADEPoint(const ADEPoint *target_point,
                                               ADEPoint &_point) {
  int different = (adecfg->mvector_unique == ADE_MVECTOR_DIFF_TARGET) ? 0 : 1;
  int nmutantcoords = CrossoverApplyMask(
      *target_point->Data()->X(), *_point.DataMutable()->Xmutable(),
      _point.ADEInfoMutable()->mask, different);

  _point.ADEInfoMutable()->CR =
      (nmutantcoords > 0) ? nmutantcoords / ((double)nfreeparams) : 0.0;

  return (different) ? 0 : -8;
}

int asyncde::CDEIterator::InitVetoVector() {
  int ivalue = -1;
  int ivaluetrial = -1;

  vetovector.clear();

  // remove fixed vectors
  switch (adecfg->targetvector) {
  case ADE_VECTOR_BEST:
    ivaluetrial = 0;
    break;
  case ADE_VECTOR_WORST:
    ivaluetrial = adecfg->MinPopSize() - 1;
    break;
  }
  if (ivaluetrial >= 0)
    InsertIntoVetoVector(ivaluetrial);

  switch (adecfg->basevector) {
  case ADE_VECTOR_BEST:
    ivalue = 0;
    break;
  case ADE_VECTOR_WORST:
    ivalue = adecfg->MinPopSize() - 1;
    break;
  }
  if (ivalue >= 0 && ivalue != ivaluetrial)
    InsertIntoVetoVector(ivalue);

#ifdef LC_DEBUG
  if (cfg->verbose >= 3) {
    fprintf(messages_stream, "InitVetoVector():\n");
    for (unsigned int ip = 0; ip < vetovector.size(); ip++)
      fprintf(messages_stream, "vetovector[%i] pid=%li\n", ip,
              pid_key[pid_index_for_value_key[vetovector[ip]]]);
  }
#endif // LC_DEBUG

  return adecfg->MinPopSize() - vetovector.size();
}

int asyncde::CDEIterator::InsertIntoVetoVector(int valueindex) {
  int ipos = 0;
  int vetovectorsize = vetovector.size();

  for (ipos = 0; ipos < vetovectorsize; ipos++)
    if (valueindex <= vetovector[ipos]) {
      if (valueindex == vetovector[ipos])
        return 0;

      vetovector.push_back(valueindex);
      if (ipos < (int)vetovector.size() + 1) {
        for (int ind = vetovector.size() - 1; ind > ipos; ind--)
          vetovector[ind] = vetovector[ind - 1];
        vetovector[ipos] = valueindex;
      }
      return 0;
    }

  vetovector.push_back(valueindex);

  return 0;
}

int asyncde::CDEIterator::FillInMutantIntADEPoint(const ADEPoint *target_point,
                                                  ADEPoint &_point) {
  const std::vector<double> *target_vector = target_point->Data()->X();
  const std::vector<double> *base_vector = nullptr;
  const std::vector<double> *diff_vector1 = nullptr;
  const std::vector<double> *diff_vector2 = nullptr;
  const std::vector<double> *best_vector =
      base_population[pid_index_for_value_key[0]]->Data()->X();

  std::vector<double> *x = _point.DataMutable()->Xmutable();

  const ADEPoint *diff_point1 = nullptr;
  const ADEPoint *diff_point2 = nullptr;

  if (adecfg->basevector != ADE_VECTOR_CURRENTTOBEST &&
      adecfg->basevector != ADE_VECTOR_TARGETTOBEST)
    AssignVector(adecfg->basevector, base_vector, &diff_point1);

#ifdef LC_DEBUG
  if (cfg->verbose >= 3) {
    fprintf(messages_stream,
            "target vector (%li):", target_point->ADEInfo()->ParentId());
    for (double xcoord : *target_vector)
      fprintf(messages_stream, " %.3e", xcoord);
    if (base_vector) {
      fprintf(messages_stream,
              "\nbase vector (%li):", diff_point1->ADEInfo()->ParentId());
      for (double xcoord : *base_vector)
        fprintf(messages_stream, " %.3e", xcoord);
      fprintf(messages_stream, "\n");
    }
  }
#endif // LC_DEBUG

  switch (adecfg->Fscaletype) {
  case ADE_FSCALE_FCauchy:
    _point.ADEInfoMutable()->F =
        (adecfg->Fmax > adecfg->Fmin)
            ? cfg->rnd->randCauchyTruncated(adecfg->Fmu, adecfg->Fsigma,
                                            adecfg->Fmin, adecfg->Fmax)
            : adecfg->Fmu;
    break;
  case ADE_FSCALE_jDE:
    // jDE adaptation
    if (target_point->ADEInfo()->F >= adecfg->Fmin &&
        target_point->ADEInfo()->F <= adecfg->Fmax &&
        cfg->rnd->next(0.0, 1.0) > adecfg->tauF)
      _point.ADEInfoMutable()->F = target_point->ADEInfo()->F;
    else
      _point.ADEInfoMutable()->F =
          (adecfg->Fmax - adecfg->Fmin > std::numeric_limits<double>::min())
              ? cfg->rnd->next(adecfg->Fmin, adecfg->Fmax)
              : adecfg->Fmin;
    break;
  case ADE_FSCALE_JADE:
  default:
    _point.ADEInfoMutable()->F =
        (adecfg->Fmax > adecfg->Fmin)
            ? cfg->rnd->randCauchyTruncated(FmuJADE, adecfg->Fsigma,
                                            adecfg->Fmin, adecfg->Fmax)
            : adecfg->Fmu;
    /*
    fprintf(stderr, "%.3e = Cauchy(%.3e, %.3e) in [%.3e, %.3e]\n",
            _point.ADEInfo()->F, FmuJADE, adecfg->Fsigma, adecfg->Fmin,
            adecfg->Fmax);
    */
  }

  // specific strategies (rand-to-best, current-to-rand) modify the base
  // vector
  switch (adecfg->basevector) {
  case ADE_VECTOR_RANDTOBEST:
    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
      (*x)[ivar] =
          (*base_vector)[ivar] +
          _point.ADEInfo()->F * ((*best_vector)[ivar] - (*base_vector)[ivar]);
    break;
  case ADE_VECTOR_CURRENTTORAND:
  case ADE_VECTOR_TARGETTORAND:
  case ADE_VECTOR_CURRENTTOPBEST:
  case ADE_VECTOR_TARGETTOPBEST:
    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
      (*x)[ivar] =
          (*target_vector)[ivar] +
          _point.ADEInfo()->F * ((*base_vector)[ivar] - (*target_vector)[ivar]);
    break;
  case ADE_VECTOR_CURRENTTOBEST:
  case ADE_VECTOR_TARGETTOBEST:
    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++)
      (*x)[ivar] =
          (*target_vector)[ivar] +
          _point.ADEInfo()->F * ((*best_vector)[ivar] - (*target_vector)[ivar]);
    break;
  default:
    *x = *base_vector;
    break;
  }

  for (unsigned int idiff = 0; idiff < adecfg->ndifferences; idiff++) {
    AssignVector(ADE_VECTOR_RAND, diff_vector1, &diff_point1);
#ifdef LC_DEBUG
    if (cfg->verbose >= 3) {
      fprintf(messages_stream, "r1(%2i) vector (%li):", idiff,
              diff_point1->ADEInfo()->ParentId());
      for (double xcoord : *diff_vector1)
        fprintf(messages_stream, " %.3e", xcoord);
      fprintf(messages_stream, "\n");
    }
#endif // LC_DEBUG
    if (adecfg->archive)
      AssignVector(ADE_VECTOR_ARAND, diff_vector2, &diff_point2);
    else
      AssignVector(ADE_VECTOR_RAND, diff_vector2, &diff_point2);
#ifdef LC_DEBUG
    if (cfg->verbose >= 3) {
      fprintf(messages_stream, "r2(%2i) vector (%li):", idiff,
              diff_point2->ADEInfo()->ParentId());
      for (double xcoord : *diff_vector2)
        fprintf(messages_stream, " %.3e", xcoord);
      fprintf(messages_stream, "\n");
    }
#endif // LC_DEBUG

    const double F = _point.ADEInfo()->F;

    int similarvar = -1;
    for (unsigned int ivar = 0; ivar < nfreeparams; ivar++) {
      (*x)[ivar] += F * ((*diff_vector2)[ivar] - (*diff_vector1)[ivar]);
      if (similarvar == -1 && _point.ADEInfo()->mask[ivar] &&
          0 == fuzzy_cmp((*diff_vector1)[ivar], (*diff_vector2)[ivar],
                         cfg->criteriarestart.xepsilon))
        similarvar = ivar;
    }

    if (similarvar != -1)
      status->Incr_nxdiff();
  }

  _point.InfoMutable()->SetStatus(POINTSTATUS_INITIALIZED);

  return 0;
}

int asyncde::CDEIterator::FillInTrialIntADEPoint(ADEPoint &_point) {
  _point.Reset();

  if (!base_population[adecfg->MinPopSize() - 1] ||
      base_population[adecfg->MinPopSize() - 1]->Info()->Id() ==
          PointInfo::NOT_INITIALIZED_POINT)
    return -1;

  const ADEPoint *target_point = nullptr;
  const std::vector<double> *target_vector = nullptr;

  const std::vector<double> *xfeas;

  status->nprojfeas = 0;
  unsigned int nmax_attempts = std::numeric_limits<unsigned int>::max();
  if ((cfg->criteriastop.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriastop.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriastop.nprojfeas;
  if ((cfg->criteriarestart.statusbits & CRITERION_nprojfeas) != 0)
    if (cfg->criteriarestart.nprojfeas < nmax_attempts)
      nmax_attempts = cfg->criteriarestart.nprojfeas;
  long int tag;

#ifdef LC_DEBUG
  if (cfg->verbose >= 3) {
    fprintf(messages_stream, "asyncde::CDEIterator::generate_ade_point():\n");
    if (cfg->verbose > 1) {
      fprintf(messages_stream, "initial population:\n");
      for (unsigned int ip = 0; ip < base_population_actual_size; ip++)
        if (base_population[ip])
          base_population[ip]->Print(messages_stream);
      fprintf(messages_stream, "----------\n");
    }
  }
#endif // LC_DEBUG

  do {
    InitVetoVector();

    AssignVector(adecfg->targetvector, target_vector, &target_point);
    _point.ADEInfoMutable()->SetParentId(target_point->ADEInfo()->ParentId());

    CrossoverMask(*target_point, _point);

    FillInMutantIntADEPoint(target_point, _point);

    tag = PointInfo::NOT_INITIALIZED_POINT;
    if (CrossoverIntADEPoint(target_point, _point) < 0)
      tag = 0;

    // project the trial vector into the feasible region
    if (1 != problem->IsXIntFeasible(*_point.Data()->X(), xexttmpvector)) {
      xfeas = (1 == problem->IsXIntFeasible(*target_vector, xexttmpvector))
                  ? target_vector
                  : BestIntPointCurrentPopulation()->Data()->X();
#ifdef LC_DEBUG
      if (cfg->verbose >= 3) {
        fprintf(messages_stream, "trial vector is infeasible\nproject point to "
                                 "feasibility, select xfeas:\n");
        for (double xcoord : *xfeas)
          fprintf(messages_stream, " %.3e", xcoord);
        fprintf(messages_stream, "\n");
      }
#endif // LC_DEBUG

      if (0 != ProjectPoint2Feasibility(*_point.DataMutable()->Xmutable(),
                                        *xfeas, xexttmpvector))
        tag = 0;
      //#ifdef LC_DEBUG
      if (cfg->verbose >= 3) {
        fprintf(messages_stream, "CDEIterator::FillInTrialIntADEPoint(): "
                                 "failed to generate feasible trial vector:\n");
        _point.Print(messages_stream);
      }
      //#endif // LC_DEBUG
    }
#ifdef LC_DEBUG
    else if (cfg->verbose >= 3)
      fprintf(messages_stream, "trial vector is feasible\n");
#endif // LC_DEBUG

    if (adecfg->mvector_unique == ADE_MVECTOR_UNIQUE_INPOP)
      IsContainsX(*_point.Data()->X(), tag);
    status->Incr_nProjFeas();
  } while ((tag != PointInfo::NOT_INITIALIZED_POINT) &&
           status->nprojfeas < nmax_attempts);

  if (tag == PointInfo::NOT_INITIALIZED_POINT) {
    _point.ADEInfoMutable()->SwitchToNewId();
    _point.ADEInfoMutable()->SetStatus(POINTSTATUS_INITIALIZED);

#ifdef LC_DEBUG
    if (cfg->verbose >= 3) {
      fprintf(
          messages_stream,
          "asyncde::CDEIterator::generate_ade_point() successfully finished\n");

      fprintf(messages_stream, "new target point:\n");
      _point.Print(messages_stream);
    }
#endif // LC_DEBUG

  } else {
    _point.ADEInfoMutable()->SetStatus(POINTSTATUS_INFEASIBLE);
    //#ifdef LC_DEBUG
    fprintf(errors_stream,
            "asyncde::CDEIterator::generate_ade_point() failed\n");
    fprintf(errors_stream, "failed target point:\n");
    _point.Print(errors_stream);
    //#endif // LC_DEBUG
  }

  return 0;
}

int asyncde::CDEIterator::UpdateStatSums(const ADEPoint *new_point,
                                         const ADEPoint *old_point) {
  int retvalue = 0;
  double Flocal, CRlocal;

  statsumsupdates_counter++;
  if (statsumsupdates_counter > statsumsupdates_counter_max) {
    statsumsupdates_counter = 0;

    Fsum = 0.0;
    F2sum = 0.0;
    CRn = 0;
    CRsum = 0.0;
    CR2sum = 0.0;

    for (unsigned int ip = 0; ip < base_population_actual_size; ip++) {
      const ADEPoint *pointptr = base_population[ip];
      Flocal = pointptr->ADEInfo()->F;
      if (Flocal >= adecfg->Fmin) {
        Fsum += Flocal;
        F2sum += Flocal * Flocal;
      }

      CRlocal = pointptr->ADEInfo()->CR;
      if (CRlocal >= adecfg->CRmin) {
        CRn++;
        CRsum += CRlocal;
        CR2sum += CRlocal * CRlocal;
      }
    }
  }

  if (old_point) {
    Flocal = old_point->ADEInfo()->F;
    if (Flocal >= adecfg->Fmin) {
      Fsum -= Flocal;
      F2sum -= Flocal * Flocal;
    }
    CRlocal = old_point->ADEInfo()->CR;
    if (CRlocal >= adecfg->CRmin) {
      CRn--;
      CRsum -= CRlocal;
      CR2sum -= CRlocal * CRlocal;
    }
  }

  if (new_point) {
    Flocal = new_point->ADEInfo()->F;
    if (Flocal >= adecfg->Fmin) {
      Fsum += Flocal;
      F2sum += Flocal * Flocal;

      FmuJADE += adecfg->Fc * (F2sum / Fsum - FmuJADE);
      if (FmuJADE < adecfg->Fmin)
        FmuJADE = adecfg->Fmin;
      else if (FmuJADE > adecfg->Fmax)
        FmuJADE = adecfg->Fmax;
    }
    CRlocal = new_point->ADEInfo()->CR;
    if (CRlocal >= adecfg->CRmin) {
      CRn++;
      CRsum += CRlocal;
      CR2sum += CRlocal * CRlocal;

      CRmu += adecfg->CRc * (CRsum / CRn - CRmu);
      if (CRmu < adecfg->CRmin)
        CRmu = adecfg->CRmin;
      else if (CRmu > adecfg->CRmax)
        CRmu = adecfg->CRmax;
    }
  }

  return retvalue;
}

int asyncde::CDEIterator::Print(FILE *stream) const {
  int retvalue = adecfg->Print(stream);

  // print status (progress/stopped)

  // print best so far value

  return retvalue;
}
