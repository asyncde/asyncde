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

#include <chrono>
#include <limits>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// #define _GPERF_PROFILER

#ifdef _GPERF_PROFILER
#include <gperftools/profiler.h>
#endif // _GPERF_PROFILER

#include "asyncde/ADEConfig.h"
#include "asyncde/ADEIterator.h"
#include "asyncde/ADEPoint.h"
#include "asyncde/ADEPointInfo.h"
#include "asyncde/AsyncIterator.h"
#include "asyncde/Function.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/ProblemMapping.h"
#include "asyncde/Rnd.h"
#include "asyncde/Variables.h"

// test functions
#include "asyncde/composite_rosenbrock.h"
#include "asyncde/rosenbrock.h"
#include "asyncde/sqr.h"

asyncde::Variables *init_variables(unsigned int nunknowns, asyncde::Rnd &rnd) {
  if (nunknowns < 1)
    return 0;

  const double xmax = 10.0;
  const double xmax_initial = xmax;

  asyncde::Variables *vars = new asyncde::Variables(nunknowns);
  for (unsigned int ip = 0; ip < nunknowns; ip++) {
    vars->SetName(ip, "");
    vars->SetCoordinate(ip, rnd.next(-xmax_initial, xmax_initial));
    //    vars->SetLowerLimit(ip, -xmax);
    //    vars->SetUpperLimit(ip, xmax);
    vars->SetInitialRange(ip, -xmax_initial, xmax_initial);
    vars->ReleaseVariable(ip);
  }

#ifdef TEST_FIXED_VARS
  int index = vars->AddVariable("fixed", 0.0);
  vars->FixVariable(index);
#endif

  return vars;
}

int ade_example(asyncde::Functor1D &fitnessfunctor, unsigned int dim,
                double &vtr, double &maxfunevals, unsigned int &pop_size) {
  int retvalue = 0;
  double iter;

  time_t t0 = time(0);
  asyncde::Rnd rnd(t0);

  asyncde::PointInfo::ResetIdCounter();
  asyncde::Variables *vars = init_variables(dim, rnd);
  asyncde::ProblemMapping problem;
  std::vector<double> f(1);
  problem.SetYsize(f.size());
  problem.SetExternalVariables(vars);
  problem.Print(stdout);

  int nparents = pop_size;
  asyncde::ADEConfig adesettings(&rnd, nparents);
  adesettings.SetStrategy("DE/rand/rand/1/acm");

  adesettings.criteriastop.SetCriterion_vtr(vtr);
  adesettings.criteriarestart.Reset();
  adesettings.criteriarestart.SetCriterion_ydelta(1.e-11); // 1.e-3 * 1.e-8
  adesettings.criteriarestart.SetCriterion_nxdiff(5);
  adesettings.criteriarestart.SetCriterion_PopSizeMax(1280);
  adesettings.criteriarestart.SetCriterion_nProjFeas(100);
  //  adesettings.criteriarecover.nrecover = 1;

  adesettings.Print(stdout);

  asyncde::AsyncIterator *iterator = adesettings.NewIterator(problem);

  asyncde::Point *tmp_point = iterator->NewExtPoint();
  const asyncde::Point *bestpoint = nullptr;

  for (iter = 0.; iter < maxfunevals; iter++) {
    iterator->FillInTrialExtPoint(*tmp_point);

    if ((retvalue |= iterator->StatusStopBits()))
      break;

    fitnessfunctor(*tmp_point->Data()->X(), f);
    tmp_point->SetY(&f);
    iterator->AddExtPoint(*tmp_point);

    if ((retvalue |= iterator->StatusStopBits()))
      break;
  }

  bestpoint = iterator->BestIntPoint();
  vtr = (bestpoint) ? (*bestpoint->Data()->Y())[0]
                    : std::numeric_limits<double>::max();
  maxfunevals = iterator->NFE();

  fprintf(stdout,
          "ade_example():\nretvalue = %i  nevals = %li  bestvalue = %.15e\n",
          retvalue, iterator->NFE(), vtr);

  if (bestpoint) {
    problem.ConvertInt2Ext(*bestpoint->Data()->X(),
                           *tmp_point->DataMutable()->Xmutable());
    for (std::vector<double>::const_iterator it =
             tmp_point->Data()->X()->begin();
         it != tmp_point->Data()->X()->end(); it++) {
      if (it != tmp_point->Data()->X()->begin())
        fprintf(stdout, "  ");
      fprintf(stdout, "%.15e", *it);
    }
    fprintf(stdout, "\n");
    for (std::vector<double>::const_iterator it =
             bestpoint->Data()->Y()->begin();
         it != bestpoint->Data()->Y()->end(); it++) {
      if (it != bestpoint->Data()->Y()->begin())
        fprintf(stdout, "  ");
      fprintf(stdout, "%.15e", *it);
    }
    fprintf(stdout, "\n");

    asyncde::ADEIterator *ADEiterator =
        dynamic_cast<asyncde::ADEIterator *>(iterator);
    if (ADEiterator)
      ADEiterator->PrintCorrMatrix();
  }

  delete iterator;
  delete vars;
  delete tmp_point;

  return retvalue;
}

int stat_test() {
  const unsigned int ndimensions = 7;
  const unsigned int Ntrials = 100;
  int ndim[ndimensions] = {2, 3, 5, 10, 20, 50, 100};
  const unsigned int initial_pop_size = 10;
  const double ftarget = 120.0;
  const double vtr = ftarget + 1.e-8;
  const double maxfunevals = 1.e6;

  unsigned int pop_size;
  double fmin, nevals;
  asyncde::composite_rosenbrock crosenbrock(ftarget);

  for (unsigned int idim = 0; idim < 6; idim++) {
    double stat_nevals = 0;
    double stat_nevals2 = 0;
    int nsuccess = 0;
    int retvalue;
    for (unsigned int itrial = 0; itrial < Ntrials; itrial++) {
      fmin = vtr;
      nevals = maxfunevals * ndim[idim];
      pop_size = initial_pop_size;
      retvalue = ade_example(crosenbrock, ndim[idim], fmin, nevals, pop_size);
      if (retvalue == 1) {
        nsuccess++;
        stat_nevals += nevals;
        stat_nevals2 += nevals * nevals;
      }
    }

    fprintf(stdout,
            "Stat. properties for dim=%i: Psucc=%.3e  MeanN=%.3e +- %.3e\n",
            ndim[idim], nsuccess / (1.0 * Ntrials),
            (nsuccess > 0) ? stat_nevals / nsuccess : 0.0,
            (nsuccess > 1)
                ? sqrt((nsuccess * stat_nevals2 - stat_nevals * stat_nevals) /
                       (nsuccess * (nsuccess - 1.0)))
                : 0.0);
  }

  return 0;
}

int main() {
  std::chrono::steady_clock timer;
  std::chrono::time_point<std::chrono::steady_clock> t0 = timer.now();

  unsigned int nvars = 100;
  unsigned int pop_size = 10;
  const double ftarget = 120.0;
  double vtr = ftarget + 1.e-8;
  double maxfunevals = 1.e7;
  //  asyncde::sqr fitnessfunctor(ftarget);
  asyncde::rosenbrock fitnessfunctor(ftarget);
  //  asyncde::composite_rosenbrock fitnessfunctor(ftarget);

#ifdef _GPERF_PROFILER
  ProfilerStart("ade_example.prof");
#endif
  int retvalue = ade_example(fitnessfunctor, nvars, vtr, maxfunevals, pop_size);
#ifdef _GPERF_PROFILER
  ProfilerStop();
#endif

  std::chrono::nanoseconds dt = timer.now() - t0;
  double tduration = dt.count(); // in nsec

  fprintf(stdout,
          "ade_example():\nretvalue = %i  nevals = %g  bestvalue = %.15e\n",
          retvalue, maxfunevals, vtr);
  fprintf(stdout, "final pop_size=%u\n", pop_size);
  fprintf(stdout, "timing: %.3e sec   %.3e ns / 1 feval\n", 1.e-9 * tduration,
          tduration / maxfunevals);

  /*
    stat_test();
  */
  return 0;
}
