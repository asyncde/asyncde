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

#include <limits>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "asyncde/ADEConfig.h"
#include "asyncde/Function.h"
#include "asyncde/PointInfo.h"
#include "asyncde/ProblemMapping.h"
#include "asyncde/Rnd.h"
#include "asyncde/Variables.h"

#include "asyncde/multithreaded_optimizer.h"

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

int threads_example(asyncde::Functor1D &fitnessfunctor, unsigned int dim,
                    double &vtr, double &maxfunevals, unsigned int &pop_size) {
  int retvalue = 0;

  time_t t0 = time(0);
  asyncde::Rnd rnd(t0);

  asyncde::PointInfo::ResetIdCounter();
  asyncde::Variables *vars = init_variables(dim, rnd);
  asyncde::ProblemMapping problem;
  std::vector<double> f(1);
  problem.SetYsize(f.size());
  problem.SetExternalVariables(vars);
  problem.Print(stdout);
  problem.SetObjectiveFunctor1D(&fitnessfunctor);
  std::vector<double> bestX(dim);
  std::vector<double> parerrlow(dim);
  std::vector<double> parerrup(dim);

  int nparents = pop_size;
  asyncde::ADEConfig adesettings(&rnd, nparents);
  adesettings.SetStrategy("DE/rand/rand/1/acm");

  adesettings.criteriastop.SetCriterion_vtr(vtr);
  adesettings.criteriastop.SetCriterion_nFE((long unsigned int)maxfunevals);
  adesettings.criteriarestart.Reset();
  adesettings.criteriarestart.SetCriterion_ydelta(1.e-11); // 1.e-3 * 1.e-8
  adesettings.criteriarestart.SetCriterion_nxdiff(5);
  adesettings.criteriarestart.SetCriterion_PopSizeMax(800);

  adesettings.Print(stdout);

  retvalue = asyncde::multithreaded_optimization_cycle(
      problem, adesettings, maxfunevals, vtr, &bestX, &parerrlow, &parerrup);

  delete vars;

  return retvalue;
}

int main() {
  unsigned int nvars = 100;
  unsigned int pop_size = 10;
  const double ftarget = -120.;
  double vtr = ftarget + 1.e-8;
  double maxfunevals = 1.e7;
  asyncde::sqr fitnessfunctor(ftarget);
  //  asyncde::rosenbrock fitnessfunctor(ftarget);
  //  asyncde::composite_rosenbrock fitnessfunctor(ftarget);

  int retvalue =
      threads_example(fitnessfunctor, nvars, vtr, maxfunevals, pop_size);

  fprintf(stdout,
          "threads_example():\nretvalue = %i  nevals = %g  bestvalue = %.15e\n",
          retvalue, maxfunevals, vtr);
  fprintf(stdout, "final pop_size=%u\n", pop_size);

  return EXIT_SUCCESS;
}
