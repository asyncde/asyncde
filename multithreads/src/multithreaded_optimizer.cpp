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
#include <mutex>
#include <thread>

#include "asyncde/AsyncIterator.h"
#include "asyncde/Function.h"
#include "asyncde/IteratorConfig.h"
#include "asyncde/Point.h"
#include "asyncde/Problem.h"
#include "asyncde/StopCriteria.h"

#include "asyncde/multithreaded_optimizer.h"

int asyncde::multithreaded_optimization_cycle(
    const Problem &problem, const IteratorConfig &cfg, double &maxfunevals,
    double &bestvalue, std::vector<double> *bestX,
    std::vector<double> * /*parerrlow*/, std::vector<double> * /*parerrup*/) {
  int retvalue = 0;
  bestvalue = std::numeric_limits<double>::max();

  const Functor1D *fitnessfunctor = problem.ObjectiveFunctor1D();
  if (!fitnessfunctor)
    return -1;

  AsyncIterator *iterator = cfg.NewIterator(problem);
  retvalue = iterator->StatusStopBits();

  const size_t nthreads = std::thread::hardware_concurrency();
  fprintf(stderr, "nthreads = %li\n", nthreads);

  std::vector<std::thread> threadpool;
  threadpool.reserve(nthreads);
  std::mutex mutx;

  const Point *bestpoint = nullptr;

  for (size_t ithr = 0; ithr < nthreads; ++ithr)
    threadpool.emplace_back([&]() {
      Point *tmp_point = nullptr;
      std::vector<double> *Y = nullptr;
      {
        std::unique_lock<std::mutex> lck(mutx);
        tmp_point = iterator->NewExtPoint();
        Y = new std::vector<double>(problem.Y()->size());
        (*Y)[0] = std::numeric_limits<double>::max();
      }

      while (retvalue == 0) {
        {
          std::unique_lock<std::mutex> lck(mutx);
          iterator->FillInTrialExtPoint(*tmp_point);
        }

        if (tmp_point->Info()->Id() != POINTSTATUS_UNDEFINED) {
          (*fitnessfunctor)(*tmp_point->Data()->X(), *Y);
          tmp_point->SetY(Y);
        }

        if (retvalue == 0 && tmp_point->Info()->Id() != POINTSTATUS_UNDEFINED) {
          std::unique_lock<std::mutex> lck(mutx);
          iterator->AddExtPoint(*tmp_point);
          retvalue |= iterator->StatusStopBits();
        }
      }

      {
        std::unique_lock<std::mutex> lck(mutx);
        delete tmp_point;
        delete Y;
      }
    });

  for (size_t ithr = 0; ithr < threadpool.size(); ++ithr)
    threadpool[ithr].join();

  retvalue |= iterator->StatusStopBits();
  maxfunevals = iterator->NFE();

  if ((bestpoint = iterator->BestIntPoint())) {
    bestvalue = (*bestpoint->Data()->Y())[0];
    if (bestX)
      problem.ConvertInt2Ext(*bestpoint->Data()->X(), *bestX);
    /*
    if (parerrlow)
      *parerrlow = iterator->parerrlow;
    if (parerrup)
      *parerrup = iterator->parerrup;
      */
  }

  delete iterator;

  return retvalue;
}
