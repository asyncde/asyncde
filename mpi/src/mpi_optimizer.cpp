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
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#include "asyncde/mpi_optimizer.h"

#include "asyncde/AsyncIterator.h"
#include "asyncde/IteratorConfig.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Problem.h"
#include "asyncde/Variables.h"

#define OBSOLETEVOIDPTR
// #define OBSOLETEVOIDPTR (void *) // for old MPI_Send

// #define LC_DEBUG

const time_t time_report = 60;

#ifdef LC_DEBUG
const int CHAR_BUF_SIZE = 1024;
#endif // LC_DEBUG

int asyncde::mpi_worker_cycle(const unsigned int Xsize,
                              const unsigned int Ysize,
                              const asyncde::Functor1D &fitnessfunctor) {
  int retvalue = 0;
  int rank, mpi_count;
  std::vector<double> x(Xsize);
  std::vector<double> y(Ysize);
#ifdef LC_DEBUG
  char bufstr[CHAR_BUF_SIZE];
  int len, dlen;
#endif // LC_DEBUG

  MPI_Status Stat;

  if (Xsize < 1 || Ysize < 1)
    return -1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef LC_DEBUG
  fprintf(stderr, "mpi_worker_cycle(): rank = %d; xsize=%i ysize=%i\n", rank,
          Xsize, Ysize);
#endif // LC_DEBUG

  for (;;) {
    retvalue |= MPI_Recv(OBSOLETEVOIDPTR x.data(), x.size(), MPI_DOUBLE, 0,
                         MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);

    if (Stat.MPI_TAG == asyncde::MPI_TAG_STOP) {
#ifdef LC_DEBUG
      fprintf(stderr, "%d <- %d: MPI_Recv() STOP\n", rank, 0);
#endif // LC_DEBUG
      break;
    } else {
      retvalue |= MPI_Get_count(&Stat, MPI_DOUBLE, &mpi_count);

#ifdef LC_DEBUG
      len = snprintf(bufstr, sizeof(bufstr), "%d <- %d: MPI_Recv() size = %d",
                     rank, 0, mpi_count);
      for (double xcoord : x) {
        dlen = snprintf(&bufstr[len], sizeof(bufstr) - len, " %.3e", xcoord);
        len += dlen;
      }
      fprintf(stderr, "%s\n", bufstr);
#endif // LC_DEBUG

      if ((int)x.size() != mpi_count) {
        fprintf(
            stderr,
            "ERROR: %d <- %d: MPI_Recv() received %d doubles instead of %li \n",
            rank, 0, mpi_count, x.size());
        break;
      }

      fitnessfunctor(x, y);

      retvalue |= MPI_Send(y.data(), y.size(), MPI_DOUBLE, 0,
                           asyncde::MPI_TAG_DOUBLE_ARRAY, MPI_COMM_WORLD);

#ifdef LC_DEBUG
      fprintf(stderr, "%d -> %d: MPI_Send() y = %.5e\n", rank, 0, y[0]);
#endif // LC_DEBUG
    }
  }

  return retvalue;
}

int asyncde::mpi_master_cycle(const Problem &problem, const IteratorConfig &cfg,
                              long int &nevals, double &bestvalue,
                              std::vector<double> *bestX,
                              std::vector<double> * /*parerrlow*/,
                              std::vector<double> * /*parerrup*/) {
  int retvalue = 0;

  nevals = -1;
  bestvalue = std::numeric_limits<double>::max();
  const Point *bestpoint = nullptr;

  int numtasks, rank, source, rc;
  std::vector<double> y(problem.Y()->size());
#ifdef LC_DEBUG
  char bufstr[CHAR_BUF_SIZE];
  int len, dlen;
  fprintf(stderr, "mpi_master_cycle()\n");
#endif // LC_DEBUG

  time_t t0, t1;

  MPI_Status Stat;

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<int> rank_active(numtasks);
  int active;

  if (numtasks < 2) {
    fprintf(stderr, "number of processes should be >= 2\n");
    return -2;
  }
  if (rank != 0) {
    fprintf(stderr, "Internal error: not a master (rank=%i)\n", rank);
    return -4;
  }

  AsyncIterator *iterator = cfg.NewIterator(problem);
  asyncde::Point *tmp_point = iterator->NewExtPoint();
  retvalue = iterator->StatusStopBits();

  std::vector<Point *> buf_point;
  buf_point.resize(numtasks, nullptr);
  for (int dest = 1; dest < numtasks; dest++)
    buf_point[dest] = iterator->NewExtPoint();

  // initialize
  for (int dest = 1; dest < numtasks; dest++) {
    iterator->FillInTrialExtPoint(*buf_point[dest]);
    rank_active[dest] = 1;
  }

  for (int dest_worker = 1; dest_worker < numtasks; dest_worker++) {
    problem.ConvertInt2Ext(*buf_point[dest_worker]->Data()->X(),
                           *tmp_point->DataMutable()->Xmutable());
    rc = MPI_Send(OBSOLETEVOIDPTR tmp_point->Data()->X()->data(),
                  tmp_point->Data()->X()->size(), MPI_DOUBLE, dest_worker,
                  asyncde::MPI_TAG_DOUBLE_ARRAY, MPI_COMM_WORLD);

#ifdef LC_DEBUG
    len = snprintf(bufstr, sizeof(bufstr),
                   "%d -> %d: MPI_Send() of point: ", rank, dest_worker);
    for (double xcoord : *buf_point[dest_worker]->Data()->X()) {
      dlen = snprintf(&bufstr[len], sizeof(bufstr) - len, " %.3e", xcoord);
      len += dlen;
    }
    fprintf(stderr, "%s\n", bufstr);
#endif // LC_DEBUG
  }

#ifdef LC_DEBUG
  fprintf(stderr,
          "mpi_master_cycle(): initialization: after MPI_Send to workers\n");
  fprintf(stderr, "mpi_master_cycle(): wait for MPI_Recv from any worker\n");
#endif // LC_DEBUG

  t0 = time(0);

  // main cycle
  for (;;) {
    rc = MPI_Recv(y.data(), y.size(), MPI_DOUBLE, MPI_ANY_SOURCE,
                  asyncde::MPI_TAG_DOUBLE_ARRAY, MPI_COMM_WORLD, &Stat);
    source = Stat.MPI_SOURCE;

#ifdef LC_DEBUG
    fprintf(stderr, "mpi_master_cycle(): got MPI_Recv from worker %i\n",
            source);
#endif // LC_DEBUG

    if (source > 0 && source < numtasks) {
#ifdef LC_DEBUG
      fprintf(stderr, "%d <- %d MPI_Recv() y = %.5e\n", rank, source, y[0]);
#endif // LC_DEBUG

      buf_point[source]->SetY(&y);
      iterator->AddExtPoint(*buf_point[source]);
      t1 = time(0);

      if (t1 - t0 >= time_report) {
        t0 = t1;
        if ((bestpoint = iterator->BestIntPoint())) {
          bestvalue = (*bestpoint->Data()->Y())[0];
          problem.ConvertInt2Ext(*bestpoint->Data()->X(),
                                 *tmp_point->DataMutable()->Xmutable());
          printf("NFE = %li   Best value = %.15e\n", iterator->NFE(),
                 bestvalue);
          for (std::vector<double>::const_iterator it =
                   tmp_point->Data()->X()->begin();
               it != tmp_point->Data()->X()->end(); it++) {
            if (it != tmp_point->Data()->X()->begin())
              printf("  ");
            printf("%.15e", *it);
          }
          printf("\n");
        }
      }
    } else {
      fprintf(stderr, "Error: unknown source %d\n", source);
      break;
    }

    retvalue |= iterator->StatusStopBits();
    if (retvalue == 0) {
      retvalue = iterator->FillInTrialExtPoint(*buf_point[source]);
      if (retvalue)
        fprintf(stderr, "Failed to generate trial point\n");
    }

    if (retvalue == 0) {
      // send next trial point
      problem.ConvertInt2Ext(*buf_point[source]->Data()->X(),
                             *tmp_point->DataMutable()->Xmutable());
#ifdef LC_DEBUG
      fprintf(stderr, "mpi_master_cycle(): MPI_Send next point: 0 -> %i\n",
              source);
#endif // LC_DEBUG
      rc = MPI_Send(OBSOLETEVOIDPTR tmp_point->Data()->X()->data(),
                    tmp_point->Data()->X()->size(), MPI_DOUBLE, source,
                    asyncde::MPI_TAG_DOUBLE_ARRAY, MPI_COMM_WORLD);
    } else {
      // waiting for all workers to finish their last evaluation
      rank_active[source] = 0;
      active = 0;
      for (int dest_worker = 1; dest_worker < numtasks; dest_worker++)
        if (rank_active[dest_worker]) {
          active = 1;
          break;
        }

      if (!active)
        break;
    }
  }

  nevals = iterator->NFE();

  if (rc != 0) {
    // TODO:
  }

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

  for (auto point: buf_point)
    delete point;

  delete iterator;
  delete tmp_point;

  return retvalue;
}
