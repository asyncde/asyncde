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

#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#include "asyncde/ADEConfig.h"
#include "asyncde/mpi_optimizer.h"
#include "asyncde/ProblemMapping.h"
#include "asyncde/Rnd.h"
#include "asyncde/Variables.h"

// test functions
#include "asyncde/composite_rosenbrock.h"
#include "asyncde/rosenbrock.h"
#include "asyncde/sqr.h"

const unsigned int Ndim = 100;
const unsigned int Ysize = 1;

const double ftarget = -120.0;

int mpi_master_cycle() {
  int retvalue = 0;
  const double xmax = 10.0;
  const double xmax_initial = xmax;

  time_t t0 = time(0);
  asyncde::Rnd rnd(t0);
  
  asyncde::Variables *vars = new asyncde::Variables(Ndim);
  for (unsigned int ip = 0; ip < Ndim; ip++) {
    vars->SetName(ip, "");
    vars->SetCoordinate(ip, rnd.next(-xmax_initial, xmax_initial));
    //    vars->SetLowerLimit(ip, -xmax);
    //    vars->SetUpperLimit(ip, xmax);
    vars->SetInitialRange(ip, -xmax_initial, xmax_initial);
    vars->ReleaseVariable(ip);
  }
  vars->Print(stderr);
  asyncde::ProblemMapping problem;
  problem.SetExternalVariables(vars);
  problem.SetYsize(Ysize);

  // ADE:
  int Nparents = 10;
  asyncde::ADEConfig iterator_settings(&rnd, Nparents);
  iterator_settings.criteriastop.SetCriterion_nFE(1000000);
  iterator_settings.criteriastop.SetCriterion_vtr(ftarget + 1.e-6);

  long int nevals = 0;
  double bestvalue = std::numeric_limits<double>::max();
  double xdefault = -std::numeric_limits<double>::max();
  std::vector<double> bestX(Ndim, xdefault);

  retvalue = mpi_master_cycle(problem, iterator_settings, nevals, bestvalue,
                              &bestX);

  printf(
      "mpi_master_cycle():\nretvalue = %i  nevals = %li  bestvalue = %.15e\n",
      retvalue, nevals, bestvalue);

  int numtasks;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  for (int dest_worker = 1; dest_worker < numtasks; dest_worker++)
    retvalue |= MPI_Send(bestX.data(), 0, MPI_DOUBLE, dest_worker,
                         asyncde::MPI_TAG_STOP, MPI_COMM_WORLD);

  // bestX
  for (std::vector<double>::const_iterator itx = bestX.begin();
       itx != bestX.end(); itx++) {
    if (itx != bestX.begin())
      printf("  ");
    printf("%.15e", *itx);
  }
  printf("\n");

  delete vars;

  return retvalue;
}

int mpi_worker_cycle_call() {
  //  const asyncde::rosenbrock fitnessfunctor(ftarget);
  const asyncde::sqr fitnessfunctor(ftarget);
  return asyncde::mpi_worker_cycle(Ndim, Ysize, fitnessfunctor);
}

int main(int argc, char *argv[]) {
  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
    mpi_master_cycle();
  else
    mpi_worker_cycle_call();

  MPI_Finalize();

  return EXIT_SUCCESS;
}
