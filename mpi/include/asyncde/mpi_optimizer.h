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

#ifndef ASYNCDE_MPI_OPTIMIZER_H
#define ASYNCDE_MPI_OPTIMIZER_H

#include <vector>

#include "asyncde/Function.h"

namespace asyncde {

// forward declarations
class IteratorConfig;
class Problem;

enum MPI_TAG_NUM {
  MPI_TAG_STOP = 0,
  MPI_TAG_LONG_INT,
  MPI_TAG_DOUBLE_ARRAY,
};

int mpi_worker_cycle(const unsigned int Xsize, const unsigned int Ysize,
                     const asyncde::Functor1D &fitnessfunctor);

int mpi_master_cycle(const Problem &problem, const IteratorConfig &cfg,
                     long int &nevals, double &bestvalue,
                     std::vector<double> *bestX, std::vector<double> *parerrlow = 0,
                     std::vector<double> *parerrup = 0);

} // namespace asyncde

#endif
