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

#ifndef ASYNCDE_MULTITHREADED_OPTIMIZER_H
#define ASYNCDE_MULTITHREADED_OPTIMIZER_H

#include <vector>

namespace asyncde {

// forward declarations
class IteratorConfig;
class Problem;

int multithreaded_optimization_cycle(const Problem &problem,
                                     const IteratorConfig &cfg,
                                     double &maxfunevals, double &bestvalue,
                                     std::vector<double> *bestX = nullptr,
                                     std::vector<double> *parerrlow = nullptr,
                                     std::vector<double> *parerrup = nullptr);

} // namespace asyncde

#endif // ASYNCDE_MULTITHREADED_OPTIMIZER_H
