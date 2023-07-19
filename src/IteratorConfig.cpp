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

#include "asyncde/IteratorConfig.h"

int asyncde::IteratorConfig::Print(FILE *stream) const {
  fprintf(stream, "asyncde::IteratorConfig:\n");

  if (criteriastop.statusbits) {
    fprintf(stream, "Stop criteria: ");
    criteriastop.Print(stream);
  }

  if (criteriarestart.statusbits) {
    fprintf(stream, "Restart criteria: ");
    criteriarestart.Print(stream);
  }

  if (criteriarecover.statusbits) {
    fprintf(stream, "Recover criteria: ");
    criteriarecover.Print(stream);
  }

  if (criteriastop.statusbits || criteriarestart.statusbits ||
      criteriarecover.statusbits)
    fprintf(stream,
            "The leading %.2e part of the population is taken into account for "
            "CRITERION_yepsilon and CRITERION_ydelta\n",
            pybest);

  return 0;
}
