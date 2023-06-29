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

#include "asyncde/StopCriteria.h"

int asyncde::StopCriteria::Print(FILE *stream, unsigned int mask) const {
  if (!stream)
    return -1;

  unsigned int activebits = statusbits & mask;
  fprintf(stream, "asyncde::StopCriteria: 0x%x\n", statusbits);
  if (0 != (activebits & CRITERION_vtr))
    fprintf(stream, "vtr = %.10e\n", vtr);
  if (0 != (activebits & CRITERION_nFE))
    fprintf(stream, "nFE = %li\n", nFE);
  if (0 != (activebits & CRITERION_popsizemax))
    fprintf(stream, "popsizemax = %i\n", popsizemax);
  if (0 != (activebits & CRITERION_nprojfeas))
    fprintf(stream, "nprojfeas = %i\n", nprojfeas);

  if (0 != (activebits & CRITERION_xepsilon))
    fprintf(stream, "xepsilon = %.3e\n", xepsilon);
  if (0 != (activebits & CRITERION_nxdiff))
    fprintf(stream, "nxdiff = %i\n", nxdiff);
  if (0 != (activebits & CRITERION_yepsilon))
    fprintf(stream, "yepsilon = %.3e\n", yepsilon);
  if (0 != (activebits & CRITERION_ydelta))
    fprintf(stream, "ydelta = %.3e\n", ydelta);

  if (0 != (activebits & CRITERION_failedtrialpoint))
    fprintf(stream, "failed to generate a trial point %u time(s)\n",
            nfailedtrialpoint);
  if (0 != (activebits & CRITERION_error))
    fprintf(stream, "error = %i\n", error);

  return 0;
}
