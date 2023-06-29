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

#include "asyncde/ADEPointInfo.h"

#include <stdio.h>

int asyncde::ADEPointInfo::Print(FILE *stream) const {
  int retvalue = fprintf(stream, "asyncde::ADEPointInfo: ");
  retvalue |= PointInfo::Print(stream);
  if (retvalue < 0)
    return retvalue;

  retvalue = fprintf(stream, "parent_id=%li\n", parent_id);
  retvalue |= fprintf(stream, " mask:");
  for (auto imask : mask)
    retvalue |= fprintf(stream, " %i", (int)imask);
  retvalue |= fprintf(stream, "\n");

  return (retvalue < 0) ? -2 : 0;
}
