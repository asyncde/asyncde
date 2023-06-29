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

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "asyncde/ADERand1Config.h"
#include "asyncde/ADERand1Iterator.h"
#include "asyncde/Rnd.h"

asyncde::AsyncIterator *
asyncde::ADERand1Config::NewIterator(const Problem &_problem) const {
  return new ADERand1Iterator(_problem, this);
}

double asyncde::ADERand1Config::GetF() {
  return (Fmax > Fmin) ? rnd->randCauchyTruncated(Fmu, Fsigma, Fmin, Fmax)
                       : Fmu;
}

int asyncde::ADERand1Config::Print(FILE *stream) const {
  int retvalue;

  if (!stream)
    return -1;

  fprintf(stream, "asyncde::ADERand1Config: ");
  retvalue = IteratorConfig::Print(stream);

  fprintf(stream,
          "minpopsize=%i  F=Cauchy(%.3e, %.3e) in [%.3e, %.3e]  verbose = %i\n",
          minpopsize, Fmu, Fsigma, Fmin, Fmax, verbose);

  return retvalue;
}
