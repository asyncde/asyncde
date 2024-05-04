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

#include "asyncde/ADEArchive.h"

#include <algorithm>
#include <stdlib.h>

#include "asyncde/ADEConfig.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Rnd.h"

asyncde::ADEArchive::ADEArchive(unsigned int _basepopsize, ADEConfig &_cfg)
    : actualsize(0), cfg(&_cfg) {
  population.resize(_basepopsize * cfg->archivesizefactor, 0);
}

asyncde::ADEArchive::~ADEArchive() {
  for (std::vector<Point *>::iterator it = population.begin();
       it != population.end(); it++)
    delete *it;
}

int asyncde::ADEArchive::Resize(unsigned int _basepopsize) {
  unsigned int newsize = _basepopsize * cfg->archivesizefactor;

  if (newsize < population.size()) {
    for (std::vector<Point *>::iterator it = population.begin() + newsize;
         it != population.end(); it++)
      delete *it;
    if (actualsize > newsize)
      actualsize = newsize;
  }

  population.resize(newsize, 0);

  return 0;
}

// insert a point into the archive
int asyncde::ADEArchive::AddPoint(const asyncde::Point &_point) {
  Point *tmp_point_ptr;
  if (actualsize == population.size()) {
    // remove random point
    unsigned int irndpos = cfg->rnd->next_uniuint(actualsize);
    if (irndpos < actualsize - 1) {
      tmp_point_ptr = population[irndpos];
      std::move(population.begin() + irndpos + 1,
                population.begin() + actualsize, population.begin() + irndpos);
      population[actualsize - 1] = tmp_point_ptr;
    }
    actualsize--;
  }

  const long int newpointid = _point.Info()->Id();
  unsigned int inewpointpos = 0;
  for (; inewpointpos < actualsize; inewpointpos++)
    if (newpointid > population[inewpointpos]->Info()->Id())
      break;

  if (inewpointpos < actualsize) {
    tmp_point_ptr = population[actualsize];
    std::move_backward(population.begin() + inewpointpos,
                       population.begin() + actualsize,
                       population.begin() + actualsize + 1);
    population[inewpointpos] = tmp_point_ptr;
  }

  if (0 != population[inewpointpos])
    population[inewpointpos]->Set(_point);
  else
    population[inewpointpos] = _point.Clone();

  actualsize++;

  return 0;
}

// position = [0, ..., actualsize - 1]
const asyncde::Point *
asyncde::ADEArchive::GetPoint(unsigned int position) const {
  return population[std::min(position, actualsize)];
}
