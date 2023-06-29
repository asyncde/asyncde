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

#ifndef ASYNCDE_ADEARCHIVE_H
#define ASYNCDE_ADEARCHIVE_H

#include <vector>

namespace asyncde {

// forward declarations:
class ADEConfig;
class Point;

class ADEArchive {
protected:
  /// population sorted by id (descending order)
  std::vector<Point *> population;

  /// actual size of the archive population (number of activated points)
  unsigned int actualsize;

  ADEConfig *cfg; // not owner

public:
  ADEArchive(unsigned int _basepopsize, ADEConfig &_cfg);

  ~ADEArchive();

  void Reset() { actualsize = 0; }

  int Size() const { return actualsize; }

  int Resize(unsigned int _basepopsize);

  // insert a point into the archive
  int AddPoint(const Point &_point);

  // position = [0, ..., actualsize - 1]
  const Point *GetPoint(unsigned int position) const;
};

} // namespace asyncde

#endif
