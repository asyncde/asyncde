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

#ifndef ASYNCDE_ADERAND1ITERATOR_H
#define ASYNCDE_ADERAND1ITERATOR_H

#include "asyncde/AsyncIterator.h"

namespace asyncde {

// forward declarations:
class ADERand1Config;
class IteratorConfig;
class Point;
class Problem;

/*
  Class to implement the Asynchronous Differential Evolution with restart:
  - Population size Np is doubled with each restart
  - fixed F
  - No Crossover (CR=1)

  This class serves as the performance reference
  (the minimal time per iteration)
*/

class ADERand1Iterator : public AsyncIterator {
protected:
  /// link to AsyncIterator::cfg
  ADERand1Config *aderand1cfg;

  /// actual size of the base population (number of activated points)
  unsigned int base_population_actual_size;

  /// population (unsorted)
  std::vector<Point *> base_population;

  /// pointer to the best point in the current population
  Point *bestpointcpop;

  /// last point ID before restart
  long int lastidbeforerestart;

public:
  /// array [nfreeparams] which contains minimal coordinates in population
  std::vector<double> parlow;

  /// array [nfreeparams] which contains maximal coordinates in population
  std::vector<double> parup;

  /// array [nfreeparams] which contains an index of the point with the minimal
  /// coordinate in the population
  std::vector<int> parlowindex;

  /// array [nfreeparams] which contains an index of the point with the maximal
  /// coordinate in the population
  std::vector<int> parupindex;

protected:
  /// temporary vector
  std::vector<double> xexttmpvector; // [nparams]

  /// temporary vectors to store pseudo-random ints
  std::vector<unsigned int> nrand;
  std::vector<unsigned int> nrand_sorted;

public:
  ADERand1Iterator(const Problem &_problem, const IteratorConfig *_cfg);

  virtual ~ADERand1Iterator();

protected:
  /// resize local arrays to accomodate _nparents
  virtual int lcResize(unsigned int _nparents);

  /// Add point to Iterator; X-coordinates according to internal variables
  virtual int AddIntPoint(Point &_point) override;

  /// Iterator fills in a trial point;
  /// X-coordinates according to internal variables
  virtual int FillInTrialIntPoint(Point &_point, int maxnrestarts = 1) override;

  /// generate a point according to the uniform distribution (only if lower and
  /// upper bounds are specified)
  int GenerateUniformRandomPoint(Point &_point);

  /// generate a trial point according to the specified DE strategy
  virtual int FillInTrialIntADEPoint(Point &_point);

  /// return true and assign the tag if x coordinate "equals" to coordinates of
  /// a point in this list.
  int IsContainsX(const std::vector<double> &x, long int &_id,
                  Point **pointptr = nullptr) const;

public:
  /// return population size
  unsigned int PopulationSize(int *restartcounter = nullptr) const;

  /// restart algorithm with new population
  virtual int Restart(unsigned int _nparents) override;

  /// recovery algorithm for the population
  virtual int Recovery() override;

  /// update minimal x(y) epsilons in the population
  /// newpoint_index - the index of the replaced population member
  void UpdatePopulationSpread(const long int newpoint_index);

  /// return true if finished
  virtual int UpdateStatus() override;

  /// Return a read-only reference to the best point in the current population
  const Point *BestPointCurrentPopulation() const { return bestpointcpop; };
};

} // namespace asyncde

#endif
