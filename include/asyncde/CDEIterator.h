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

/**
 CDEIterator implements asynchronous variants
 of (Classic) Differential Evolution with the binomial crossover
 */

#ifndef ASYNCDE_CDEITERATOR_H
#define ASYNCDE_CDEITERATOR_H

#include <stdio.h>
#include <vector>

#include "asyncde/AsyncIterator.h"

namespace asyncde {

// forward declarations:
class ADEArchive;
class ADEConfig;
class ADEPoint;
class Point;
class Problem;
class Status;

class CDEIterator : public AsyncIterator {
protected:
  /// link to AsyncIterator::cfg
  ADEConfig *adecfg;

protected:
  /// number of free parameters (variables)
  unsigned int nfreeparams;

  /// actual size of the base population (number of activated points)
  unsigned int base_population_actual_size;

  /// population sorted by parent id (pid_key)
  std::vector<ADEPoint *> base_population;

  /// sorted array of parent ID keys (ascending)
  std::vector<long int> pid_key;

  /// sorted array of value keys (best is at [0])
  std::vector<double> value_key;

  /// parent ID (= base population) index for a value key
  std::vector<int> pid_index_for_value_key;

  /// value key index for a parent ID key
  std::vector<int> value_index_for_pid_key;

  /// last point ID before restart
  long int lastidbeforerestart;

  ADEArchive *archive;

public:
  /// array [nfreeparams] which contains minimal coordinates in population
  std::vector<double> parlow;

  /// array [nfreeparams] which contains maximal coordinates in population
  std::vector<double> parup;

  /// array [nfreeparams] which contains parent ID for a point with the minimal
  /// coordinate in the population
  std::vector<long int> parlowpid;

  /// array [nfreeparams] which contains parent ID for a point with the maximal
  /// coordinate in the population
  std::vector<long int> paruppid;

public:
  ///
  long int statsumsupdates_counter; // number of differential updates after the
                                    // last complete recalculation
  long int statsumsupdates_counter_max; // maximal allowed number of
                                        // differential updates

  /// for adaptive F
  double Fsum;
  double F2sum;
  double FmuJADE;

  /// for adaptive CR
  int CRn;
  double CRsum;
  double CR2sum;
  double CRmu;

protected:
  /// a temporary vector to keep used indices (value_key), sorted; used within
  /// generate_ade_trial_point()
  std::vector<int> vetovector;

  /// temporary vector
  std::vector<double> xexttmpvector; // [nparams]

  /// temporary vector. Used by CrossoverMaskUniform
  std::vector<int> tmpindices; // [nfreeparams]

public:
  CDEIterator(const Problem &_problem, const IteratorConfig *_cfg);

  virtual ~CDEIterator();

protected:
  /// resize local arrays to accomodate _nparents
  virtual int lcResize(unsigned int _nparents);

  virtual int FCrDefaultSettings();

protected:
  /// Iterator allocates a point specific to this iterator (internal variables)
  virtual Point *NewIntPoint() const override;

  /// return 1 if the _point is selected
  virtual int Selection(ADEPoint &_adepoint, int &ipidpos2replace);

  /// Add point to Iterator; X-coordinates according to internal variables
  /*
     return 1: if point is added to the population
     return 0: if point is not selected for the population
     return -1: if error
   */
  virtual int AddIntPoint(Point &_point) override;

  /// Iterator fills in a trial point; X-coordinates according to internal
  /// variables
  virtual int FillInTrialIntPoint(Point &_point, int maxnrestarts = 1) override;

public:
  /// Iterator allocates a point specific to this iterator (external variables)
  virtual Point *NewExtPoint() const override;

  /// return population size
  unsigned int PopulationSize(int *restartcounter = 0) const;

  /// restart algorithm with new population
  virtual int Restart(unsigned int _nparents) override;

  /// recovery algorithm for the population
  virtual int Recovery() override;

  /// Return a read-only reference to the best point in the current population
  const ADEPoint *BestIntPointCurrentPopulation() const;

  /// update minimal x(y) epsilons in the population
  /// iposition: index of the selected point in the population
  void UpdatePopulationSpread(const long int iposition);

  /// return true if finished
  virtual int UpdateStatus() override;

  const std::vector<ADEPoint *> *BasePopulation() const {
    return &base_population;
  }

  virtual int UpdateStatSums(const ADEPoint *new_point,
                             const ADEPoint *old_point);

  virtual int Print(FILE *stream = stdout) const override;

protected:
  /// generate a point according to the uniform distribution (only if lower and
  /// upper bounds are specified)
  int GenerateUniformRandomPoint(ADEPoint &_point);

  /// uniform (binomial) crossover, returns number of free coordinates selected
  /// from the mutant vector
  int CrossoverMaskUniform(const double CR, std::vector<char> &mask);

  /// set crossover mask for _point
  virtual int CrossoverMask(const ADEPoint &target_point, ADEPoint &_point);

  /// applies crossover operator according to mask,
  /// returns number of free coordinates selected from the mutant vector
  /// different == 1 if mutant vector is different from the target vector
  int CrossoverApplyMask(const std::vector<double> &target_vector,
                         std::vector<double> &mutant_vector,
                         const std::vector<char> &mask, int &different) const;

  /// Fills in the _point by the trial vector
  virtual int CrossoverIntADEPoint(const ADEPoint *target_point,
                                   ADEPoint &_point);

  virtual int AssignVector(const int vector_type,
                           const std::vector<double> *&vector,
                           const ADEPoint **_point = 0);

  int InitVetoVector();

  int InsertIntoVetoVector(int valueindex);

  /// fills in the _point by the mutant vector
  virtual int FillInMutantIntADEPoint(const ADEPoint *target_point,
                                      ADEPoint &_point);

  /// generates a trial point according to the specified DE strategy
  virtual int FillInTrialIntADEPoint(ADEPoint &_point);

  /// returns true and assign the tag if x coordinate "equals" to coordinates of
  /// a point in this list.
  int IsContainsX(const std::vector<double> &x, long int &_id,
                  Point **pointptr = 0) const;
};

} // namespace asyncde

#endif
