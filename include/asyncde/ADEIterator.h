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
 ADEIterator implements Asynchronous Differential Evolution
 with Adaptive Correlation Matrix
 according to our work:
 Zhabitsky, M., Zhabitskaya, E. 2013.
 Asynchronous differential evolution with adaptive correlation matrix.
 In Proceedings of the 15th annual conference on Genetic and evolutionary
 computation (GECCO '13).
 Association for Computing Machinery, New York, NY, USA, 455–462.
 https://doi.org/10.1145/2463372.2463428
*/

#ifndef ASYNCDE_ADEITERATOR_H
#define ASYNCDE_ADEITERATOR_H

#include <stdio.h>
#include <vector>

#include "asyncde/ADEConfig.h"
#include "asyncde/AsyncIterator.h"
#include "asyncde/CDEIterator.h"

namespace asyncde {

// forward declarations:
class ADEConfig;
class ADEPoint;
class Point;
class Problem;

class ADEIterator : public CDEIterator {
public:
  /// correlation matrix components
  int corr_n;
  int corr_counter; // number of CM updates after the last complete update
  int corr_counter_maxfactor;   // times population size - maximal number of CM
                                // updates after the last complete update
  std::vector<double> corr_xi;  // [nfreeparams]
  std::vector<double> corr_xi2; // [nfreeparams]
  std::vector<double> corr_xij; // [nfreeparams, nfreeparams] use corr_index

  /// sample CM
  std::vector<double> scm_xii; // [nfreeparams]
  std::vector<double> scm_xij; // [nfreeparams, nfreeparams] use corr_index

  /// adaptive CM
  int acm_first;
  std::vector<double> acm_xii; // [nfreeparams]
  std::vector<double> acm_xij; // [nfreeparams, nfreeparams] use corr_index

  int corr_reliable;

public:
  ADEIterator(const Problem &_problem, const IteratorConfig *_cfg);

  virtual ~ADEIterator() {}

protected:
  virtual int FCrDefaultSettings() override;

public:
  unsigned int CorrXijIndex(unsigned int i, unsigned int j) const {
    return (i < j) ? nfreeparams * i + j - ((i + 1) * (i + 2)) / 2
                   : nfreeparams * j + i - ((j + 1) * (j + 2)) / 2;
  }

  virtual int UpdateStatSums(const ADEPoint *new_point,
                             const ADEPoint *old_point) override;

  int ResetSampleCorrMatrix();
  int ResetAdaptiveCorrMatrix();
  int UpdateSampleCorrMatrixSums(const std::vector<double> &x, double weight);
  int UpdateSampleCorrMatrixSumsByDiff(const std::vector<double> &xnew,
                                       const std::vector<double> &xold);
  int CalcSampleCorrMatrixFromSums();

  int CalculateSampleCorrMatrix();
  /// update Sample and Adaptive Correlation Matrices
  int UpdateCorrMatrices();

  /// calculate Adaptive Correlation Matrix
  int CalculateAdaptiveCorrMatrix();

  int PrintCorrMatrix() const;
  int TestCorrReliable() const;

  virtual int Print(FILE *fout = stdout) const override;

protected:
  /// crossover based on correlation matrix, returns number of free coordinates
  /// selected from the mutant vector
  int CrossoverMaskCorrMatrix(double corr_thr, unsigned int corr_index,
                              std::vector<unsigned char> &mask) const;

  /// set crossover mask for _point, returns number of coordinates
  /// selected from the mutant vector
  virtual int CrossoverMask(const ADEPoint &target_point,
                            ADEPoint &_point) override;
};

} // namespace asyncde

#endif
