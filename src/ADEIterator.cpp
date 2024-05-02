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

#include <stdio.h>
#include <stdlib.h>

#include "asyncde/ADEIterator.h"

#include "asyncde/ADEConfig.h"
#include "asyncde/ADEPoint.h"
#include "asyncde/ADEPointInfo.h"
#include "asyncde/Point.h"
#include "asyncde/PointData.h"
#include "asyncde/PointInfo.h"
#include "asyncde/Problem.h"
#include "asyncde/Rnd.h"
#include "asyncde/Variables.h"

//#define LC_DEBUG

//#define messages_stream stdout
#define messages_stream stderr
#define errors_stream stderr

asyncde::ADEIterator::ADEIterator(const Problem &_problem,
                                  const IteratorConfig *_cfg)
    : CDEIterator(_problem, _cfg), corr_n(0), corr_counter(0),
      corr_counter_maxfactor(2), acm_first(1), corr_reliable(0) {
  if (nfreeparams > 0) {
    corr_xi.resize(nfreeparams);
    corr_xi2.resize(nfreeparams);

    scm_xii.resize(nfreeparams);
    acm_xii.resize(nfreeparams);
    unsigned int maxsize = (nfreeparams * (nfreeparams - 1)) / 2;
    if (maxsize > 0) {
      corr_xij.resize(maxsize);
      scm_xij.resize(maxsize);
      acm_xij.resize(maxsize);
    }

    ADEIterator::FCrDefaultSettings();
  }
}

int asyncde::ADEIterator::FCrDefaultSettings() {
  CDEIterator::FCrDefaultSettings();

  ResetAdaptiveCorrMatrix();
  corr_reliable = 0;

  return 0;
}

int asyncde::ADEIterator::CrossoverMaskCorrMatrix(
    double corr_thr, unsigned int corr_index, std::vector<char> &mask) const {
  if (nfreeparams < 1)
    return -1;

  int nmutantcoords = nfreeparams;

  const std::vector<double> *cm_xij =
      (adecfg->crossovertype == ADE_CROSSOVER_SCM) ? &scm_xij : &acm_xij;

  for (unsigned int ivar = 0; ivar < nfreeparams; ivar++) {
    if (ivar == corr_index) {
      mask[ivar] = 1;
      continue;
    }

    if (fabs((*cm_xij)[CorrXijIndex(ivar, corr_index)]) < corr_thr) {
      mask[ivar] = 0;
      nmutantcoords--;
    } else
      mask[ivar] = 1;
  }

  return nmutantcoords;
}

int asyncde::ADEIterator::CrossoverMask(const ADEPoint &target_point,
                                        ADEPoint &_point) {
  if (!corr_reliable)
    corr_reliable = TestCorrReliable() > 0;

  if (adecfg->CRupdatetype != ADE_CROSSOVER_UPDATE_ACM || !corr_reliable)
    return CDEIterator::CrossoverMask(target_point, _point);

  _point.ADEInfoMutable()->corr_thr =
      (target_point.ADEInfo()->corr_thr >= 0.0 &&
       target_point.ADEInfo()->corr_thr <= 1.0 &&
       cfg->rnd->next(0.0, 1.0) > adecfg->tauCorr)
          ? target_point.ADEInfo()->corr_thr
          : cfg->rnd->next(0.0, 1.0);

  int corr_indexlocal = cfg->rnd->next_ulong(nfreeparams);

  int nmutantcoords =
      CrossoverMaskCorrMatrix(_point.ADEInfoMutable()->corr_thr,
                              corr_indexlocal, _point.ADEInfoMutable()->mask);

  _point.ADEInfoMutable()->CR =
      (nmutantcoords > 0) ? nmutantcoords / ((double)nfreeparams) : 0.0;

  return nmutantcoords;
}

int asyncde::ADEIterator::ResetSampleCorrMatrix() {
  unsigned int maxsize = (nfreeparams * (nfreeparams - 1)) / 2;

  corr_n = 0;
  corr_counter = 0;

  for (unsigned int i = 0; i < nfreeparams; i++) {
    corr_xi[i] = 0.0;
    corr_xi2[i] = 0.0;
    scm_xii[i] = 0.0;
  }
  for (unsigned int j = 0; j < maxsize; j++) {
    corr_xij[j] = 0.0;
    scm_xij[j] = 0.0;
  }

  return 0;
}

int asyncde::ADEIterator::ResetAdaptiveCorrMatrix() {
  unsigned int maxsize = (nfreeparams * (nfreeparams - 1)) / 2;

  acm_first = 1;
  for (unsigned int i = 0; i < nfreeparams; i++)
    acm_xii[i] = 0.0;

  for (unsigned int j = 0; j < maxsize; j++)
    acm_xij[j] = 0.0;

  return ResetSampleCorrMatrix();
}

int asyncde::ADEIterator::UpdateStatSums(const ADEPoint *new_point,
                                         const ADEPoint *old_point) {
  int retvalue = CDEIterator::UpdateStatSums(new_point, old_point);

  if (new_point)
    retvalue |= (old_point != 0)
                    ? UpdateSampleCorrMatrixSumsByDiff(*new_point->Data()->X(),
                                                       *old_point->Data()->X())
                    : UpdateSampleCorrMatrixSums(*new_point->Data()->X(), 1.0);

  if (old_point) // population is complete
    UpdateCorrMatrices();

  return retvalue;
}

int asyncde::ADEIterator::UpdateSampleCorrMatrixSums(
    const std::vector<double> &x, double weight) {
  for (unsigned int i = 0; i < nfreeparams; i++) {
    double v = x[i];
    double vw = v * weight;
    corr_xi[i] += vw;
    corr_xi2[i] += v * vw;

    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      corr_xij[index] += x[j] * vw;
      index++;
    }
  }

  if (weight < 0.0)
    corr_n--;
  else
    corr_n++;

  corr_counter++;

  return 0;
}

int asyncde::ADEIterator::UpdateSampleCorrMatrixSumsByDiff(
    const std::vector<double> &xnew, const std::vector<double> &xold) {
  for (unsigned int i = 0; i < nfreeparams; i++) {
    double vnew = xnew[i];
    double vold = xold[i];
    double d = vnew - vold;
    corr_xi[i] += d;
    corr_xi2[i] += d * (vnew + vold);

    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      corr_xij[index] += xnew[j] * vnew - xold[j] * vold;
      index++;
    }
  }

  corr_counter++;

  return 0;
}

int asyncde::ADEIterator::CalcSampleCorrMatrixFromSums() {
  if (corr_n < 2) {
    for (unsigned int i = 0; i < nfreeparams; i++) {
      scm_xii[i] = 0.0;
      unsigned int index = CorrXijIndex(i, i + 1);
      for (unsigned j = i + 1; j < nfreeparams; j++) {
        scm_xij[index] = 0.0;
        index++;
      }
    }
    return -1;
  }

  double denom = 1.0 / sqrt(corr_n * (corr_n - 1));
  for (unsigned int i = 0; i < nfreeparams; i++) {
    const double v = corr_n * corr_xi2[i] - corr_xi[i] * corr_xi[i];
    scm_xii[i] = (v > 0.0) ? denom * sqrt(v) : 0.0;
  }

  for (unsigned int i = 0; i < nfreeparams; i++) {
    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      double v = (corr_n - 1.0) * scm_xii[i] * scm_xii[j];
      if (v > 0.0) {
        v = (corr_xij[index] - corr_xi[i] * corr_xi[j] / corr_n) / v;
        if (v > 1.0)
          v = 1.0;
        else if (v < -1.0)
          v = -1.0;
      } else
        v = 0.0;
      scm_xij[index] = v;
      index++;
    }
  }

  return 0;
}

int asyncde::ADEIterator::CalculateSampleCorrMatrix() {
  ResetSampleCorrMatrix();

  for (unsigned int ipop = 0; ipop < base_population_actual_size; ipop++)
    UpdateSampleCorrMatrixSums(*base_population[ipop]->Data()->X(), 1.0);

  return CalcSampleCorrMatrixFromSums();
}

int asyncde::ADEIterator::PrintCorrMatrix() const {
  printf("SCM:\n");
  for (unsigned int i = 0; i < nfreeparams; i++) {
    for (unsigned int j = 0; j < i; j++)
      printf("------------ ");

    printf("%12.5e ", scm_xii[i]);

    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      printf("%12.5e ", scm_xij[index]);
      index++;
    }
    printf("\n");
  }

  printf("ACM:\n");
  for (unsigned int i = 0; i < nfreeparams; i++) {
    for (unsigned int j = 0; j < i; j++)
      printf("------------ ");

    printf("%12.5e ", acm_xii[i]);

    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      printf("%12.5e ", acm_xij[index]);
      index++;
    }
    printf("\n");
  }

  return 0;
}

int asyncde::ADEIterator::UpdateCorrMatrices() {
  int retvalue = (corr_counter <
                  (int)(corr_counter_maxfactor * base_population_actual_size))
                     ? CalcSampleCorrMatrixFromSums()
                     : CalculateSampleCorrMatrix();

  retvalue |= CalculateAdaptiveCorrMatrix();

  return retvalue;
}

int asyncde::ADEIterator::CalculateAdaptiveCorrMatrix() {
  int retvalue = 0;

  const double wCR = (acm_first) ? 1.0 : adecfg->wCR;

  for (unsigned int i = 0; i < nfreeparams; i++) {
    acm_xii[i] += wCR * (scm_xii[i] - acm_xii[i]);

    unsigned int index = CorrXijIndex(i, i + 1);
    for (unsigned int j = i + 1; j < nfreeparams; j++) {
      acm_xij[index] += wCR * (scm_xij[index] - acm_xij[index]);
      index++;
    }
  }

  acm_first = 0;

  return retvalue;
}

int asyncde::ADEIterator::TestCorrReliable() const {
  unsigned int ipop;
  unsigned int modified = 0;

  if (base_population_actual_size != adecfg->MinPopSize())
    return -1;

  for (ipop = 0; ipop < base_population_actual_size; ipop++)
    if (base_population[ipop]->Info()->Id() >
        base_population[ipop]->ADEInfo()->ParentId())
      modified++;

  return (modified >= base_population_actual_size * adecfg->pcorr_reliable) ? 1
                                                                            : 0;
}

int asyncde::ADEIterator::Print(FILE *stream) const {
  int retvalue = adecfg->Print(stream);

  // print status (progress/stopped)

  // print best so far value

  return retvalue;
}
