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

#include "asyncde/ADEConfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <string.h>

#include "asyncde/ADEIterator.h"

asyncde::AsyncIterator *
asyncde::ADEConfig::NewIterator(const Problem &_problem) const {
  return new ADEIterator(_problem, this);
}

int asyncde::ADEConfig::SetStrategy(const std::string &strategyname) {
  int target_vector_choice = -1;
  int base_vector_choice = -1;
  int ndifferences_choice = -1;
  ade_crossover_types crossover_choice = ADE_CROSSOVER_UNDEFINED;

  int retvalue = 0;
  int slashcount = 1;
  const char *start = strategyname.c_str();
  const char *next;
  int len;

  if (!(next = strchr(start, '/')))
    return -2;

  if ((next - start != 2) || (strncasecmp(start, "DE", 2) != 0))
    return -2;

  next++;
  start = next;
  while (*next) {
    if (*next == '/')
      slashcount++;
    next++;
  }
  if (slashcount < 3 || slashcount > 4)
    return -2;

  while (slashcount > 0) {
    if (!(next = (slashcount > 1) ? strchr(start, '/') : strchr(start, '\0')))
      return -2;

    switch (slashcount) {
    case 4:
      // target vector (asynchronous DE)
      do {
#undef ADE_VECTOR
#define ADE_VECTOR(ID, NAME)                                                   \
  len = strlen(NAME);                                                          \
  if ((next - start == len) && (strncasecmp(start, NAME, len) == 0)) {         \
    target_vector_choice = ADE_VECTOR_##ID;                                    \
    break;                                                                     \
  }

        ADE_VECTORS

        break;
      } while (0);

      if (target_vector_choice < 0)
        retvalue = -4;
      break;
    case 3:
      // base vector
      do {
#undef ADE_VECTOR
#define ADE_VECTOR(ID, NAME)                                                   \
  len = strlen(NAME);                                                          \
  if ((next - start == len) && (strncasecmp(start, NAME, len) == 0)) {         \
    base_vector_choice = ADE_VECTOR_##ID;                                      \
    break;                                                                     \
  }

        ADE_VECTORS

        break;
      } while (0);

      if (base_vector_choice < 0)
        retvalue = -4;
      break;
    case 2:
      // number of difference vectors
      ndifferences_choice = atoi(start);
      if (ndifferences_choice < 1)
        retvalue = -4;
      break;
    case 1:
      // crossover type
      do {
#undef ADE_CROSSOVER
#define ADE_CROSSOVER(ID, NAME)                                                \
  len = strlen(NAME);                                                          \
  if (strncasecmp(start, NAME, len) == 0) {                                    \
    crossover_choice = ADE_CROSSOVER_##ID;                                     \
    break;                                                                     \
  }

        ADE_CROSSOVER_TYPES

        break;
      } while (0);

      if (crossover_choice <= ADE_CROSSOVER_UNDEFINED)
        retvalue = -4;
      break;
    default:
      return -2;
    }

    slashcount--;
    next++;
    start = next;
  }

  if (retvalue == 0) {
    targetvector = target_vector_choice;
    if (targetvector == -1)
      targetvector = ADE_VECTOR_RAND;
    basevector = base_vector_choice;
    ndifferences = ndifferences_choice;
    crossovertype = crossover_choice;
    switch (crossovertype) {
    case ADE_CROSSOVER_ACM:
    case ADE_CROSSOVER_SCM:
      CRupdatetype = ADE_CROSSOVER_UPDATE_ACM;
      break;
    case ADE_CROSSOVER_BIN:
    default:
      CRupdatetype = ADE_CROSSOVER_UPDATE_Cauchy;
    }

    if (MinPopSize() < MinimalADEPopulationSize())
      SetMinPopSize(MinimalADEPopulationSize());
  }

  return retvalue;
}

std::string asyncde::ADEConfig::ComposeStrategyName() const {
  std::string strategyname = "DE/";

  // target vector (asynchronous DE)
  for (;;) {
#undef ADE_VECTOR
#define ADE_VECTOR(ID, NAME)                                                   \
  if (targetvector == ADE_VECTOR_##ID) {                                       \
    strategyname += NAME;                                                      \
    strategyname += "/";                                                       \
    break;                                                                     \
  }

    ADE_VECTORS

    strategyname += "UNKNOWN/";
    break;
  }

  // base vector
  for (;;) {
#undef ADE_VECTOR
#define ADE_VECTOR(ID, NAME)                                                   \
  if (basevector == ADE_VECTOR_##ID) {                                         \
    strategyname += NAME;                                                      \
    strategyname += "/";                                                       \
    break;                                                                     \
  }

    ADE_VECTORS

    strategyname += "UNKNOWN/";
    break;
  }

  // number of difference vectors
  strategyname += std::to_string(ndifferences);
  strategyname += "/";

  // crossover type
  for (;;) {
#undef ADE_CROSSOVER
#define ADE_CROSSOVER(ID, NAME)                                                \
  if (crossovertype == ADE_CROSSOVER_##ID) {                                   \
    strategyname += NAME;                                                      \
    break;                                                                     \
  }

    ADE_CROSSOVER_TYPES

    strategyname += "UNKNOWN";
    break;
  }

  return strategyname;
}

int asyncde::ADEConfig::Print(FILE *stream) const {
  int retvalue;

  if (!stream)
    return -1;

  fprintf(stream, "asyncde::ADEConfig: ");
  retvalue = IteratorConfig::Print(stream);

  fprintf(stream, "minpopsize = %i  verbose = %i\n", minpopsize, verbose);

  std::string strategyname = ComposeStrategyName();
  retvalue |= fprintf(stream, "Strategy: %s,", strategyname.c_str());
  if (basevector == ADE_VECTOR_CURRENTTOPBEST)
    retvalue |= fprintf(stream, " pb=%.3e,", pbest);
  if (targetvector == ADE_VECTOR_PWORST)
    retvalue |= fprintf(stream, " pw=%.3e,", pworst);
  retvalue |= fprintf(stream, " Npmin=%i\n", minpopsize);

  if (0 > retvalue)
    return retvalue;

  switch (Fscaletype) {
  case ADE_FSCALE_jDE:
    if (0 > (retvalue = fprintf(stream, "F type=%i in [%.5e, %.5e] tauF=%.5e\n",
                                Fscaletype, Fmin, Fmax, tauF)))
      return retvalue;
    break;
  case ADE_FSCALE_JADE:
    if (0 >
        (retvalue = fprintf(
             stream,
             "F=Cauchy(%.5e, %.5e) in [%.5e, %.5e] Fmu adapted with cF=%.5e\n",
             Fmu, Fsigma, Fmin, Fmax, Fc)))
      return retvalue;
    break;
  case ADE_FSCALE_Cauchy:
  default:
    if (0 >
        (retvalue = fprintf(stream, "F=Cauchy(%.5e, %.5e) in [%.5e, %.5e]\n",
                            Fmu, Fsigma, Fmin, Fmax)))
      return retvalue;
  }

  switch (CRupdatetype) {
  case ADE_CROSSOVER_UPDATE_ACM:
    if (0 > (retvalue = fprintf(stream, "CR CM: tauCorr=%.5e wCR=%.5e\n",
                                tauCorr, wCR)))
      return retvalue;
    break;
  case ADE_CROSSOVER_UPDATE_JADE:
    if (0 > (retvalue =
                 fprintf(stream, "CR=N(%.5e, %.5e) in [%.5e, %.5e] cCR=%.5e\n",
                         CRmu, CRsigma, CRmin, CRmax, CRc)))
      return retvalue;
    break;
  case ADE_CROSSOVER_UPDATE_jDE:
    if (0 > (retvalue = fprintf(stream, "CR in [%.5e, %.5e] tauCR=%.5e\n",
                                CRmin, CRmax, tauCR)))
      return retvalue;
    break;
  case ADE_CROSSOVER_UPDATE_Cauchy:
  default:
    if (0 >
        (retvalue = fprintf(stream, "CR=Cauchy(%.5e, %.5e) in [%.5e, %.5e]\n",
                            CRmu, CRsigma, CRmin, CRmax)))
      return retvalue;
    break;
  }

  switch (selectiontype) {
  case ADE_SELECTION_RANDOM:
    if (0 > (retvalue = fprintf(
                 stream, "Darwinian (1+1)-ES=(random+trial)-ES selection\n")))
      return retvalue;
    break;
  case ADE_SELECTION_WORST:
    if (0 > (retvalue = fprintf(
                 stream, "Darwinian (Np+1)-ES=(Np+trial)-ES selection\n")))
      return retvalue;
    break;
  case ADE_SELECTION_PARENT:
  default:
    if (0 > (retvalue = fprintf(
                 stream, "Darwinian (1+1)-ES=(parent+trial)-ES selection\n")))
      return retvalue;
  }

  if (archive)
    if (0 >
        (retvalue = fprintf(stream, "With archive, archivesizefactor=%.2e\n",
                            archivesizefactor)))
      return retvalue;

  switch (mvector_unique) {
  case ADE_MVECTOR_UNIQUE_NOTEST:
    if (0 >
        (retvalue = fprintf(stream, "Mutant vector: not tested as unique\n")))
      return retvalue;
    break;
  case ADE_MVECTOR_DIFF_TARGET:
    if (0 > (retvalue = fprintf(
                 stream,
                 "Mutant vector: different from base and target vectors\n")))
      return retvalue;
    break;
  case ADE_MVECTOR_UNIQUE_INPOP:
    if (0 >
        (retvalue = fprintf(stream, "Mutant vector: unique in population\n")))
      return retvalue;
    break;
  default:
    if (0 > (retvalue = fprintf(stream, "Mutant vector: %i\n", mvector_unique)))
      return retvalue;
  }

  return 0;
}
