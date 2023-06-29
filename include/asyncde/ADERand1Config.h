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

#ifndef ASYNCDE_ADERAND1CONFIG_H
#define ASYNCDE_ADERAND1CONFIG_H

#include "asyncde/IteratorConfig.h"

namespace asyncde {

class ADERand1Config : public IteratorConfig {
protected:
  /// size of the population
  unsigned int minpopsize;

public:
  /// F according to a truncated Cauchy distrbution
  double Fmu;
  double Fsigma;
  double Fmin;
  double Fmax;

protected:
  void Init() {
    minpopsize = 10;
    Fmu = 0.9;
    Fsigma = 0.1;
    Fmin = 0.5;
    Fmax = 0.95;
  }

public:
  ADERand1Config(Rnd *_rnd, unsigned int _minpopsize) : IteratorConfig(_rnd) {
    Init();
    minpopsize = _minpopsize;
  }

  ADERand1Config(const asyncde::IteratorConfig &_cfg) : IteratorConfig(_cfg) {
    Init();
  }

  virtual ~ADERand1Config() {}

  virtual void Reset(unsigned int _minpopsize) override {
    IteratorConfig::Reset(_minpopsize);
    Init();
    minpopsize = _minpopsize;
  }

  virtual AsyncIterator *NewIterator(const Problem &_problem) const override;

  virtual void SetMinPopSize(unsigned int _minpopsize) {
    minpopsize = _minpopsize;
  }

  unsigned int MinPopSize() const { return minpopsize; }

  double GetF();

  /// return minimal population size for ADE
  unsigned int MinimalADEPopulationSize() const { return 3; }

  virtual int Print(FILE *stream = stdout) const override;
};

} // namespace asyncde

#endif
