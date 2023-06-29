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

#ifndef ASYNCDE_VARIABLES_H
#define ASYNCDE_VARIABLES_H

#include <limits>
#include <stdio.h>
#include <string>
#include <vector>

/**
 class to store input information for minimization:
   - vector of variables
   - limits on variables
 */

namespace asyncde {

class Variables {
protected:
  /// Maximal Variables'id used up to now (used to generate unique Variables'id)
  static long int id_counter;

  /// Unique Variables' identifier
  long int id;

  /// list of names
  std::vector<std::string> name;

  /// vector of variables
  std::vector<double> x;

  /// status of variables (free = 1 / fixed = 0)
  std::vector<int> freevar;

  /// lower limits for variables
  std::vector<double> minlimit;

  /// upper limits for variables
  std::vector<double> maxlimit;

  /// initial lower limits
  std::vector<double> mininitial;

  /// initial upper limits
  std::vector<double> maxinitial;

public:
  Variables(unsigned int nvariables = 0);

  Variables(const Variables &);

  virtual ~Variables() {}

  /// Return the highest Variables'id already used
  long int IdCounter() const { return id_counter; }

  /// Variables' id
  long int Id() const { return id; }

  /// Switch to a new unique ID
  void SwitchToNewId() {
    id_counter++;
    id = id_counter;
  }

public:
  /// Reset into an empty set of variables
  void Reset() {
    id = -1;
    Resize(0);
  }

  /// Reserve memory to contain nvars_newmaxsize variables
  void Reserve(unsigned int nvars_newmaxsize);

  /// Resize Variables to contain nvars_newsize variables
  /// added variables filled in by DefaultVariable()
  void Resize(unsigned int nvars_newsize);

  unsigned int NVariables() const { return x.size(); }

  /// return index of a new variable
  int AddVariable(const std::string _name, double value);

  /// copy variable with index extindex from extvars into variable with index
  /// intindex of this
  int CopyVariable(unsigned int intindex, const Variables &extvars,
                   int extindex);

  const std::vector<std::string> *Names() const { return &name; }
  const std::string Name(unsigned int ivar) const {
    if (ivar >= NVariables())
      return 0;
    return name[ivar];
  }
  int SetName(unsigned int ivar, const std::string _name);
  /// return index of the first variable named _name; -1 if not found
  int Find(const std::string _name) const;

  const std::vector<double> *X() const { return &x; }
  int SetCoordinate(unsigned int ivar, double value);
  int SetX(const std::vector<double> &xnew);

  const std::vector<int> *Mask() const { return &freevar; }
  int IsFreeVariable(unsigned int ivar) const {
    if (ivar >= NVariables())
      return 0;
    return freevar[ivar];
  }

  int FixVariable(unsigned int ivar) {
    if (ivar >= NVariables())
      return -1;
    freevar[ivar] = 0;
    return 0;
  }
  int ReleaseVariable(unsigned int ivar) {
    if (ivar >= NVariables())
      return -1;
    freevar[ivar] = 1;
    return 0;
  }

  unsigned int FindNFreeVariables() const;

  double DefaultLowerLimit() const {
    return -std::numeric_limits<double>::max();
  }

  double DefaultUpperLimit() const {
    return std::numeric_limits<double>::max();
  }

  const std::vector<double> *LowerLimits() const { return &minlimit; }
  int SetLowerLimit(unsigned int ivar, double value);

  const std::vector<double> *UpperLimits() const { return &maxlimit; }
  int SetUpperLimit(unsigned int ivar, double value);

  int SetLimits(const std::vector<double> &lowerlimits,
                const std::vector<double> &upperlimits);

  const std::vector<double> *InitialLowerLimits() const { return &mininitial; }

  const std::vector<double> *InitialUpperLimits() const { return &maxinitial; }

  int SetInitialRange(unsigned int ivar, double min, double max);

  int SetInitialLimits(const std::vector<double> &lowerlimits,
                       const std::vector<double> &upperlimits);

  int IsVariableConsistent(unsigned int ivar) const;

  // if failed (returns 0) index_failed is a first failed index; if failed_index
  // == NVariables() - wrong size of arrays
  int IsConsistent(int *index_failed = 0) const;

  // return 1 if both sets of variables are equivalent within tolerance
  int IsEquivalent(const Variables &_vars, double epsilon) const;

  int Print(FILE *stream = stdout) const;

protected:
  int DefaultVariable(unsigned int ivar);
};

} // namespace asyncde

#endif
