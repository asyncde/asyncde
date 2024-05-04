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

#ifndef ASYNCDE_ADEPOINTINFO_H
#define ASYNCDE_ADEPOINTINFO_H

#include <vector>

#include "asyncde/ADEConfig.h"
#include "asyncde/PointInfo.h"

namespace asyncde {

class ADEPointInfo : public PointInfo {
protected:
  /// Point's parent identifier
  long int parent_id;

public:
  /// F value used to generate this point (> 0); default = -1. Used by jDE
  double F;

  /// CR value used to generate this point (> 0); default = -1. Used by jDE
  double CR;

  /// threshold in correlation matrix [0; 1]; default = -1. Used by SCM/ACM
  double corr_thr;

  /// mask (blocks of coordinates inherited from a mutant vector). Used by
  /// ACM/Sum2
  std::vector<unsigned char> mask;

  /// choice of a competitor for the Darwinian selection
  ade_selection_types selectiontype;

public:
  ADEPointInfo(unsigned int nfreevars)
      : PointInfo(), parent_id(NOT_INITIALIZED_POINT), F(-1.0), CR(-1.0),
        corr_thr(-1.0), selectiontype(ADE_SELECTION_PARENT) {
    SetParentId(id);
    mask.resize(nfreevars);
    ResetMasks();
  }

  virtual ~ADEPointInfo() {}

  /// Set this point by contents of _pointinfo
  virtual int Set(const PointInfo &_pointinfo) override {
    int retvalue;

    if (0 != (retvalue = PointInfo::Set(_pointinfo)))
      return retvalue;

    const ADEPointInfo *adepointinfo =
        dynamic_cast<const ADEPointInfo *>(&_pointinfo);
    if (adepointinfo) {
      parent_id = adepointinfo->ParentId();
      F = adepointinfo->F;
      CR = adepointinfo->CR;
      corr_thr = adepointinfo->corr_thr;
      const std::vector<unsigned char> *_mask = adepointinfo->Mask();
      std::copy(_mask->begin(), _mask->end(), mask.begin());
      selectiontype = adepointinfo->selectiontype;
    } else
      ResetMasks();

    return retvalue;
  }

  /// Allocate a copy of this
  virtual PointInfo *Clone() const override {
    ADEPointInfo *newpointinfo = new ADEPointInfo(mask.size());
    newpointinfo->SetId(Id());
    newpointinfo->SetStatus(Status());
    newpointinfo->SetParentId(ParentId());
    newpointinfo->F = F;
    newpointinfo->CR = CR;
    newpointinfo->corr_thr = corr_thr;
    std::copy(mask.begin(), mask.end(), newpointinfo->MaskMutable()->begin());
    newpointinfo->selectiontype = selectiontype;

    return newpointinfo;
  }

  /// Reset point (set its fields to default values)
  virtual void Reset() override {
    PointInfo::Reset();
    parent_id = NOT_INITIALIZED_POINT;
    F = -1.0;
    CR = -1.0;
    corr_thr = -1.0;
    ResetMasks();
    selectiontype = ADE_SELECTION_PARENT;
  }

  /// Switch to a new unique ID
  virtual void SwitchToNewId() override {
    PointInfo::SwitchToNewId();
    if (parent_id < 0)
      parent_id = id;
  }

  /// Return the parent ID of the point
  long int ParentId() const { return parent_id; }

  /// Set the parent ID of the point
  void SetParentId(long int new_parent_id = NOT_INITIALIZED_POINT) {
    parent_id = (new_parent_id > -1) ? new_parent_id : id;
  }

  const std::vector<unsigned char> *Mask() const { return &mask; }

  std::vector<unsigned char> *MaskMutable() { return &mask; }

  /// reset masks to zeros
  void ResetMasks() { std::fill(mask.begin(), mask.end(), (unsigned char)0); }

  double CRfromMask(const std::vector<unsigned char> &corr_mask) const {
    size_t mutated = 0;
    size_t size = corr_mask.size();
    for (auto imask : corr_mask)
      if (imask > 0)
        mutated++;

    if (size < 1)
      return -1.0;
    return mutated / (double)size;
  }

  virtual int Print(FILE *stream = stdout) const override;
};

} // namespace asyncde

#endif
