// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_OPTIQUANTALGORITHM_H
#define OPENMS_ANALYSIS_QUANTITATION_OPTIQUANTALGORITHM_H

#include <utility>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>

namespace OpenMS
{

/**
  @brief Perform label-free quantification.

  Perform label-free quantification.

  @htmlinclude OpenMS_OptiQuant.parameters

  @ingroup Quantitation
*/

class OPENMS_DLLAPI OptiQuantAlgorithm :
    public DefaultParamHandler,
    public ProgressLogger
{
public:

  class FeatureHypothesis
  {

  public:

    FeatureHypothesis(Int charge) :
      charge_(charge),
      masstraces_()
    {
    }

    /// Copy constructor
    FeatureHypothesis(const FeatureHypothesis& rhs) :
      charge_(rhs.charge_),
      masstraces_(rhs.masstraces_)
    {
    }

    /// Assignment operator
    FeatureHypothesis& operator=(FeatureHypothesis const& rhs)
    {
      charge_ = rhs.charge_;
      masstraces_ = rhs.masstraces_;

      return *this;
    }

    /// Destructor
    virtual ~FeatureHypothesis()
    {
    }

    const std::vector<std::pair<Size, Size> >& getMassTraces() const
    {
      return masstraces_;
    }

    void addMassTrace(const std::pair<Size, Size>& mt)
    {
      masstraces_.push_back(mt);
    }

    Int getCharge() const
    {
      return charge_;
    }

    Size size() const
    {
      return masstraces_.size();
    }

  protected:

    Int charge_;
    std::vector<std::pair<Size, Size> > masstraces_; // pairs iso_pos / mass trace index

  };

  /// Constructor
  OptiQuantAlgorithm(const ConsensusMap& input_map = ConsensusMap(), Int num_threads = 1);

  /// Default destructor
  virtual ~OptiQuantAlgorithm();

  /// Main method of the alorithm
  void run(ConsensusMap& output_map);

protected:

  /// Maximum number of parallel threads CPLEX is allowed to use
  Int threads_;

  /// Number of maps
  Size num_maps_;

  /// TODO
  double rt_tol_secs_;

  /// TODO
  double mz_tol_;

  /// TODO
  bool mz_ppm_;

  /// TODO
  Int charge_low_;

  /// TODO
  Int charge_high_;

  /// TODO
  Size max_nr_traces_;

  /// TODO
  double min_int_score_thresh_;

  /// TODO
  Size require_first_n_traces_;

  /// TODO
  Size min_nr_traces_per_map_;

  /// TODO
  bool require_monoiso_;

  /// TODO
  bool use_ids_;

  /// TODO
  Int solver_time_limit_;

  /// TODO
  bool include_unassembled_traces_;

  /// TODO
  bool include_unidentified_unassembled_traces_;

  /// TODO
  UInt score_size_exp_;

  /// TODO
  UInt score_int_exp_;

  /// TODO
  double score_int_weight_;

  /// TODO
  UInt score_mz_exp_;

  /// TODO
  double score_mz_weight_;

  /// TODO
  UInt score_rt_exp_;

  /// TODO
  double score_rt_weight_;

  /// TODO
  double score_denom_;

  /// TODO
  double score_id_weight_;

  /// TODO
  bool score_int_ignore_missing_;

  /// TODO
  Size quantify_top_;

  /// TODO
  KDTreeFeatureMaps kd_data_;

  /// TODO
  const ConsensusMap* input_map_;

  /// TODO
  std::vector<Int> mt_assembled_;

  /// TODO
  void assembleFeatures_(std::vector<FeatureHypothesis>& features);

  /// TODO
  void addHypotheses_(Size mono_iso_mt_index, const std::vector<Size>& candidate_mts, std::vector<FeatureHypothesis>& hypos, std::vector<std::vector<Size> >& hypos_for_mt);

  /// TODO
  void resolveConflicts_(const std::vector<FeatureHypothesis>& hypos, const std::vector<std::vector<Size> >& hypos_for_mt, std::vector<FeatureHypothesis>& result);

  /// TODO
  void resolveHypothesisCluster_(const std::vector<FeatureHypothesis>& hypos, const std::vector<std::vector<Size> >& hypos_for_mt, const std::set<Size>& hypo_cluster_indices, std::vector<FeatureHypothesis>& result);

  /// TODO
  double averagineCorrelation_(const std::vector<std::pair<Size, double> >& hypo_int_pairs, const double& mol_weight, bool ignore_missing = false) const;

  /// TODO
  double computeMZScore_(const FeatureHypothesis& hypo) const;

  /// TODO
  double computeRTScore_(const FeatureHypothesis& hypo) const;

  /// TODO
  double computeIntensityScore_(const FeatureHypothesis& hypo, bool consensus_only = false) const;

  /// TODO
  double computeScore_(const FeatureHypothesis& hypo) const;

  /// TODO
  void compileResults_(const std::vector<FeatureHypothesis>& features, ConsensusMap& output_map);

  /// TODO
  void outputStatistics_(const ConsensusMap& cmap) const;

  virtual void updateMembers_();

  /// TODO
  struct IsoTraceTuple
  {
    Size iso_pos;
    Size size;
    double intensity;

    bool operator<(const IsoTraceTuple& rhs) const
    {
      if (size > rhs.size)
      {
        return true;
      }
      else if (size == rhs.size)
      {
        return intensity > rhs.intensity;
      }

      return false;
    }
  };
};

}

#endif // OPENMS_ANALYSIS_QUANTITATION_OPTIQUANTALGORITHM_H
