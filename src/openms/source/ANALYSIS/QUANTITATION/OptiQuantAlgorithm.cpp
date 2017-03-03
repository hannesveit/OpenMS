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

#include <OpenMS/ANALYSIS/QUANTITATION/OptiQuantAlgorithm.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <assert.h>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/container/vector.hpp>
#include <boost/unordered_set.hpp>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using namespace std;

namespace OpenMS
{

OptiQuantAlgorithm::OptiQuantAlgorithm(const ConsensusMap& input_map, Int num_threads) :
  DefaultParamHandler("OptiQuantAlgorithm"),
  ProgressLogger()
{
  input_map_ = &input_map;
  threads_ = num_threads;
  num_maps_ = input_map.getFileDescriptions().size();
  kd_data_.addFeatures(0, input_map, true);
  mt_assembled_ = vector<Int>(input_map.size(), false);

  defaults_.setValue("mz_tol", 10.0, "m/z tolerance (in Da or ppm)");
  defaults_.setMinFloat("mz_tol", 0.0);

  defaults_.setValue("mz_unit", "ppm", "unit of m/z tolerance");
  defaults_.setValidStrings("mz_unit", ListUtils::create<String>("ppm,Da"));

  defaults_.setValue("rt_tol", 5.0, "RT tolerance (in seconds)");
  defaults_.setMinFloat("rt_tol", 0.0);

  defaults_.setValue("charge_low", 1, "Lowest charge state to consider");
  defaults_.setMinInt("charge_low", 1);

  defaults_.setValue("charge_high", 5, "Highest charge state to consider");
  defaults_.setMinInt("charge_high", 1);

  defaults_.setValue("min_averagine_score", 0.85, "Minimum averagine similarity score a hypothesis must achieve in order to be considered.");
  defaults_.setMinFloat("min_averagine_score", 0.0);
  defaults_.setMaxFloat("min_averagine_score", 1.0);

  defaults_.setValue("max_nr_traces", 6, "Consider only the first max_nr_traces isotope traces");
  defaults_.setMinInt("max_nr_traces", 1);

  defaults_.setValue("quantify_top", 3, "In final intensity calculation, quantify only the top n most abundant (found across most maps) isotopic consensus traces per hypothesis. Ties are broken by abundance. If fewer traces are present, all of them are quantified. If set to 0, all traces are quantified.)");
  defaults_.setMinInt("quantify_top", 0);

  defaults_.setValue("require_first_n_traces", 3, "Do not consider consensus feature hypotheses in which any of the first n isotope traces are missing across all maps (including the monoisotopic trace)");
  defaults_.setMinInt("require_first_n_traces", 1);

  defaults_.setValue("min_nr_traces_per_map", 2, "Ignore subfeatures with less than this many detected mass traces");
  defaults_.setMinInt("min_nr_traces_per_map", 1);

  defaults_.setValue("keep_unassembled_traces", "identified", "Include unassembled traces in the results? When set to 'identified', keep only those traces annotated with at least one peptide identification");
  defaults_.setValidStrings("keep_unassembled_traces", ListUtils::create<String>("all,identified,none"));

  // advanced:

  defaults_.setValue("adaptive_iso_mass_diff", "true", "Re-estimate isotopic mass difference while collecting isotopic traces", ListUtils::create<String>("advanced"));
  defaults_.setValidStrings("adaptive_iso_mass_diff", ListUtils::create<String>("true,false"));

  defaults_.setValue("require_monoiso", "true", "Include subfeature for map i only if the monoisotopic trace of this feature has been detected in map i", ListUtils::create<String>("advanced"));
  defaults_.setValidStrings("require_monoiso", ListUtils::create<String>("true,false"));

  defaults_.setValue("use_ids", "true", "If a mass trace has identifications attached, generate only hypotheses for the charge state(s) found in these peptide IDs when generating hypotheses for this (monoisotopic) mass trace.", ListUtils::create<String>("advanced"));
  defaults_.setValidStrings("use_ids", ListUtils::create<String>("true,false"));

  defaults_.setValue("score:size_exp", 2, "Exponent of the size component (number of isotope traces) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinInt("score:size_exp", 1);

  defaults_.setValue("score:int_exp", 2, "Exponent of the intensity component (averagine similarity) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinInt("score:int_exp", 1);

  defaults_.setValue("score:int_weight", 2.0, "Weighting factor of the intensity component (averagine similarity) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinFloat("score:int_weight", 0.0);

  defaults_.setValue("score:mz_exp", 1, "Exponent of the m/z component (relative MAD of reconstructed monoisotopic m/z values based on all isotope traces) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinInt("score:mz_exp", 1);

  defaults_.setValue("score:mz_weight", 1.0, "Weighting factor of the m/z component (relative MAD of reconstructed monoisotopic m/z values based on all isotope traces) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinFloat("score:mz_weight", 0.0);

  defaults_.setValue("score:rt_exp", 1, "Exponent of the RT component (relative MAD of RT values of all isotope traces) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinInt("score:rt_exp", 1);

  defaults_.setValue("score:rt_weight", 1.0, "Weighting factor of the RT component (relative MAD of RT values of all isotope traces) in the overall hypothesis score", ListUtils::create<String>("advanced"));
  defaults_.setMinFloat("score:rt_weight", 0.0);

  defaults_.setValue("score:id_weight", 1.0, "Weighting factor applied to final score for hypotheses in which the monoisotopic mass trace is annotated with an identification. When set to 1, the ID is ignored and the hypothesis is scored as if it weren't there.", ListUtils::create<String>("advanced"));
  defaults_.setMinFloat("score:id_weight", 1.0);

  defaults_.setValue("score:int_ignore_missing", "false", "If true, averagine correlation is computed using only the isotopic traces that were found. If false, zero values are inserted at the missing positions.", ListUtils::create<String>("advanced"));
  defaults_.setValidStrings("score:int_ignore_missing", ListUtils::create<String>("true,false"));

  defaults_.setSectionDescription("score", "Parameters of the hypothesis scoring function: s = (s_size)^(e_size) * (w_int * (s_int)^(e_int) + w_mz * (s_mz)^(e_mz) + w_rt * (s_rt)^(e_rt)) / (w_int + w_mz + w_rt).");

  defaults_.setValue("solver_time_limit", -1, "CPLEX time limit (in seconds) for solving one cluster of contiguous hypotheses. No time limit when set to -1.", ListUtils::create<String>("advanced"));

  defaultsToParam_();

  this->setLogType(CMD);
}

OptiQuantAlgorithm::~OptiQuantAlgorithm()
{
}

void OptiQuantAlgorithm::run(ConsensusMap& output_map)
{
  // set parameters here instead of constructor (param_ is now set)
  kd_data_.setParameters(param_);

  // assemble consensus features
  vector<FeatureHypothesis> features;
  assembleFeatures_(features);

  // fill result map
  compileResults_(features, output_map);

  // ouput some statistics
  outputStatistics_(output_map);
}

void OptiQuantAlgorithm::assembleFeatures_(vector<FeatureHypothesis>& features)
{
  // somewhat hacky (see below)
  double epsilon = 10e-7;

  // compute m/z window size
  double mz_win_height = Constants::C13C12_MASSDIFF_U * (double)max_nr_traces_ / (double)charge_low_;

  // collect feature hypotheses here
  vector<FeatureHypothesis> hypos;
  vector<vector<Size> > hypos_for_mt(kd_data_.size(), vector<Size>());

  UInt progress = 0;
  startProgress(0, kd_data_.size(), "generating hypotheses");
  for (Size i = 0; i < kd_data_.size(); ++i)
  {
    // reference m/z and RT of monoisotopic mass trace
    double ref_mz = kd_data_.mz(i);
    double ref_rt = kd_data_.rt(i);

    double rt_low = ref_rt - rt_tol_secs_ / 2;
    double rt_high = ref_rt + rt_tol_secs_ / 2;
    double mz_low = ref_mz + epsilon; // hacky: exclude reference masstrace itself by adding epsilon
    double mz_high = ref_mz + mz_win_height;
    vector<Size> candidate_indices;
    kd_data_.queryRegion(rt_low, rt_high, mz_low, mz_high, candidate_indices);

    // add all hypotheses for monoisotopic masstrace i to hypos
    addHypotheses_(i, candidate_indices, hypos, hypos_for_mt);
    setProgress(++progress);
  }
  endProgress();

  resolveConflicts_(hypos, hypos_for_mt, features);
}

void OptiQuantAlgorithm::addHypotheses_(Size mono_iso_mt_index, const vector<Size>& candidate_indices, vector<FeatureHypothesis>& hypos, vector<vector<Size> >& hypos_for_mt)
{
  // generate the set of charge states to consider for this monoisotopic mass trace
  set<Int> charges;
  const vector<PeptideIdentification>& ids = kd_data_.feature(mono_iso_mt_index)->getPeptideIdentifications();
  if (use_ids_ && ids.size() && ids[0].getHits().size())
  {
    // consider only charge states supported by the attached IDs
    for (vector<PeptideIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
      const PeptideIdentification& pep_id = *it;
      const vector<PeptideHit>& hits = pep_id.getHits();
      for (vector<PeptideHit>::const_iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
      {
        const PeptideHit& hit = *hit_it;
        charges.insert(hit.getCharge());
      }
    }
  }
  else
  {
    // consider all charge states in the user-specified range
    for (Int c = charge_low_; c <= charge_high_; ++c)
    {
      charges.insert(c);
    }
  }

  // generate hypotheses for all considered charge states
  double mz = kd_data_.mz(mono_iso_mt_index);
  vector<FeatureHypothesis> tmp_hypos;
  for (set<Int>::const_iterator c_it = charges.begin(); c_it != charges.end(); ++c_it)
  {
    Int z = *c_it;
    vector<FeatureHypothesis> hypos_for_charge;
    FeatureHypothesis hypo(z);
    hypo.addMassTrace(make_pair(0, mono_iso_mt_index));
    hypos_for_charge.push_back(hypo);

    for (Size iso_pos = 1; iso_pos < max_nr_traces_; ++iso_pos)
    {
      vector<FeatureHypothesis> hypos_for_iso_pos;

      for (vector<FeatureHypothesis>::iterator h_it = hypos_for_charge.begin();
           h_it != hypos_for_charge.end();
           ++h_it)
      {
        double iso_mz = mz + h_it->getIsoMassDiff() * (double)iso_pos / (double)z;
        pair<double, double> mz_win = Math::getTolWindow(iso_mz, mz_tol_, mz_ppm_);

        for (vector<Size>::const_iterator it = candidate_indices.begin();
             it != candidate_indices.end();
             ++it)
        {
          if (mz_win.first <= kd_data_.mz(*it) && kd_data_.mz(*it) <= mz_win.second)
          {
            // make new hypothesis
            FeatureHypothesis h(*h_it);
            h.addMassTrace(make_pair(iso_pos, *it));

            // adapt isotopic mass difference for this hypo
            if (adaptive_iso_mass_diff_)
            {
              // weighted average for existing and new iso traces
              double new_diff = (h.size() - 1) * h.getIsoMassDiff();
              new_diff += (kd_data_.mz(*it) - kd_data_.mz(h.getMassTraces()[0].second)) / (double)iso_pos * (double)z;
              new_diff /= h.size();
              h.setIsoMassDiff(new_diff);
            }

            // append this compatible mass trace
            hypos_for_iso_pos.push_back(h);
          }
        }
      }

      for (vector<FeatureHypothesis>::iterator it = hypos_for_iso_pos.begin(); it != hypos_for_iso_pos.end(); ++it)
      {
        hypos_for_charge.push_back(*it);
      }
    }

    for (vector<FeatureHypothesis>::iterator it = hypos_for_charge.begin(); it != hypos_for_charge.end(); ++it)
    {
      tmp_hypos.push_back(*it);
    }
  }

  //compile final set of hypotheses
  for (vector<FeatureHypothesis>::iterator it = tmp_hypos.begin(); it != tmp_hypos.end(); ++it)
  {
    const FeatureHypothesis& h = *it;

    // filter out bad hypotheses
    if (h.size() < require_first_n_traces_)
    {
      // too small
      continue;
    }
    for (Size i = 0; i < require_first_n_traces_; ++i)
    {
      if (h.getMassTraces()[i].first != i)
      {
        // one of the first n isotope traces is missing => ignore this hypo
        continue;
      }
    }
    // check correlation with averagine model
    if (computeIntensityScore_(h, true) < min_int_score_thresh_)
    {
      // too bad
      continue;
    }

    // this will be h's index:
    Size hypo_index = hypos.size();
    // ... after this:
    hypos.push_back(h);
    // add this hypothesis to the sets of hypos these mass traces are involved in
    for (vector<pair<Size, Size> >::const_iterator mt_it = h.getMassTraces().begin(); mt_it != h.getMassTraces().end(); ++mt_it)
    {
      hypos_for_mt[mt_it->second].push_back(hypo_index);
    }
  }
}

void OptiQuantAlgorithm::resolveConflicts_(const vector<FeatureHypothesis>& hypos, const vector<vector<Size> >& hypos_for_mt, vector<FeatureHypothesis>& result)
{
  // - find all connected components (HCCs) on hypergraph with nodes = MTs and hyperedges = hypotheses
  // - once a HCC is complete, resolve it and move on to finding the next one
  Size num_nodes = kd_data_.size();
  boost::container::vector<bool> bfs_visited;
  bfs_visited.resize(num_nodes, false);
  queue<Size> bfs_queue;
  Size search_pos = 0;

  // BFS
  UInt progress = 0;
  startProgress(0, hypos.size(), "assembling consensus features");
  while (true)
  {
    bool finished = true;
    for (Size i = search_pos; i < num_nodes; ++i)
    {
      if (!bfs_visited[i])
      {
        bfs_queue.push(i);
        bfs_visited[i] = true;
        finished = false;
        search_pos = i + 1;
        break;
      }
    }
    if (finished) break;

    set<Size> current_hypo_cluster_indices;

    while (!bfs_queue.empty())
    {
      Size i = bfs_queue.front();
      bfs_queue.pop();

      for (vector<Size>::const_iterator it = hypos_for_mt[i].begin(); it != hypos_for_mt[i].end(); ++it)
      {
        current_hypo_cluster_indices.insert(*it);
        const FeatureHypothesis& current_hypo = hypos[*it];

        for (vector<pair<Size, Size> >::const_iterator it2 = current_hypo.getMassTraces().begin(); it2 != current_hypo.getMassTraces().end(); ++it2)
        {
          Size j = it2->second;
          if (!bfs_visited[j])
          {
            bfs_queue.push(j);
            bfs_visited[j] = true;
          }
        }
      }
    }

    // resolve conflicts in this cluster and add final features to result
    resolveHypothesisCluster_(hypos, hypos_for_mt, current_hypo_cluster_indices, result);
    progress += current_hypo_cluster_indices.size();
    setProgress(progress);
  }
  endProgress();
}

void OptiQuantAlgorithm::resolveHypothesisCluster_(const vector<FeatureHypothesis>& hypos, const vector<vector<Size> >& hypos_for_mt, const set<Size>& hypo_cluster_indices, vector<FeatureHypothesis>& result)
{
  Size n = hypo_cluster_indices.size();
  if (n == 0)
  {
    return;
  }
  if (n == 1)
  {
    // solution is trivial
    result.push_back(hypos[*(hypo_cluster_indices.begin())]);
    return;
  }

//  // ---------------------------------
//  double rt_max = 0;
//  double mz_max = 0;
//  double rt_min = 1e10;
//  double mz_min = 1e10;

//  for (set<Size>::const_iterator it = hypo_cluster_indices.begin(); it != hypo_cluster_indices.end(); ++it)
//  {
//    double mz = kd_data_.mz(hypos[*it].getMassTraces()[0].second);
//    double rt = kd_data_.rt(hypos[*it].getMassTraces()[0].second);
//    if (mz < mz_min) mz_min = mz;
//    if (rt < rt_min) rt_min = rt;
//    if (mz > mz_max) mz_max = mz;
//    if (rt > rt_max) rt_max = rt;
//  }
//  cout << "### Resolving hypothesis cluster of size " << n
//       << fixed << setprecision(2)
//       << " (mz = " << mz_min << " - " << mz_max
//       << "; rt = " << rt_min << " - " << rt_max
//       << ") ###" << endl;
//  // ---------------------------------

  vector<Size> hypo_indices(hypo_cluster_indices.begin(), hypo_cluster_indices.end());

  // Solve using CPLEX LP
  IloEnv env;
  try
  {
    IloModel model(env);
    IloNumVarArray x(env);
    //IloRangeArray c(env);

    // compute scores
    vector<double> scores(hypo_indices.size());
    for (Size i = 0; i < hypo_indices.size(); ++i)
    {
      scores[i] = computeScore_(hypos[hypo_indices[i]]);
    }

    // objective function
    IloExpr ilo_var_sum(env);
    for (Size i = 0; i < hypo_indices.size(); ++i)
    {
      IloNumVar ilo_var(env, 0.0, 1.0, IloNumVar::ILOBOOL);
      x.add(ilo_var);
      ilo_var_sum += scores[i] * ilo_var;
    }
    model.add(IloMaximize(env, ilo_var_sum));

//    // constraints
//    for (Size i = 0; i < hypo_indices.size() - 1; ++i)
//    {
//      for (Size j = i + 1; j < hypo_indices.size(); ++j)
//      {
//        bool compatible = true;

//        const FeatureHypothesis& hypo_i = hypos[hypo_indices[i]];
//        const FeatureHypothesis& hypo_j = hypos[hypo_indices[j]];

//        for (Size k = 0; k < hypo_i.size(); ++k)
//        {
//          for (Size l = 0; l < hypo_j.size(); ++l)
//          {
//            if (hypo_i.getMassTraces()[k].second == hypo_j.getMassTraces()[l].second)
//            {
//              // i and j share mass traces => incompatible
//              compatible = false;
//              c.add(x[i] + x[j] <= 1.0);
//              break;
//            }
//          }
//          if (!compatible)
//          {
//            break;
//          }
//        }
//      }
//    }
//    model.add(c);

    // TODO: eventually refactor this (so the following hack and similar stunts are not needed anymore)
    map<Size, Size> hypo_idx_reverse_map;
    for (Size i = 0; i < hypo_indices.size(); ++i)
    {
      hypo_idx_reverse_map[hypo_indices[i]] = i;
    }

    // compile a list of all involved mass trace indices
    set<Size> involved_mass_traces;
    for (Size i = 0; i < hypo_indices.size(); ++i)
    {
      const FeatureHypothesis& h = hypos[hypo_indices[i]];
      const vector<pair<Size, Size> >& h_mts = h.getMassTraces();
      for (Size j = 0; j < h_mts.size(); ++j)
      {
        involved_mass_traces.insert(h_mts[j].second);
      }
    }

    // for each involved mass trace, create a SOS1 constraint so that only one of the hypotheses
    // that this mass trace is involved in can be selected
    IloSOS1Array sos_constraints(env);

    for (set<Size>::const_iterator it = involved_mass_traces.begin(); it != involved_mass_traces.end(); ++it)
    {
      // TODO: speed up by not adding duplicate SOS1s? or does CPLEX eliminate these before solving?

      const vector<Size>& mt_hypos = hypos_for_mt[*it];

      if (mt_hypos.size() < 2)
      {
        // no constraint needed (this mass trace belongs to only a single hypothesis)
        continue;
      }

      // extract scores for hypotheses, add epsilon if scores are identical (HACK)
      // (SOS weights must be unique for CPLEX)
      set<UInt64> unique_weights_set;
      vector<double> unique_weights(mt_hypos.size());
      for (Size i = 0; i < mt_hypos.size(); ++i)
      {
        // get original score for this hypothesis
        Size local_idx = hypo_idx_reverse_map[mt_hypos[i]];
        double weight = scores[local_idx];
        // compute final score, making sure weights are unique:
        // if weight already taken, add small epsilons until unique
        for (;unique_weights_set.count((UInt64)(1e8 * weight)); weight += 1e-6);
        unique_weights_set.insert((UInt64)(1e8 * weight));
        unique_weights[i] = weight;
      }

      // construct SOS1 constraint using the now unique weights
      IloNumVarArray vars(env);
      IloNumArray vals(env);
      for (Size i = 0; i < mt_hypos.size(); ++i)
      {
        Size local_idx = hypo_idx_reverse_map[mt_hypos[i]];
        vars.add(x[local_idx]);
        vals.add(unique_weights[i]);
      }
      sos_constraints.add(IloSOS1(env, vars, vals));
    }
    model.add(sos_constraints);

    // optimize
    IloCplex cplex(model);
    cplex.setParam(IloCplex::Threads, threads_);

    // silence solver
    cplex.setOut(env.getNullStream());

    if (solver_time_limit_ > 0)
    {
      cplex.setParam(IloCplex::TiLim, (double)solver_time_limit_);
    }
    //cplex.exportModel("cplex.lp");

    if (!cplex.solve())
    {
      env.error() << "Failed to optimize LP" << endl;
      throw(-1);
    }
    //env.out() << "Solution status = " << cplex.getStatus() << endl;
    //env.out() << "Solution value  = " << cplex.getObjValue() << endl;

    // read solution values
    IloNumArray vals(env);
    cplex.getValues(vals, x);
    for (Int i = 0; i < vals.getSize(); ++i)
    {
      if (vals[i] > 0.5)
      {
        // add final consensus feature to results
        const FeatureHypothesis& h = hypos[hypo_indices[i]];
        result.push_back(h);
      }
    }
  }
  catch (IloException& e)
  {
    cerr << "Concert exception caught: " << e << endl;
  }
  catch (...)
  {
    cerr << "Unknown exception caught" << endl;
  }

  env.end();
}

double OptiQuantAlgorithm::averagineCorrelation_(const vector<pair<Size, double> >& hypo_int_pairs, const double& mol_weight, bool ignore_missing) const
{
  IsotopeDistribution isodist(max_nr_traces_);
  isodist.estimateFromPeptideWeight(mol_weight);
  vector<pair<Size, double> > averagine_dist = isodist.getContainer();

  vector<double> hypo_ints;
  vector<double> averagine_ints;

  if (ignore_missing)
  {
    hypo_ints.reserve(hypo_int_pairs.size());
    averagine_ints.reserve(hypo_int_pairs.size());
    for (vector<pair<Size, double> >::const_iterator it = hypo_int_pairs.begin(); it != hypo_int_pairs.end(); ++it)
    {
      Size iso_index = it->first;
      double hypo_int = it->second;
      hypo_ints.push_back(hypo_int);
      averagine_ints.push_back(averagine_dist[iso_index].second);
    }
  }
  else
  {
    hypo_ints.resize(max_nr_traces_, 0.0);
    for (vector<pair<Size, double> >::const_iterator it = hypo_int_pairs.begin(); it != hypo_int_pairs.end(); ++it)
    {
      Size iso_index = it->first;
      double hypo_int = it->second;
      hypo_ints[iso_index] = hypo_int;
    }

    averagine_ints.reserve(max_nr_traces_);
    for (vector<pair<Size, double> >::const_iterator it = averagine_dist.begin(); it != averagine_dist.end(); ++it)
    {
      averagine_ints.push_back(it->second);
    }
  }

  return Math::pearsonCorrelationCoefficient(hypo_ints.begin(), hypo_ints.end(),
                                             averagine_ints.begin(), averagine_ints.end());
}

double OptiQuantAlgorithm::computeMZScore_(const FeatureHypothesis& hypo) const
{
  Int z = hypo.getCharge();
  const vector<pair<Size, Size> >& masstraces = hypo.getMassTraces();

  if (masstraces.size() < 2)
  {
    return 0.0;
  }

  double mono_mz = kd_data_.mz(masstraces[0].second);
  vector<double> mz_diffs;

  for (Size i = 1; i < masstraces.size(); ++i)
  {
    Size iso_pos = masstraces[i].first;
    Size mt_index = masstraces[i].second;
    double mz = kd_data_.mz(mt_index);
    double diff = (mz - mono_mz) / (double)iso_pos;
    mz_diffs.push_back(diff);
  }

  double median_diff = Math::median(mz_diffs.begin(), mz_diffs.end());
  double mad = Math::MAD(mz_diffs.begin(), mz_diffs.end(), median_diff);

  // approximate maximum possible diff for scaling to [0,1]
  double max_possible_diff = 2.0 * (mz_ppm_ ? mz_tol_ * mono_mz / 1e6 : mz_tol_);
  double mz_score = 1 - mad / max_possible_diff;
  // rare, but possible (if isotopic mass difference re-estimation is enabled)
  if (mz_score < 0.0) mz_score = 0.0;

  return mz_score;
}

double OptiQuantAlgorithm::computeRTScore_(const FeatureHypothesis& hypo) const
{
  const vector<pair<Size, Size> >& masstraces = hypo.getMassTraces();

  if (masstraces.size() < 2)
  {
    return 0.0;
  }

  vector<double> rts;
  for (Size i = 0; i < masstraces.size(); ++i)
  {
    Size mt_index = masstraces[i].second;
    rts.push_back(kd_data_.rt(mt_index));
  }

  double median_rt = Math::median(rts.begin(), rts.end());
  double mad = Math::MAD(rts.begin(), rts.end(), median_rt);

  double max_possible_deviation = rt_tol_secs_;
  double rt_score = 1 - mad / max_possible_deviation;

  return rt_score;
}

double OptiQuantAlgorithm::computeIntensityScore_(const FeatureHypothesis& hypo, bool consensus_only) const
{
  double z = hypo.getCharge();
  double mz = kd_data_.mz(hypo.getMassTraces()[0].second);
  double mol_weight = z * mz;
  double final_score(0.0);

  if (consensus_only)
  {
    vector<pair<Size, double> > iso_ints;
    for (vector<pair<Size, Size> >::const_iterator it = hypo.getMassTraces().begin(); it != hypo.getMassTraces().end(); ++it)
    {
      Size iso_pos = it->first;
      double intensity = kd_data_.intensity(it->second);
      iso_ints.push_back(make_pair(iso_pos, intensity));
    }
    final_score = (1.0 + averagineCorrelation_(iso_ints, mol_weight, score_int_ignore_missing_));
  }

  else // compute weighted average averagine correlation for all subfeatures from different maps
  {
    // precompute subfeature trace intensities for all maps
    map<Size, vector<pair<Size, double> > > intensities_for_map_idx;
    for (vector<pair<Size, Size> >::const_iterator it = hypo.getMassTraces().begin(); it != hypo.getMassTraces().end(); ++it)
    {
      Size iso_pos = it->first;
      Size mt_index = it->second;

      const ConsensusFeature& mt_cf = (*input_map_)[mt_index];
      const ConsensusFeature::HandleSetType& handles = mt_cf.getFeatures();

      for (ConsensusFeature::HandleSetType::iterator fh_it = handles.begin(); fh_it != handles.end(); ++fh_it)
      {
        Size map_idx = fh_it->getMapIndex();
        if (!intensities_for_map_idx.count(map_idx))
        {
          intensities_for_map_idx[map_idx] = vector<pair<Size, double> >();
        }
        intensities_for_map_idx[map_idx].push_back(make_pair(iso_pos, fh_it->getIntensity()));
      }
    }

    // compute scores for all potential subfeatures
    double summed_score = 0.0;
    Size total_nr_traces = 0;
    for (Size i = 0; i < num_maps_; ++i)
    {
      if (!intensities_for_map_idx.count(i))
      {
        // no traces detected => this subfeature does not contribute to score
        continue;
      }

      // iso-positions and intensities for map i
      const vector<pair<Size, double> >& iso_ints = intensities_for_map_idx[i];

      // number of traces for map i
      Size nr_traces = iso_ints.size();

      if (nr_traces < min_nr_traces_per_map_)
      {
        // too few traces found => this subfeature does not contribute to score
        continue;
      }

      if (require_monoiso_ && iso_ints[0].first != 0)
      {
        // monoisotopic trace missing => this subfeature does not contribute to score
        continue;
      }

      // compute intensity score
      double averagine_score = (1.0 + averagineCorrelation_(iso_ints, mol_weight, score_int_ignore_missing_)) / 2.0;

      // add to combined score
      summed_score += (double)nr_traces * averagine_score;
      total_nr_traces += nr_traces;
    }
    // weighted average averagine similarity score (in [0,1])
    final_score = total_nr_traces ? summed_score / (double)total_nr_traces : 0.0;
  }
  return final_score;
}

double OptiQuantAlgorithm::computeScore_(const FeatureHypothesis& hypo) const
{
  // size score
  double size_score = Math::pown((double)hypo.size() / (double)max_nr_traces_, score_size_exp_);

  // Intensity score
  double int_score = score_int_weight_ * Math::pown(computeIntensityScore_(hypo), score_int_exp_);

  // m/z score
  double mz_score = score_mz_weight_ * Math::pown(computeMZScore_(hypo), score_mz_exp_);

  // RT score
  double rt_score = score_rt_weight_ * Math::pown(computeRTScore_(hypo), score_rt_exp_);

  // combined score
  double combined_score = size_score * (int_score + mz_score + rt_score) / score_denom_;

  // upweight features with identified monoisotopic trace
  double final_score = combined_score;
  if ((*input_map_)[hypo.getMassTraces()[0].second].getPeptideIdentifications().size())
  {
    final_score *= score_id_weight_;
  }

  //
//  cout << "M/Z:\t\t" << kd_data_.mz(hypo.getMassTraces()[0].second) << endl;
//  cout << "RT:\t\t" << kd_data_.rt(hypo.getMassTraces()[0].second) << endl;
//  cout << "Int score:\t" << int_score << endl;
//  cout << "M/Z score:\t" << mz_score << endl;
//  cout << "RT score:\t" << rt_score << endl;
//  cout << "Size score:\t" << size_score << endl;
//  cout << "Final score:\t" << final_score << endl << endl;
  //

  return final_score;
}

void OptiQuantAlgorithm::compileResults_(const vector<FeatureHypothesis>& features, ConsensusMap& output_map)
{
  for (vector<FeatureHypothesis>::const_iterator hypo_it = features.begin(); hypo_it != features.end(); ++hypo_it)
  {
    const FeatureHypothesis& hypo = *hypo_it;

    // precompute feature m/z, RT, intensities
    map<Size, vector<double> > intensities_for_map_idx;
    map<Size, double> mz_for_map_idx;
    map<Size, double> rt_for_map_idx;
    vector<IsoTraceTuple> trace_picker;

    for (vector<pair<Size, Size> >::const_iterator it = hypo.getMassTraces().begin(); it != hypo.getMassTraces().end(); ++it)
    {
      Size iso_pos = it->first;
      Size mt_index = it->second;
      const ConsensusFeature& mt_cf = (*input_map_)[mt_index];
      const ConsensusFeature::HandleSetType& handles = mt_cf.getFeatures();

      IsoTraceTuple itt = {iso_pos, handles.size(), mt_cf.getIntensity()};
      trace_picker.push_back(itt);

      for (ConsensusFeature::HandleSetType::iterator fh_it = handles.begin(); fh_it != handles.end(); ++fh_it)
      {
        Size map_idx = fh_it->getMapIndex();
        if (!intensities_for_map_idx.count(map_idx))
        {
          intensities_for_map_idx[map_idx] = vector<double>(max_nr_traces_, 0.0);
        }
        intensities_for_map_idx[map_idx][iso_pos] = fh_it->getIntensity();

        if (iso_pos == 0)
        {
          // set subfeature m/z and RT based on monoisotopic trace for this map index
          mz_for_map_idx[map_idx] = fh_it->getMZ();
          rt_for_map_idx[map_idx] = fh_it->getRT();
        }
      }
    }

    // determine which consensus traces to use for quantification
    sort(trace_picker.begin(), trace_picker.end());
    vector<Size> quant_iso_positions;
    for (Size i = 0; i < trace_picker.size() && (i < quantify_top_ || quantify_top_ == 0); ++i)
    {
      quant_iso_positions.push_back(trace_picker[i].iso_pos);
    }

    // construct consensus feature
    Size mono_mt_index = hypo.getMassTraces()[0].second;
    const ConsensusFeature& mono_mt_cf = (*input_map_)[mono_mt_index];

    // copy-construct from monoisotopic consensus mass trace to set m/z, RT
    ConsensusFeature final_cf(mono_mt_cf);
    final_cf.setCharge(hypo.getCharge());

    // clear subfeatures
    final_cf.clear();

    // compute subfeatures
    for (Size i = 0; i < num_maps_; ++i)
    {
      if (!intensities_for_map_idx.count(i))
      {
        continue;
      }

      const vector<double>& iso_ints = intensities_for_map_idx[i];
      if (require_monoiso_ && iso_ints[0] == 0.0)
      {
        // monoisotopic trace missing => do not add subfeature for this map
        continue;
      }

      Size nr_detected_traces = 0;
      for (vector<double>::const_iterator it = iso_ints.begin(); it != iso_ints.end(); ++it)
      {
        if (*it > 0.0)
        {
          ++nr_detected_traces;
        }
      }

      if (nr_detected_traces < min_nr_traces_per_map_)
      {
        // too few traces found => do not add subfeature for this map
        continue;
      }

      // compute overall intensity using selected traces
      double summed_int = 0.0;
      for (vector<Size>::const_iterator it = quant_iso_positions.begin(); it != quant_iso_positions.end(); ++it)
      {
        summed_int += iso_ints[*it];
      }

      BaseFeature f;
      if (mz_for_map_idx.count(i) && rt_for_map_idx.count(i))
      {
        // monoisotopic trace found for map i, use its RT and m/z values
        f.setMZ(mz_for_map_idx[i]);
        f.setRT(rt_for_map_idx[i]);
      }
      else
      {
        // monoisotopic trace not found for map i (and require_monoiso_ == false)
        // ==> use values from monoisotopic consensus trace
        f.setMZ(mono_mt_cf.getMZ());
        f.setRT(mono_mt_cf.getRT());
      }
      f.setCharge(hypo.getCharge());
      f.setIntensity(summed_int);

      final_cf.insert(i, f);
    }

    if (final_cf.size())
    {
      // add feature to output map
      output_map.push_back(final_cf);

      // mark its mass traces assembled
      const vector<pair<Size, Size> >& h_mts = hypo.getMassTraces();
      for (Size j = 0; j < h_mts.size(); ++j)
      {
        mt_assembled_[h_mts[j].second] = true;
      }
    }
  }

  // add unassembled mass traces if requested
  if (include_unassembled_traces_)
  {
    for (Size i = 0; i < input_map_->size(); ++i)
    {
      if (!mt_assembled_[i])
      {
        bool id_attached = false;
        Int charge = 0;

        const vector<PeptideIdentification>& pep_ids = (*input_map_)[i].getPeptideIdentifications();
        if (pep_ids.size())
        {
          const vector<PeptideHit>& pep_hits = pep_ids[0].getHits();
          if (pep_hits.size())
          {
            id_attached = true;
            charge = pep_hits[0].getCharge();
          }
        }

        if (include_unidentified_unassembled_traces_ || id_attached)
        {
          ConsensusFeature mt_cf((*input_map_)[i]);
//          // transfer charge from peptide ID if present
//          if (charge != 0)
//          {
//            mt_cf.setCharge(charge);
//            typedef ConsensusFeature::HandleSetType HST;
//            const HST& subfeatures = mt_cf.getFeatures();
//            for (HST::iterator it = subfeatures.begin(); it != subfeatures.end(); ++it)
//            {
//              const_cast<FeatureHandle&>(*it).setCharge(charge);
//            }
//          }
          output_map.push_back(mt_cf);
        }
      }
    }
  }
}

void OptiQuantAlgorithm::outputStatistics_(const ConsensusMap& cmap) const
{
  // some statistics
  map<Size, UInt> num_consfeat_of_size;
  map<Size, UInt> num_consfeat_with_id_of_size;
  map<Size, UInt> num_unassembled_of_size;
  Size total_num_id = 0;
  Size total_num_unassembled = 0;

  for (ConsensusMap::const_iterator cmit = cmap.begin();
       cmit != cmap.end(); ++cmit)
  {
    Size size = cmit->size();
    ++num_consfeat_of_size[size];
    if (cmit->getPeptideIdentifications().size())
    {
      ++num_consfeat_with_id_of_size[size];
      ++total_num_id;
    }
    if (cmit->getCharge() == 0)
    {
      ++num_unassembled_of_size[size];
      ++total_num_unassembled;
    }
  }

  LOG_INFO << "Number of consensus features:" << endl;
  LOG_INFO << "                " << setw(6) << "total " << setw(6) << "ID'ed " << setw(6) << "unassembled" << endl;
  for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin();
       i != num_consfeat_of_size.rend(); ++i)
  {
    Size size = i->first;
    LOG_INFO << "  of size " << setw(2) << size << ": " << setw(6)
             << i->second << " " << setw(6)
             << num_consfeat_with_id_of_size[size] << " "
             << setw(6) << num_unassembled_of_size[size]
             << endl;
  }
  LOG_INFO << "  total:               " << setw(6) << cmap.size() << endl;
  LOG_INFO << "  total identified:    " << setw(6) << total_num_id << endl;
  LOG_INFO << "  total unassembled:   " << setw(6) << total_num_unassembled << endl;
}

void OptiQuantAlgorithm::updateMembers_()
{
  rt_tol_secs_ = (double)(param_.getValue("rt_tol"));
  mz_tol_ = (double)(param_.getValue("mz_tol"));
  mz_ppm_ = (param_.getValue("mz_unit").toString() == "ppm");
  charge_low_ = (Int)(param_.getValue("charge_low"));
  charge_high_ = (Int)(param_.getValue("charge_high"));
  min_int_score_thresh_ = (double)(param_.getValue("min_averagine_score"));
  require_first_n_traces_ = (Size)(param_.getValue("require_first_n_traces"));
  min_nr_traces_per_map_ = (Size)(param_.getValue("min_nr_traces_per_map"));
  max_nr_traces_ = (Size)(param_.getValue("max_nr_traces"));
  use_ids_ = (param_.getValue("use_ids").toString() == "true");
  require_monoiso_ = (param_.getValue("require_monoiso").toString() == "true");
  solver_time_limit_ = param_.getValue("solver_time_limit");
  include_unassembled_traces_ = (param_.getValue("keep_unassembled_traces").toString() != "none");
  include_unidentified_unassembled_traces_ = (param_.getValue("keep_unassembled_traces").toString() == "all");
  score_id_weight_ = (double)(param_.getValue("score:id_weight"));
  score_size_exp_ = (UInt)(param_.getValue("score:size_exp"));
  score_int_exp_ = (UInt)(param_.getValue("score:int_exp"));
  score_int_weight_ = (double)(param_.getValue("score:int_weight"));
  score_mz_exp_ = (UInt)(param_.getValue("score:mz_exp"));
  score_mz_weight_ = (double)(param_.getValue("score:mz_weight"));
  score_rt_exp_ = (UInt)(param_.getValue("score:rt_exp"));
  score_rt_weight_ = (double)(param_.getValue("score:rt_weight"));
  score_denom_ = score_int_weight_ + score_mz_weight_ + score_rt_weight_;
  score_int_ignore_missing_ = (param_.getValue("score:int_ignore_missing").toString() == "true");
  quantify_top_ = (Size)(param_.getValue("quantify_top"));
  adaptive_iso_mass_diff_ = (param_.getValue("adaptive_iso_mass_diff").toString() == "true");
}

}
