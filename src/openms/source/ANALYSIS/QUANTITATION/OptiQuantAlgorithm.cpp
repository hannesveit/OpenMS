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

OptiQuantAlgorithm::OptiQuantAlgorithm(const ConsensusMap& input_map) :
  DefaultParamHandler("OptiQuantAlgorithm"),
  ProgressLogger()
{
  input_map_ = &input_map;
  num_maps_ = input_map.getFileDescriptions().size();
  kd_data_.addFeatures(0, input_map, true);
  mt_assembled_ = vector<Int>(input_map.size(), false);

  defaults_.setValue("mz_tol", 10.0, "m/z tolerance (in Da or ppm)");
  defaults_.setMinFloat("mz_tol", 0.0);

  defaults_.setValue("mz_unit", "ppm", "unit of m/z tolerance");
  defaults_.setValidStrings("mz_unit", ListUtils::create<String>("ppm,Da"));

  defaults_.setValue("rt_tol", 20.0, "RT tolerance (in seconds)");
  defaults_.setMinFloat("rt_tol", 0.0);

  defaults_.setValue("charge_low", 2, "Lowest charge state to consider");
  defaults_.setMinInt("charge_low", 1);

  defaults_.setValue("charge_high", 5, "Highest charge state to consider");
  defaults_.setMinInt("charge_high", 1);

  defaults_.setValue("require_first_n_traces", 3, "Do not consider consensus feature hypotheses in which any of the first n isotope traces are missing across all maps (including the monoisotopic trace)");
  defaults_.setMinInt("require_first_n_traces", 1);

  defaults_.setValue("min_nr_traces_per_map", 2, "Ignore subfeatures with less than this many detected mass traces");
  defaults_.setMinInt("min_nr_traces_per_map", 1);

  defaults_.setValue("max_nr_traces", 6, "Search only for the first max_nr_traces isotope traces");
  defaults_.setMinInt("max_nr_traces", 1);

  defaults_.setValue("require_monoiso", "true", "Include subfeature for map i only if the monoisotopic trace of this feature has been detected in map i");
  defaults_.setValidStrings("require_monoiso", ListUtils::create<String>("true,false"));

  defaults_.setValue("use_ids", "true", "If a mass trace has identifications attached, generate only hypotheses for the charge state(s) found in these peptide IDs when generating hypotheses for this (monoisotopic) mass trace.");
  defaults_.setValidStrings("use_ids", ListUtils::create<String>("true,false"));

  defaults_.setValue("keep_unassembled_traces", "identified", "Include unassembled traces in the results? When set to 'identified', keep only those traces annotated with at least one peptide identification");
  defaults_.setValidStrings("keep_unassembled_traces", ListUtils::create<String>("all,identified,none"));

  defaults_.setValue("solver_time_limit", -1, "CPLEX time limit (in seconds) for solving one cluster of contiguous hypotheses. No time limit when set to -1.");

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
}

void OptiQuantAlgorithm::assembleFeatures_(vector<FeatureHypothesis>& features)
{
  // somewhat hacky (see below)
  double epsilon = 10e-7;

  // compute m/z window size
  double mz_win_height = (1.000857*(double)max_nr_traces_ + 0.001091) / (double)charge_low_; // TODO: confirm values, make variable

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
  vector<FeatureHypothesis> new_hypos;
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

      // TODO: speed up using tolerance window query instead of looping over all candidates?
      double iso_mz = mz + (1.000857*(double)iso_pos + 0.001091) / (double)z; // TODO: confirm values, make variable
      pair<double, double> mz_win = Math::getTolWindow(iso_mz, mz_tol_, mz_ppm_);
      for (vector<Size>::const_iterator it = candidate_indices.begin();
                                        it != candidate_indices.end();
                                        ++it)
      {
        if (mz_win.first <= kd_data_.mz(*it) && kd_data_.mz(*it) <= mz_win.second)
        {
          for (vector<FeatureHypothesis>::iterator h_it = hypos_for_charge.begin();
                                                   h_it != hypos_for_charge.end();
                                                   ++h_it)
          {
            FeatureHypothesis h(*h_it); // make new hypothesis
            h.addMassTrace(make_pair(iso_pos, *it)); // append this compatible mass trace
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
      new_hypos.push_back(*it);
    }
  }

  //compile final set of hypotheses
  for (vector<FeatureHypothesis>::iterator it = new_hypos.begin(); it != new_hypos.end(); ++it)
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
      set<double> unique_weights_set;
      vector<double> unique_weights(mt_hypos.size());
      for (Size i = 0; i < mt_hypos.size(); ++i)
      {
        // get original score for this hypothesis
        Size local_idx = hypo_idx_reverse_map[mt_hypos[i]];
        double weight = scores[local_idx];
        // compute final score, making sure weights are unique:
        // if weight already taken, add small epsilons until unique
        for (;unique_weights_set.count(weight); weight += 1e-6);
        unique_weights_set.insert(weight);
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

double OptiQuantAlgorithm::computeIntensityScore_(const std::vector<double>& hypo_ints, const double& mol_weight) const
{
  IsotopeDistribution isodist(hypo_ints.size());
  isodist.estimateFromPeptideWeight(mol_weight);

  std::vector<std::pair<Size, double> > averagine_dist = isodist.getContainer();

  double max_int(0.0);
  double theo_max_int(0.0);

  for (Size i = 0; i < hypo_ints.size(); ++i)
  {
    if (hypo_ints[i] > max_int)
    {
      max_int = hypo_ints[i];
    }

    if (averagine_dist[i].second > theo_max_int)
    {
      theo_max_int = averagine_dist[i].second;
    }
  }

  std::vector<double> averagine_ratios;
  std::vector<double> detected_ratios;

  for (Size i = 0; i < hypo_ints.size(); ++i)
  {
    averagine_ratios.push_back(averagine_dist[i].second / theo_max_int);
    detected_ratios.push_back(hypo_ints[i] / max_int);
  }

  return computeCosineSim_(averagine_ratios, detected_ratios);
}

double OptiQuantAlgorithm::computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const
{
  if (x.size() != y.size())
  {
    return 0.0;
  }

  double mixed_sum(0.0);
  double x_squared_sum(0.0);
  double y_squared_sum(0.0);


  for (Size i = 0; i < x.size(); ++i)
  {
    mixed_sum += x[i] * y[i];
    x_squared_sum += x[i] * x[i];
    y_squared_sum += y[i] * y[i];
  }

  double denom(std::sqrt(x_squared_sum) * std::sqrt(y_squared_sum));

  return (denom > 0.0) ? mixed_sum / denom : 0.0;
}

double OptiQuantAlgorithm::computeScore_(const FeatureHypothesis& hypo) const
{
  double z = hypo.getCharge();
  double mz = kd_data_.mz(hypo.getMassTraces()[0].second);
  double mol_weight = z * mz;
  map<Size, vector<double> > intensities_for_map_idx;

  // precompute feature intensities
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
        intensities_for_map_idx[map_idx] = vector<double>(max_nr_traces_, 0.0);
      }
      intensities_for_map_idx[map_idx][iso_pos] = fh_it->getIntensity();
    }
  }

  double summed_score = 0.0;
  // compute scores for all potential subfeatures
  for (Size i = 0; i < num_maps_; ++i)
  {
    if (!intensities_for_map_idx.count(i))
    {
      // no traces detected => this subfeature does not contribute to score
      continue;
    }

    const vector<double>& iso_ints = intensities_for_map_idx[i];
    if (require_monoiso_ && iso_ints[0] == 0.0)
    {
      // monoisotopic trace missing => this subfeature does not contribute to score
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
      // too few traces found => this subfeature does not contribute to score
      continue;
    }

    // compute averagine score
    double averagine_score = computeIntensityScore_(iso_ints, mol_weight);

    // weight by number of detected traces in map i
    averagine_score *= (double)nr_detected_traces;

    // add to combined score
    summed_score += averagine_score;
  }

  return summed_score;
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

      double summed_int = accumulate(iso_ints.begin(), iso_ints.end(), 0.0);

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
          mt_cf.setCharge(charge);
          output_map.push_back(mt_cf);
        }
      }
    }
  }
}

void OptiQuantAlgorithm::updateMembers_()
{
  rt_tol_secs_ = (double)(param_.getValue("rt_tol"));
  mz_tol_ = (double)(param_.getValue("mz_tol"));
  mz_ppm_ = (param_.getValue("mz_unit").toString() == "ppm");
  charge_low_ = (Int)(param_.getValue("charge_low"));
  charge_high_ = (Int)(param_.getValue("charge_high"));
  require_first_n_traces_ = (Size)(param_.getValue("require_first_n_traces"));
  min_nr_traces_per_map_ = (Size)(param_.getValue("min_nr_traces_per_map"));
  max_nr_traces_ = (Size)(param_.getValue("max_nr_traces"));
  use_ids_ = (param_.getValue("use_ids").toString() == "true");
  require_monoiso_ = (param_.getValue("require_monoiso").toString() == "true");
  solver_time_limit_ = param_.getValue("solver_time_limit");
  include_unassembled_traces_ = (param_.getValue("keep_unassembled_traces").toString() != "none");
  include_unidentified_unassembled_traces_ = (param_.getValue("keep_unassembled_traces").toString() == "all");
}

}
