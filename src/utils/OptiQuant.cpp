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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/OptiQuantAlgorithm.h>

#include <fstream>
#include <iostream>

using namespace OpenMS;
using namespace std;

//
//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_OptiQuant OptiQuant

 @brief Label-free quantification on linked mass traces using ILP.

 Label-free quantification on linked mass traces using ILP.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_OptiQuant.cli

 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPOptiQuant : public TOPPBase, public ProgressLogger
{
public:

  TOPPOptiQuant() :
    TOPPBase("OptiQuant",
             "Label-free quantification on linked mass traces using ILP.", false),
    ProgressLogger()
  {
    setLogType(CMD);
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Linked mass traces");
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out", "<file>", "", "Linked assembled features");
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));

    addEmptyLine_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    return OptiQuantAlgorithm().getParameters();
  }

  ExitCodes main_(int, const char **)
  {
    // read input
    startProgress(0, 1, "reading linked mass traces");
    String in_filename = getStringOption_("in");
    ConsensusMap input_map;
    ConsensusXMLFile().load(in_filename, input_map);
    endProgress();

    // initialize algorithm
    Param p(getParam_());
    Param p_optiquant_algo = p.copy("algorithm:", true);
    OptiQuantAlgorithm oq_algo(input_map);
    oq_algo.setParameters(p_optiquant_algo);

    // prepare output map
    ConsensusMap output_map;
    output_map.setFileDescriptions(input_map.getFileDescriptions());
    output_map.ensureUniqueId();

    // run algorithm
    oq_algo.run(output_map);

    // store results
    startProgress(0, 1, "storing results");
    String out_filename = getStringOption_("out");
    ConsensusXMLFile().store(out_filename, output_map);
    endProgress();

    return EXECUTION_OK;
  }

};

int main(int argc, const char **argv)
{
  TOPPOptiQuant tool;
  return tool.main(argc, argv);
}

/// @endcond
