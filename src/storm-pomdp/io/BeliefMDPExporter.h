//
// Created by spook on 12.01.24.
//

#ifndef STORM_BELIEFMDPEXPORTER_H
#define STORM_BELIEFMDPEXPORTER_H

#include "storm/io/GEXFExporter.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelChecker.h"
#include "storm-pomdp/transformer/BoundUnfolder.h"

namespace storm {
namespace exporter {

template<typename ValueType, typename BeliefType> // TODO are these sufficient / right?
class BeliefMDPExporter: public GEXFExporter {
   public:
    BeliefMDPExporter() = default;

     void createGEXFOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo, std::string filename);

   private:

    void determineOgColors(std::vector<std::vector<uint64_t>> &stateColors);

    void determineNumberOfEpochs(std::vector<std::string> &numbersOfEpochs, std::vector<std::string> &maxSingleStateEpochNumbers, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

    void determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

};

} // storm
} // exporter

#endif //STORM_BELIEFMDPEXPORTER_H
