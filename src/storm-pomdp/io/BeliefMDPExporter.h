//
// Created by spook on 12.01.24.
//

#pragma once

#include "storm/io/GEXFExporter.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelChecker.h"
#include "storm-pomdp/transformer/BoundUnfolder.h"

namespace storm {
namespace exporter {

template<typename ValueType, typename BeliefType> // TODO are these sufficient / right?
class BeliefMDPExporter: public GEXFExporter<ValueType, BeliefType> {
   public:
    BeliefMDPExporter() = default;

     void createGEXFOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo, std::string filename);

   private:

    void outputSummary(std::vector<std::string> numEpochs, std::vector<std::string> maxSingleStateEpochs, bool under, std::string filename);

    void determineOgColors(std::vector<std::vector<uint64_t>> &stateColors);

    void determineOgColorsExtended(std::vector<std::vector<uint64_t>> &stateColors);

    void determineNumberOfEpochs(std::vector<std::string> &numbersOfEpochs, std::vector<std::string> &maxSingleStateEpochNumbers, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

    void determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

    std::string summaryFile = "/home/spook/Documents/bunfbench_results/summary.txt";

private:
    std::tuple<uint64_t, uint64_t, uint64_t> adaptedEuclid(uint64_t i, uint64_t x, uint64_t y);
};
} // storm
} // exporter

