//
// Created by spook on 11.01.24.
//

#ifndef STORM_GEXFEXPORTER_H
#define STORM_GEXFEXPORTER_H

#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/Pomdp.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelChecker.h"
#include "storm-pomdp/transformer/BoundUnfolder.h"

namespace storm {
    namespace exporter {

        template<typename ValueType, typename BeliefType> // TODO are these sufficient / right?
        class GEXFExporter {
        public:
            enum GEXFAttributeType { GEXF_string, GEXF_integer, GEXF_long, GEXF_float, GEXF_double, GEXF_boolean, GEXF_short, GEXF_byte, GEXF_date, GEXF_anyURI };

            GEXFExporter() = default;

            void createGEXFOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo, std::string filename);
        private:
            uint64_t encodeColor(uint64_t r, uint64_t g, uint64_t b);

            std::string attributeTypeToString(GEXFAttributeType attributeType);

            void exportGEXFToStream(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp, std::ostream& outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes = std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>());

            void determineOgColors(std::vector<std::vector<uint64_t>> &stateColors);

            void determineNumberOfEpochs(std::vector<std::string> &numbersOfEpochs, std::vector<std::string> &maxSingleStateEpochNumbers, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

            void determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox);

        };

    } // storm
} // exporter

#endif //STORM_GEXFEXPORTER_H
