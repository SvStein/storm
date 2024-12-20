//
// Created by spook on 11.01.24.
//

#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/NondeterministicModel.h"
#include "storm/models/sparse/DeterministicModel.h"
#include "storm/models/sparse/Pomdp.h"
#include "storm/models/sparse/StochasticTwoPlayerGame.h"
#include "storm/models/sparse/Smg.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Dtmc.h"


namespace storm {
    namespace exporter {

        template<typename ValueType>
        class GEXFExporter {
        public:
            enum GEXFAttributeType { GEXF_string, GEXF_integer, GEXF_long, GEXF_float, GEXF_double, GEXF_boolean, GEXF_short, GEXF_byte, GEXF_date, GEXF_anyURI };

            GEXFExporter() = default;

            void exportGEXFToStream(std::shared_ptr<storm::models::sparse::Model<ValueType>> model, std::ostream& outStream, std::optional<std::vector<std::vector<uint64_t>>> colorVec, std::optional<std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>> additionalAttributesMap, bool allowDefaultColoring, bool convertRatesToProbs);

            std::vector<std::vector<uint64_t>> createPomdpObservationColoring(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> pomdp);

            std::vector<std::vector<uint64_t>> createSmgPlayerColoring(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg, const std::map<storage::PlayerIndex, uint64_t>& playerInfo);

            std::vector<std::vector<uint64_t>> createMaStateTypeColoring(std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> ma);

            std::vector<std::vector<uint64_t>> createBasicColoring(uint64_t numberOfColors);



        protected:
            uint64_t encodeColor(uint64_t r, uint64_t g, uint64_t b);

            std::string attributeTypeToString(GEXFAttributeType attributeType);

            std::tuple<uint64_t, uint64_t, uint64_t> adaptedEuclid(uint64_t i, uint64_t x, uint64_t y);

            void exportNonDetGEXFToStream(std::shared_ptr<storm::models::sparse::NondeterministicModel<ValueType>> model, std::ostream& outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes, storm::storage::SparseMatrix<ValueType> transMatrix);

            void exportDetGEXFToStream(std::shared_ptr<storm::models::sparse::DeterministicModel<ValueType>> model, std::ostream& outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes, storm::storage::SparseMatrix<ValueType> transMatrix);

            void prepPomdp(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> pomdp, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring);

            void prepMdp(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring);

            void prepDtmc(std::shared_ptr<storm::models::sparse::Dtmc<ValueType>> dtmc, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring);

            void prepCtmc(std::shared_ptr<storm::models::sparse::Ctmc<ValueType>> ctmc, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring, bool convertRatesToProbs);

            void prepMa(std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> ma, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring, bool convertRatesToProbs);

            void prepSmg(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg, std::vector<std::vector<uint64_t>>& colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>& additionalAttributes, bool allowDefaultColoring);

            std::map<storage::PlayerIndex,uint64_t> createPlayersToColorIdxMap(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg);

            std::map<storage::PlayerIndex,std::string> createPlayersToNameMap(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg, const std::map<storage::PlayerIndex, uint64_t>& playerToColorIdx);

            std::vector<std::string> collectStateLabelings(std::shared_ptr<storm::models::sparse::Model<ValueType>> model);

        };

    } // storm
} // exporter

