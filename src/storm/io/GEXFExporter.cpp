//
// Created by spook on 11.01.24.
//

#include <storm/exceptions/NotSupportedException.h>
#include "GEXFExporter.h"
#include "storm/io/file.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/NondeterministicModel.h"
#include "storm/models/sparse/DeterministicModel.h"

namespace storm {
    namespace exporter {

        template<typename ValueType>
        uint64_t GEXFExporter<ValueType>::encodeColor(uint64_t r, uint64_t g, uint64_t b) {
            STORM_LOG_ASSERT(r < 256 && g < 256 && b < 256, "Invalid color values. Need to be in [0, 255].");
            return (r*1000000 + g*1000 + b);
        }

        template<typename ValueType>
        std::string GEXFExporter<ValueType>::attributeTypeToString(GEXFAttributeType attributeType) {
            switch (attributeType) {
                case GEXF_string: return "string";
                case GEXF_integer: return "integer";
                case GEXF_long: return "long";
                case GEXF_float: return "float";
                case GEXF_double: return "double";
                case GEXF_boolean: return "boolean";
                case GEXF_short: return "short";
                case GEXF_byte: return "byte";
                case GEXF_date: return "date";
                case GEXF_anyURI: return "anyURI";
            }
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::exportGEXFToStream(
                std::shared_ptr<storm::models::sparse::Model<ValueType>> model, std::ostream &outStream,
                std::optional<std::vector<std::vector<uint64_t>>> colorVec,
                std::optional<std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>> additionalAttributesMap, bool allowDefaultColoring, bool convertRatesToProbs) {

            // for nonexisting colors + additionalAttributes, use empty vec/map bc then we don't have to deal with either constantly using .value() of the optionals but also
            // don't have to have eg colors given to be able to give addAttr as we would if we gave the params default values
            auto colors = colorVec? colorVec.value() : std::vector<std::vector<uint64_t>>();
            std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes = additionalAttributesMap? additionalAttributesMap.value() : std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>();

            if (!colors.empty() && colors.size() != model->getNumberOfStates()) {
                STORM_LOG_WARN("Color vector has weird size (not zero indicating no coloring, but also not equal to the number of states) and will be reset now.");
                colors = std::vector<std::vector<uint64_t>>();
            }

            // All model types may have labels and they seem reasonable to include as node attributes
            if (!additionalAttributes.contains("Labels")) {
                additionalAttributes["Labels"] = std::make_pair(GEXFAttributeType::GEXF_string, collectStateLabelings(model));
            } else {
                STORM_LOG_WARN("Attribute \"Labels\" is already given and will overwrite the default attribute that is usually generated.");
            }


            // TODO check correct size etc of additionalAttributes and/ or colors or maybe in actual export function?

            // distinguish between model types, call their respective prep functions (introduce coloring and some attributes if required / allowed), then export
            // for models with continuous time models, distinguish between output in "Rate Matrix" format vs. "Probability Matrix + State Rates as attributes" format
            switch (model->getType()) {
                case models::ModelType::Pomdp:
                    prepPomdp(model->template as<storm::models::sparse::Pomdp<ValueType>>(), colors, additionalAttributes, allowDefaultColoring);
                    exportNonDetGEXFToStream(model->template as<storm::models::sparse::NondeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    break;
                case models::ModelType::Mdp:
                    prepMdp(model->template as<storm::models::sparse::Mdp<ValueType>>(), colors, additionalAttributes, allowDefaultColoring);
                    exportNonDetGEXFToStream(model->template as<storm::models::sparse::NondeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    break;
                case models::ModelType::MarkovAutomaton:
                    prepMa(model->template as<storm::models::sparse::MarkovAutomaton<ValueType>>(), colors, additionalAttributes, allowDefaultColoring, convertRatesToProbs);
                    if (convertRatesToProbs) {
                        // We are assuming here (due to the MA constructors (and only those) using the turnRatesToProbabilities function whenever the components give rates) that an MA by default uses a prob matrix
                        exportNonDetGEXFToStream(model->template as<storm::models::sparse::NondeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    } else {
                        exportNonDetGEXFToStream(model->template as<storm::models::sparse::NondeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->template as<storm::models::sparse::MarkovAutomaton<ValueType>>()->computeRateMatrix());
                    }
                    break;
                case models::ModelType::S2pg:
                    // TODO special treatment for this one needed, where do i even start
                    break;
                case models::ModelType::Smg:
                    prepSmg(model->template as<storm::models::sparse::Smg<ValueType>>(), colors, additionalAttributes, allowDefaultColoring);
                    exportNonDetGEXFToStream(model->template as<storm::models::sparse::NondeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    break;
                case models::ModelType::Dtmc:
                    prepDtmc(model->template as<storm::models::sparse::Dtmc<ValueType>>(), colors, additionalAttributes, allowDefaultColoring);
                    exportDetGEXFToStream(model->template as<storm::models::sparse::DeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    break;
                case models::ModelType::Ctmc:
                    prepCtmc(model->template as<storm::models::sparse::Ctmc<ValueType>>(), colors, additionalAttributes, allowDefaultColoring, convertRatesToProbs);
                    if (convertRatesToProbs) {
                        exportDetGEXFToStream(model->template as<storm::models::sparse::DeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->template as<storm::models::sparse::Ctmc<ValueType>>()->computeProbabilityMatrix());
                    } else {
                        exportDetGEXFToStream(model->template as<storm::models::sparse::DeterministicModel<ValueType>>(), outStream, colors, additionalAttributes, model->getTransitionMatrix());
                    }
                    break;
                default:
                    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "GEXF export not supported for this model type"); // for possible future model types
            }
        }

        template<typename ValueType>
        std::vector<std::string>
        GEXFExporter<ValueType>::collectStateLabelings(std::shared_ptr<storm::models::sparse::Model<ValueType>> model) {
            auto labels = std::vector<std::string>(model->getNumberOfStates());
            for (uint64_t state = 0; state < model->getNumberOfStates(); state++) {
                auto labelsOfState = model->getStateLabeling().getLabelsOfState(state);
                std::string labelString;
                bool first = true;
                for (const auto& l: labelsOfState) {
                    if (first) {
                        labelString = l;
                        first = false;
                    } else {
                        labelString += ", " + l;
                    }
                }
                labels[state] = labelString;
            }
            return labels;
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepCtmc(std::shared_ptr<storm::models::sparse::Ctmc<ValueType>> ctmc,
                                               std::vector<std::vector<uint64_t>> &colors,
                                               std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes,
                                               bool allowDefaultColoring, bool convertRatesToProbs) {
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                colors = createBasicColoring(ctmc->getNumberOfStates());
            }

            if (convertRatesToProbs) {
                if (!additionalAttributes.contains("ExitRate")) {
                    auto stateRateVec = ctmc->getExitRateVector();
                    auto rates = std::vector<std::string>(ctmc->getNumberOfStates());
                    for (uint64_t state = 0; state < ctmc->getNumberOfStates(); state++) {
                        rates[state] = storm::utility::to_string(stateRateVec[state]);
                    }
                    additionalAttributes["ExitRate"] = std::make_pair(GEXFAttributeType::GEXF_string, std::move(
                            rates)); // using string here bc I don't know if rationals would count as valid double but also don't wanna convert unnecessarily
                } else {
                    STORM_LOG_WARN(
                            "Attribute \"ExitRate\" is already given and will overwrite the default attribute that is usually generated.");
                }
            }

        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepPomdp(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> pomdp,
                                                std::vector<std::vector<uint64_t>> &colors,
                                                std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes, bool allowDefaultColoring) {
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                colors = createPomdpObservationColoring(pomdp);
            }
            if (!additionalAttributes.contains("Observation")) {
                auto observations = std::vector<std::string>(pomdp->getNumberOfStates());
                for (uint64_t state = 0; state < pomdp->getNumberOfStates(); state++) {
                    observations[state] = pomdp->getObservation(state);
                }
                additionalAttributes["Observation"] = std::make_pair(GEXFAttributeType::GEXF_integer, std::move(observations));
            } else {
                STORM_LOG_WARN("Attribute \"Observation\" is already given and will overwrite the default attribute that is usually generated.");
            }
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepMdp(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp,
                                              std::vector<std::vector<uint64_t>> &colors,
                                              std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes,
                                              bool allowDefaultColoring) {
            // This doesn't do much but still putting it in its own function so every model type has a prep function and is thus easily expandable
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                colors = createBasicColoring(mdp->getNumberOfStates());
            }
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepDtmc(std::shared_ptr<storm::models::sparse::Dtmc<ValueType>> dtmc,
                                               std::vector<std::vector<uint64_t>> &colors,
                                               std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes,
                                               bool allowDefaultColoring) {
            // This doesn't do much but still putting it in its own function so every model type has a prep function and is thus easily expandable
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                colors = createBasicColoring(dtmc->getNumberOfStates());
            }
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepMa(std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> ma,
                                             std::vector<std::vector<uint64_t>> &colors,
                                             std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes, bool allowDefaultColoring, bool convertRatesToProbs) {
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                colors = createMaStateTypeColoring(ma);
            }

            if (!additionalAttributes.contains("StateType")) {
                auto stateTypes = std::vector<std::string>(ma->getNumberOfStates());
                for (uint64_t state = 0; state < ma->getNumberOfStates(); state++) {
                    if (ma->isHybridState(state)) {
                        stateTypes[state] = "Hybrid";
                    } else if (ma->isProbabilisticState(state)) {
                        stateTypes[state] = "Probabilistic";
                    } else {
                        assert(ma->isMarkovianState(state));
                        stateTypes[state] = "Markovian";
                    }
                }
                additionalAttributes["StateType"] = std::make_pair(GEXFAttributeType::GEXF_string, std::move(stateTypes));
            } else {
                STORM_LOG_WARN("Attribute \"StateType\" is already given and will overwrite the default attribute that is usually generated.");
            }

            if (convertRatesToProbs) {
                if (!additionalAttributes.contains("ExitRate")) {
                    auto stateRateVec = ma->getExitRates();
                    auto rates = std::vector<std::string>(ma->getNumberOfStates());
                    for (uint64_t state = 0; state < ma->getNumberOfStates(); state++) {
                        if (!ma->isProbabilisticState(state)) {
                            rates[state] = storm::utility::to_string(stateRateVec[state]);
                        } else {
                            rates[state] = "-";
                        }
                    }
                    additionalAttributes["ExitRate"] = std::make_pair(GEXFAttributeType::GEXF_string, std::move(
                            rates)); // using string here bc I don't know if rationals would count as valid double but also don't wanna convert unnecessarily
                } else {
                    STORM_LOG_WARN(
                            "Attribute \"ExitRate\" is already given and will overwrite the default attribute that is usually generated.");
                }
            }
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::prepSmg(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg,
                                              std::vector<std::vector<uint64_t>> &colors,
                                              std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> &additionalAttributes,
                                              bool allowDefaultColoring) {
            std::optional<std::map<storage::PlayerIndex, uint64_t>> playerToColorIdx;
            if (colors.empty() && allowDefaultColoring) {
                STORM_LOG_INFO("Applying default coloring.");
                playerToColorIdx = createPlayersToColorIdxMap(smg);
                colors = createSmgPlayerColoring(smg, playerToColorIdx.value());
            }

            if (!additionalAttributes.contains("PlayerIndex")) {
                auto playerIndices = std::vector<std::string>(smg->getNumberOfStates());
                for (uint64_t state = 0; state < smg->getNumberOfStates(); state++) {
                    playerIndices[state] = std::to_string(smg->getPlayerOfState(state));
                }
                additionalAttributes["PlayerIndex"] = std::make_pair(GEXFAttributeType::GEXF_integer, std::move(playerIndices));
            } else {
                STORM_LOG_WARN("Attribute \"PlayerIndex\" is already given and will overwrite the default attribute that is usually generated.");
            }
            if (!additionalAttributes.contains("PlayerName") && !smg->getPlayerNameToIndexMap().empty()) {
                if (!playerToColorIdx) {
                    playerToColorIdx = createPlayersToColorIdxMap(smg);
                }
                auto playerToName = createPlayersToNameMap(smg, playerToColorIdx.value());
                auto playerNames = std::vector<std::string>(smg->getNumberOfStates());
                for (uint64_t state = 0; state < smg->getNumberOfStates(); state++) {
                    playerNames[state] = playerToName[smg->getPlayerOfState(state)];
                }
                additionalAttributes["PlayerName"] = std::make_pair(GEXFAttributeType::GEXF_string, std::move(playerNames));
            } else {
                STORM_LOG_WARN("Attribute \"PlayerName\" is already given and will overwrite the default attribute that is usually generated.");
            }


        }

        template<typename ValueType>
        std::map<storage::PlayerIndex, std::string>
        GEXFExporter<ValueType>::createPlayersToNameMap(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg, const std::map<storage::PlayerIndex, uint64_t>& playerToColorIdx) {
            // we use the color index map here because that one for sure has all occurring players as keys, incl. invalid index or possibly just players without names
            auto playerToName = std::map<storage::PlayerIndex, std::string>();
            for (const auto& namePlayerPair : smg->getPlayerNameToIndexMap()) {
                playerToName[namePlayerPair.second] = namePlayerPair.first + " (" + std::to_string(namePlayerPair.second) + ")";
            }
            for (auto playerColorPair : playerToColorIdx) {
                if (!playerToName.contains(playerColorPair.first)) {
                    if (playerColorPair.first == storage::INVALID_PLAYER_INDEX) {
                        playerToName[playerColorPair.first] = "Invalid Player Index (" + std::to_string(playerColorPair.first) + ")";
                    } else {
                        playerToName[playerColorPair.first] = "Nameless (" + std::to_string(playerColorPair.first) + ")";
                    }
                }
            }
            STORM_LOG_ASSERT(playerToName.size() == playerToColorIdx.size(), "Uh oh!");
            return playerToName;
        }

        template<typename ValueType>
        std::map<storage::PlayerIndex, uint64_t> GEXFExporter<ValueType>::createPlayersToColorIdxMap(
                std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg) {
            // Collect players
            auto playerToColorIdx = std::map<storage::PlayerIndex,uint64_t>();
            std::vector<storage::PlayerIndex> players;
            std::copy(smg->getStatePlayerIndications().begin(), smg->getStatePlayerIndications().end(), std::back_inserter(players));
            std::sort(players.begin(), players.end());
            auto last = std::unique(players.begin(), players.end());
            players.erase(last, players.end());

            for (uint64_t i = 0; i < players.size(); i++) {
                playerToColorIdx[players[i]] = i;
            }
            return playerToColorIdx;
        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::exportNonDetGEXFToStream(std::shared_ptr<storm::models::sparse::NondeterministicModel<ValueType>> model, std::ostream &outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes, storm::storage::SparseMatrix<ValueType> transMatrix) {
            STORM_LOG_ASSERT(colors.empty() || colors.size() == model->getNumberOfStates(), "Color vector has weird size (not zero indicating no coloring, but also not equal to the number of states)!");
            bool useColors = model->getNumberOfStates() == colors.size();

            // Preamble
            outStream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                         "<gexf xmlns=\"http://gexf.net/1.3\"\n"
                         "      xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                         "      xmlns:viz=\"http://gexf.net/1.3/viz\"\n"
                         "      xsi:schemaLocation=\"http://gexf.net/1.3\n"
                         "                          http://gexf.net/1.3/gexf.xsd\"\n"
                         "      version=\"1.3\">\n"
                         "  <meta>\n"
                         "    <creator>Storm</creator>\n"
                         "    <description>Visualization of a non-deterministic model</description>\n"
                         "  </meta>\n"
                         "  <graph defaultedgetype=\"directed\">\n"
                         "    <attributes class=\"node\">\n"
                         "      <attribute id=\"0\" title=\"init\" type=\"boolean\"/>\n"
                         "      <attribute id=\"1\" title=\"target\" type=\"boolean\"/>\n"
                         "      <attribute id=\"2\" title=\"colored\" type=\"boolean\"/>\n"
                         "      <attribute id=\"3\" title=\"colorCode\" type=\"integer\"/>\n";
            uint64_t attrID = 4;
            for (const auto& attribute : additionalAttributes) {
                STORM_LOG_ASSERT(model->getNumberOfStates() == attribute.second.second.size(), "Additional node attribute vector does not have the same number of values as there are states.");
                outStream << "      <attribute id=\"" << attrID << "\" title=\"" << attribute.first << "\" type=\"" << attributeTypeToString(attribute.second.first) << "\"/>\n";
                attrID++;
            }
            outStream << "    </attributes>\n"
                         "    <nodes>\n";
            // State nodes + intermediary action nodes
            auto initStates = model->getInitialStates();
            auto targetStates = model->getStateLabeling().getStates("target");
            auto singleSelfLoopActions = model->identifySingleSelfLoopActions();
            for (uint_fast64_t state = 0; state < model->getNumberOfStates(); state++) {
                std::string label = initStates[state] ? std::to_string(state) + "_INIT" : (targetStates[state] ? std::to_string(state) + "_TARGET" : std::to_string(state));
                std::string initInfo = initStates[state] ? "true" : "false";
                std::string targetInfo = targetStates[state] ? "true" : "false";
                uint64_t colorCode = useColors? encodeColor(colors[state][0], colors[state][1], colors[state][2]) : 255255255;
                std::string coloredInfo = (colorCode!= 0 && colorCode != 255255255) ? "true" : "false";
                outStream << "      <node id=\"" << state << "\" label=\"" << label << "\"><!-- State node -->\n"
                             "        <attvalues>\n"
                             "          <attvalue for=\"0\" value=\"" << initInfo << "\"/>\n"
                             "          <attvalue for=\"1\" value=\"" << targetInfo << "\"/>\n"
                             "          <attvalue for=\"2\" value=\"" << coloredInfo << "\"/>\n"
                             "          <attvalue for=\"3\" value=\"" << colorCode << "\"/>\n";
                attrID = 4;
                for (auto attribute : additionalAttributes) {
                    outStream << "          <attvalue for=\"" << attrID << "\" value=\"" << attribute.second.second[state] << "\"/>\n";
                    attrID++;
                }
                if (useColors) {
                    outStream << "        </attvalues>"
                                 "        <viz:color r=\"" << colors[state][0] << "\" g=\"" << colors[state][1] << "\" b=\"" << colors[state][2] << "\" a=\"1\" />\n"
                                 "      </node>\n";
                } else {
                    outStream << "        </attvalues>"
                                 "        <viz:color r=\"255\" g=\"255\" b=\"255\" a=\"1\" />\n"
                                 "      </node>\n";
                }

                for (uint_fast64_t action = 0; action < transMatrix.getRowGroupSize(state); action++) {
                    if (!singleSelfLoopActions[transMatrix.getRowGroupIndices()[state] + action]) {
                        outStream << "      <node id=\"" << state << "a" << action << "\" label=\"\"> <!-- Intermediate node -->\n"
                                     "        <viz:color r=\"0\" g=\"0\" b=\"0\" a=\"1\" /> <!-- Black and small-->\n"
                                     "        <viz:size value=\"3\"/>\n"
                                     "      </node>\n";
                    }
                }
            }

            outStream << "    </nodes>\n"
                         "    <edges>\n";
            // Transitions
            for (uint_fast64_t state = 0; state < model->getNumberOfStates(); state++) {
                for (uint_fast64_t action = 0; action < transMatrix.getRowGroupSize(state); action++) {
                    uint64_t rowIndex = transMatrix.getRowGroupIndices()[state] + action;
                    std::string actionLabel;
                    if (model->hasChoiceLabeling() && model->getChoiceLabeling().getLabelsOfChoice(rowIndex).size() == 1) {
                        actionLabel = *model->getChoiceLabeling().getLabelsOfChoice(rowIndex).begin();
                    } else {
                        actionLabel = std::to_string(action);
                    }
                    if (!singleSelfLoopActions[rowIndex]){
                        auto row = transMatrix.getRow(state, action);
                        std::string intermediate = std::to_string(state) + "a" + std::to_string(action);
                        outStream << "      <edge source=\"" << state << "\" target=\"" << intermediate << "\" label=\"" << actionLabel << "\"> <!-- State to Intermediate -->\n";
                        outStream << "        <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                        outStream << "      </edge>\n";
                        for (const auto& entry : row) {
                            outStream << "      <edge source=\"" << intermediate << "\" target=\"" << entry.getColumn() << "\" label=\"" << entry.getValue() << "\"> <!-- Intermediate to Successor State -->\n";
                            outStream << "        <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                            outStream << "      </edge>\n";
                        }
                    } else {
                        outStream << "      <edge source=\"" << state << "\" target=\"" << state << "\" label=\"" << actionLabel << "\"> <!-- Almost-Sure Self Loop -->\n";
                        outStream << "        <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                        outStream << "      </edge>\n";
                    }
                }
            }

            outStream << "    </edges>\n"
                         "  </graph>\n"
                         "</gexf>";

        }

        template<typename ValueType>
        void GEXFExporter<ValueType>::exportDetGEXFToStream(std::shared_ptr<storm::models::sparse::DeterministicModel<ValueType>> model, std::ostream &outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes, storm::storage::SparseMatrix<ValueType> transMatrix) {
            STORM_LOG_ASSERT(colors.empty() || colors.size() == model->getNumberOfStates(), "Color vector has weird size (not zero indicating no coloring, but also not equal to the number of states)!");
            bool useColors = model->getNumberOfStates() == colors.size();

            // Preamble
            outStream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                         "<gexf xmlns=\"http://gexf.net/1.3\"\n"
                         "      xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                         "      xmlns:viz=\"http://gexf.net/1.3/viz\"\n"
                         "      xsi:schemaLocation=\"http://gexf.net/1.3\n"
                         "                          http://gexf.net/1.3/gexf.xsd\"\n"
                         "      version=\"1.3\">\n"
                         "  <meta>\n"
                         "    <creator>Storm</creator>\n"
                         "    <description>Visualization of a deterministic model</description>\n"
                         "  </meta>\n"
                         "  <graph defaultedgetype=\"directed\">\n"
                         "    <attributes class=\"node\">\n"
                         "      <attribute id=\"0\" title=\"init\" type=\"boolean\"/>\n"
                         "      <attribute id=\"1\" title=\"target\" type=\"boolean\"/>\n"
                         "      <attribute id=\"2\" title=\"colored\" type=\"boolean\"/>\n"
                         "      <attribute id=\"3\" title=\"colorCode\" type=\"integer\"/>\n";
            uint64_t attrID = 4;
            for (const auto& attribute : additionalAttributes) {
                STORM_LOG_ASSERT(model->getNumberOfStates() == attribute.second.second.size(), "Additional node attribute vector does not have the same number of values as there are states.");
                outStream << "      <attribute id=\"" << attrID << "\" title=\"" << attribute.first << "\" type=\"" << attributeTypeToString(attribute.second.first) << "\"/>\n";
                attrID++;
            }
            outStream << "    </attributes>\n"
                         "    <nodes>\n";
            // State nodes + intermediary action nodes
            auto initStates = model->getInitialStates();
            auto targetStates = model->getStateLabeling().getStates("target");
            for (uint_fast64_t state = 0; state < model->getNumberOfStates(); state++) {
                std::string label = initStates[state] ? std::to_string(state) + "_INIT" : (targetStates[state] ? std::to_string(state) + "_TARGET" : std::to_string(state));
                std::string initInfo = initStates[state] ? "true" : "false";
                std::string targetInfo = targetStates[state] ? "true" : "false";
                uint64_t colorCode = useColors? encodeColor(colors[state][0], colors[state][1], colors[state][2]) : 255255255;
                std::string coloredInfo = (colorCode!= 0 && colorCode != 255255255) ? "true" : "false";
                outStream << "      <node id=\"" << state << "\" label=\"" << label << "\"><!-- State node -->\n"
                             "        <attvalues>\n"
                             "          <attvalue for=\"0\" value=\"" << initInfo << "\"/>\n"
                             "          <attvalue for=\"1\" value=\"" << targetInfo << "\"/>\n"
                             "          <attvalue for=\"2\" value=\"" << coloredInfo << "\"/>\n"
                             "          <attvalue for=\"3\" value=\"" << colorCode << "\"/>\n";
                attrID = 4;
                for (auto attribute : additionalAttributes) {
                    outStream << "          <attvalue for=\"" << attrID << "\" value=\"" << attribute.second.second[state] << "\"/>\n";
                    attrID++;
                }
                if (useColors) {
                    outStream << "        </attvalues>"
                                 "        <viz:color r=\"" << colors[state][0] << "\" g=\"" << colors[state][1] << "\" b=\"" << colors[state][2] << "\" a=\"1\" />\n"
                                 "      </node>\n";
                } else {
                    outStream << "        </attvalues>"
                                 "        <viz:color r=\"255\" g=\"255\" b=\"255\" a=\"1\" />\n"
                                 "      </node>\n";
                }
            }

            outStream << "    </nodes>\n"
                         "    <edges>\n";
            // Transitions
            for (uint_fast64_t state = 0; state < model->getNumberOfStates(); state++) {
                auto row = transMatrix.getRow(state);
                for (const auto& entry : row) {
                    outStream << "      <edge source=\"" << state << "\" target=\"" << entry.getColumn() << "\" label=\"" << entry.getValue() << "\"> <!-- State to Successor State -->\n";
                    outStream << "        <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                    outStream << "      </edge>\n";
                }
            }
            outStream << "    </edges>\n"
                         "  </graph>\n"
                         "</gexf>";

        }

        template<typename ValueType>
        std::vector<std::vector<uint64_t>> GEXFExporter<ValueType>::createPomdpObservationColoring(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> pomdp) {
            auto colors = std::vector<std::vector<uint64_t>>();
            uint64_t numberOfObservations = pomdp->getNrObservations();
            auto obsColors = createBasicColoring(numberOfObservations);
            if (colors.size() == numberOfObservations) { // if the number of colors could be provided
                for (uint64_t state = 0; state < pomdp->getNumberOfStates(); state++) {
                    colors.emplace_back(3);
                    auto obs = pomdp->getObservation(state);
                    colors[state][0] = obsColors[obs][0];
                    colors[state][1] = obsColors[obs][1];
                    colors[state][2] = obsColors[obs][2];
                }
            }
            return colors;
        }

        template<typename ValueType>
        std::vector<std::vector<uint64_t>> GEXFExporter<ValueType>::createMaStateTypeColoring(std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> ma) {
            auto colors = std::vector<std::vector<uint64_t>>();
            auto baseColors = createBasicColoring(3);
            for (uint64_t state = 0; state < ma->getNumberOfStates(); state++) {
                colors.emplace_back(3);
                uint64_t stateType;
                if (ma->isHybridState(state)) {
                    stateType = 0;
                } else if (ma->isProbabilisticState(state)) {
                    stateType = 1;
                } else {
                    assert(ma->isMarkovianState(state));
                    stateType = 2;
                }
                colors[state][0] = baseColors[stateType][0];
                colors[state][1] = baseColors[stateType][1];
                colors[state][2] = baseColors[stateType][2];
            }
            return colors;
        }

        template<typename ValueType>
        std::vector<std::vector<uint64_t>>
        GEXFExporter<ValueType>::createSmgPlayerColoring(std::shared_ptr<storm::models::sparse::Smg<ValueType>> smg, const std::map<storage::PlayerIndex, uint64_t>& playerToColorIdx) {
            auto colors = std::vector<std::vector<uint64_t>>();
            auto baseColors = createBasicColoring(playerToColorIdx.size());
            for (uint64_t state = 0; state < smg->getNumberOfStates(); state++) {
                colors.emplace_back(3);
                uint64_t colorIdx = playerToColorIdx.at(smg->getPlayerOfState(state));
                colors[state][0] = baseColors[colorIdx][0];
                colors[state][1] = baseColors[colorIdx][1];
                colors[state][2] = baseColors[colorIdx][2];
            }
            return colors;
        }

        template<typename ValueType>
        std::vector<std::vector<uint64_t>> GEXFExporter<ValueType>::createBasicColoring(uint64_t numberOfColors) {
            // gives the graph a basic coloring with distinct colors that are at least somewhat evenly distributed
            auto colors = std::vector<std::vector<uint64_t>>();
            if (numberOfColors > 256 * 256 * 256 - 2) {// -2 because we want to reserve black and white
                STORM_LOG_INFO("Color Space too small for required number of colors");
                return colors; // return empty vector
            }
            for (uint64_t i = 0; i < numberOfColors; i++){
                colors.emplace_back(3);
            }

            if (numberOfColors <= 1530) {
                // Simple version with colors evenly spaced on the full-saturation-full-brightness rgb spectrum
                // (i.e. one value is always 255, one value 0, and the third is in [0,255])

                // we have six intervals when circling the above-mentioned spectrum
                // we save here where each of the rgb values is UP at 255 (u), DOWN at 0 (d), FALLING (f) or rising (r)
                // I've thought long and hard about how to do this elegantly, but I could not come up with anything much better than this
                auto colorBehaviorOnInterval = std::vector<std::vector<char>>(6);
                colorBehaviorOnInterval[0] = {'u', 'd', 'r'};
                colorBehaviorOnInterval[1] = {'f', 'd', 'u'};
                colorBehaviorOnInterval[2] = {'d', 'r', 'u'};
                colorBehaviorOnInterval[3] = {'d', 'u', 'f'};
                colorBehaviorOnInterval[4] = {'r', 'u', 'd'};
                colorBehaviorOnInterval[5] = {'u', 'f', 'd'};

                uint64_t colorDistance = 1530 / numberOfColors;
                for (uint_fast64_t i = 0; i < numberOfColors; i++) {
                    uint64_t value = i * colorDistance; // translate each color to a value within our spectrum of 1530 colors
                    uint64_t intervalIndex = value / 255; // determine the values interval of the spectrum
                    uint64_t offset = value - intervalIndex * 255; // determine offset within the interval
                    for (auto rgbIndex = 0; rgbIndex < 3; rgbIndex++) { // determine value for r, g, and b values with help of vector above
                        switch (colorBehaviorOnInterval[intervalIndex][rgbIndex]) {
                            case 'u':
                                colors[i][rgbIndex] = 255;
                                break;
                            case 'd':
                                colors[i][rgbIndex] = 0;
                                break;
                            case 'r':
                                colors[i][rgbIndex] = offset;
                                break;
                            case 'f':
                                colors[i][rgbIndex] = 255-offset;
                                break;
                        }
                    }
                }

            } else {
                // We wanna divide the color "cube" with:
                // x "cuts" orthogonal to r axis,
                // y "cuts" orthogonal to g axis,
                // z "cuts" orthogonal to b axis,
                // and take the corner points of the resulting sections (minus black and white bc they are reserved) as our different color values
                // we check value triples with x >= y >= z >= x-1
                // with x,y,z close to each other for a relatively even distribution across the color space
                // but allowing that y and z may be 1 smaller than x to provide intermediate sizes for the color space

                uint64_t x = 1;
                uint64_t y;
                uint64_t z;
                // determine the smallest x/y/z that are sufficient
                for (; x < 255; x++) {
                    if ((x + 2) * (x + 1) * (x + 1) - 2 >= numberOfColors) {
                        y = x - 1;
                        z = x - 1;
                        break;
                    }
                    if ((x + 2) * (x + 2) * (x + 1) - 2 >= numberOfColors) {
                        y = x;
                        z = x - 1;
                        break;
                    }
                    if ((x + 2) * (x + 2) * (x + 2) - 2 >= numberOfColors) {
                        y = x;
                        z = x;
                        break;
                    }
                }
                for (uint64_t state = 0; state < numberOfColors; state++) {
                    auto rgb = adaptedEuclid(state + 1, x, y); // +1 so we skip black
                    colors[state][0] = std::get<0>(rgb) * (255 / (x +
                                                                  1)); // +1 because x is the number of cuts, thus x +1 is the number of blocks in r direction
                    colors[state][1] = std::get<1>(rgb) * (255 / (y + 1));
                    colors[state][2] = std::get<2>(rgb) * (255 / (z + 1));
                }
            }
        }

        template<typename ValueType>
        std::tuple<uint64_t, uint64_t, uint64_t> GEXFExporter<ValueType>::adaptedEuclid(uint64_t i, uint64_t x, uint64_t y) {
            // Maps the i-th state to the color block that is the r-th along the r axis, the g-th along the g axis and the b-th along the b axis of the color "cube"
            // Doing this with an adaptation of the Euclidean algorithm where the number at which you go to the next "digit" varies for each digit
            // r is the lowest digit, g the next higher one, and b the highest one
            // With this, r is in [0, x+1], g is in [0, y+1], b is in [0, z+1]
            // And we have i = r + g * (x+2) + b * (x+2)(y+2)
            // (note that the r axis is divided x times, i.e. into x + 1 sections, but including the start and end, we have x + 2 distinct values of interval borders. Thus, r can take x + 2 values which is the size of the interval)
            uint64_t b = i / ((x + 2) * (y + 2));
            i = i - b * (x + 2) * (y + 2);
            uint64_t g = i / (x + 2);
            i = i - g * (x + 2);
            uint64_t r = i;
            return std::make_tuple(r, g, b);
        }

        template class GEXFExporter<double>;
        template class GEXFExporter<storm::RationalNumber>;
    } // storm
} // exporter