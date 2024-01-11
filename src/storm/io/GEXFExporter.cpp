//
// Created by spook on 11.01.24.
//

#include "GEXFExporter.h"

#include <boost/container/flat_map.hpp>
#include "storm/io/file.h"

namespace storm {
    namespace exporter {

        template<class ValueType, typename BeliefType>
        uint64_t GEXFExporter<ValueType, BeliefType>::encodeColor(uint64_t r, uint64_t g, uint64_t b) {
            return (r*1000000 + g*1000 + b);
        }

        template<class ValueType, typename BeliefType>
        std::string GEXFExporter<ValueType, BeliefType>::attributeTypeToString(GEXFAttributeType attributeType) {
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

        template<class ValueType, typename BeliefType>
        void GEXFExporter<ValueType, BeliefType>::exportGEFXToStream(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp, std::ostream &outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes) {
            // TODO make colors optional
            STORM_LOG_ASSERT(mdp->getNumberOfStates() == colors.size(), "The color vector's size does not equal the number of states.");
            // TODO maybe some basic format check if the strings given as values for the attributes adhere to the given attribute type
    
            // Preamble
            outStream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                         "<gexf xmlns=\"http://gexf.net/1.3\"\n"
                         "      xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                         "      xmlns:viz=\"http://gexf.net/1.3/viz\"\n"
                         "      xsi:schemaLocation=\"http://gexf.net/1.3\n"
                         "                          http://gexf.net/1.3/gexf.xsd\"\n"
                         "      version=\"1.3\">\n"
                         "  <meta>\n"
                         "    <creator>Spook</creator>\n"
                         "    <description>MDP Visualization</description>\n"
                         "  </meta>\n"
                         "  <graph defaultedgetype=\"directed\">\n"
                         "    <attributes class=\"node\">\n"
                         "      <attribute id=\"0\" title=\"init\" type=\"boolean\"/>\n"
                         "      <attribute id=\"1\" title=\"target\" type=\"boolean\"/>\n"
                         "      <attribute id=\"2\" title=\"colored\" type=\"boolean\"/>\n"
                         "      <attribute id=\"3\" title=\"colorCode\" type=\"integer\"/>\n";
            uint64_t attrID = 4;
            for (auto attribute : additionalAttributes) {
                STORM_LOG_ASSERT(mdp->getNumberOfStates() == attribute.second.second.size(), "Additional node attribute vector does not have the same number of values as there are states.");
                outStream << "      <attribute id=\"" << attrID << "\" title=\"" << attribute.first << "\" type=\"" << attributeTypeToString(attribute.second.first) << "\"/>\n";
                attrID++;
            }
            outStream << "    </attributes>\n"
                         "    <nodes>\n";
            // State nodes + intermediary action nodes
            auto initStates = mdp->getInitialStates();
            auto targetStates = mdp->getStateLabeling().getStates("target");
            auto almostSureSelfLoops = mdp->identifyAlmostSureSelfLoops();
            for (auto state = 0; state < mdp->getNumberOfStates(); state++) {
                std::string label = initStates[state] ? std::to_string(state) + "_INIT" : (targetStates[state] ? std::to_string(state) + "_TARGET" : std::to_string(state));
                std::string initInfo = initStates[state] ? "true" : "false";
                std::string targetInfo = targetStates[state] ? "true" : "false";
                uint64_t colorCode = encodeColor(colors[state][0], colors[state][1], colors[state][2]);
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
                outStream << "        </attvalues>"
                             "        <viz:color r=\"" << colors[state][0] << "\" g=\"" << colors[state][1] << "\" b=\"" << colors[state][2] << "\" a=\"1\" />\n"
                             "      </node>\n";
                for (auto action = 0; action < mdp->getTransitionMatrix().getRowGroupSize(state); action++) {
                    if (!almostSureSelfLoops[mdp->getTransitionMatrix().getRowGroupIndices()[state] + action]) {
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
            for (auto state = 0; state < mdp->getNumberOfStates(); state++) {
                for (auto action = 0; action < mdp->getTransitionMatrix().getRowGroupSize(state); action++) {
                    uint64_t rowIndex = mdp->getTransitionMatrix().getRowGroupIndices()[state] + action;
                    std::string actionLabel;
                    if (mdp->hasChoiceLabeling() && mdp->getChoiceLabeling().getLabelsOfChoice(rowIndex).size() == 1) {
                        actionLabel = *mdp->getChoiceLabeling().getLabelsOfChoice(rowIndex).begin();
                    } else {
                        actionLabel = std::to_string(action);
                    }
                    if (!almostSureSelfLoops[rowIndex]){
                        auto row = mdp->getTransitionMatrix().getRow(state, action);
                        std::string intermediate = std::to_string(state) + "a" + std::to_string(action);
                        outStream << "      <edge source=\"" << state << "\" target=\"" << intermediate << "\" label=\"" << actionLabel << "\"> <!-- State to Intermediate -->\n";
                        outStream << "        <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                        outStream << "      </edge>\n";
                        for (auto entry : row) {
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

        template<class ValueType, typename BeliefType>
        void GEXFExporter<ValueType, BeliefType>::determineOgColors(std::vector<std::vector<uint64_t>> &stateColors) {
            // gives the og beliefMDP graph a coloring
            // with colors evenly spaced on the full-saturation-full-brightness rgb spectrum
            // (i.e. one value is always 255, one value 0, and the third is in [0,255]
            // depending on how many colors will be needed, possibly add saturation / brightness variation too

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

            auto numberOfColors = stateColors.size();
            uint64_t colorDistance = 1530 / numberOfColors;
            for (auto i = 0; i < numberOfColors; i++) {
                uint64_t value = i * colorDistance; // translate each color to a value within our spectrum of 1530 colors
                uint64_t intervalIndex = value / 255; // determine the values interval of the spectrum
                uint64_t offset = value - intervalIndex * 255; // determine offset within the interval
                for (auto rgbIndex = 0; rgbIndex < 3; rgbIndex++) { // determine value for r, g, and b values with help of vector above
                    switch (colorBehaviorOnInterval[intervalIndex][rgbIndex]) {
                        case 'u':
                            stateColors[i][rgbIndex] = 255;
                            break;
                        case 'd':
                            stateColors[i][rgbIndex] = 0;
                            break;
                        case 'r':
                            stateColors[i][rgbIndex] = offset;
                            break;
                        case 'f':
                            stateColors[i][rgbIndex] = 255-offset;
                            break;
                    }
                }
            }
        }

        template<class ValueType, typename BeliefType>
        void GEXFExporter<ValueType, BeliefType>::determineNumberOfEpochs(std::vector<std::string> &numbersOfEpochs, std::vector<std::string> &maxSingleStateEpochNumbers, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox) {
            // TODO Determine which values make sense for 0 and 1 later
            numbersOfEpochs[0] = "0";
            numbersOfEpochs[1] = "0";
            maxSingleStateEpochNumbers[0] = "0";
            maxSingleStateEpochNumbers[1] = "0";
            for (auto unfBMDPState = 2; unfBMDPState < numbersOfEpochs.size(); unfBMDPState++) {
                boost::container::flat_map<uint64_t , BeliefType> unfBeliefUnfStateDistribution = (underApprox ? unfCheckingResult.beliefManagerUnder->getBelief(unfCheckingResult.mdpStateToBeliefIdMapUnder[unfBMDPState]) : unfCheckingResult.beliefManagerOver->getBelief(unfCheckingResult.mdpStateToBeliefIdMapOver[unfBMDPState]));
                std::set<ValueType> totalEpochs;
                std::map<uint64_t, std::set<ValueType>> epochsPerState;
                for (auto it = unfBeliefUnfStateDistribution.begin(); it != unfBeliefUnfStateDistribution.end(); it++) {
                    // find out what og state the unfState refers to if we ignore its epoch
                    auto ogStateEpochPair = unfoldingInfo.newStateToStateEpoch[it->first];
                    auto ogState = ogStateEpochPair.first;
                    auto epoch = ogStateEpochPair.second;
                    // sort it into the belief distribution over og states
                    totalEpochs.insert(epoch);
                    if (epochsPerState.find(ogState) == epochsPerState.end()) {
                        epochsPerState[ogState] = std::set<ValueType>();
                    }
                    epochsPerState[ogState].insert(epoch);
                }
                numbersOfEpochs[unfBMDPState] = std::to_string(totalEpochs.size());
                auto setSizeCmp = [] (const std::pair<uint64_t, std::set<ValueType>> &v1, const std::pair<uint64_t, std::set<ValueType>> &v2) -> bool {
                    return v1.second.size() < v2.second.size();
                };
                maxSingleStateEpochNumbers[unfBMDPState] = std::to_string(std::max_element(epochsPerState.begin(), epochsPerState.end(), setSizeCmp)->second.size());
            }
        }


        template<class ValueType, typename BeliefType>
        void GEXFExporter<ValueType, BeliefType>::determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox) {
            // this is one possible implementation of deciding when to color an unfolded state the same as a non-unfolded state, others may be added / implemented as well
            // here: unfBelief ~= ogBelief iff for all s in supp(ogBelief): ogBelief(s) = unfBelief(s)
            // where unfBelief(s) = sum_{(s,e) in supp(unfBelief)} unfBelief((s,e))

            for (auto rgbIndex = 0; rgbIndex < 3; rgbIndex++){
                unfStateColors[0][rgbIndex] = ogStateColors[0][rgbIndex];
                unfStateColors[1][rgbIndex] = ogStateColors[1][rgbIndex];
            }

            // for each unfBelief / state in the unfBeliefMdp
            for (auto unfBMDPState = 2; unfBMDPState < unfStateColors.size(); unfBMDPState++) {
                // set up distribution over the og states
                auto unfBeliefOgStateDistribution = boost::container::flat_map<uint64_t , BeliefType>(); // TODO is BeliefType here what BeliefValueType is in BeliefManager?
                // grab distribution over the unf states
                boost::container::flat_map<uint64_t , BeliefType> unfBeliefUnfStateDistribution = (underApprox ? unfCheckingResult.beliefManagerUnder->getBelief(unfCheckingResult.mdpStateToBeliefIdMapUnder[unfBMDPState]) : unfCheckingResult.beliefManagerOver->getBelief(unfCheckingResult.mdpStateToBeliefIdMapOver[unfBMDPState])); // TODO getBelief is private, how to resolve?
                // for each (unfState,prob) in the support of the unfBelief
                for (auto it = unfBeliefUnfStateDistribution.begin(); it != unfBeliefUnfStateDistribution.end(); it++) {
                    // find out what og state the unfState refers to if we ignore its epoch
                    uint64_t ogState = unfoldingInfo.newStateToStateEpoch[it->first].first;
                    // sort it into the belief distribution over og states
                    if (unfBeliefOgStateDistribution.find(ogState) == unfBeliefOgStateDistribution.end()) {
                        unfBeliefOgStateDistribution[ogState] = it->second;
                    } else {
                        unfBeliefOgStateDistribution[ogState] += it->second;
                    }
                }

                // we now have computed the unfolded belief's distribution over the og states
                // now we need to see if it matches any of the og beliefs' distributions over the og states

                bool found = false;
                // for each og BMDP state
                for (auto ogBMDPState = 2; ogBMDPState < ogStateColors.size(); ogBMDPState++){
                    // grab belief distribution
                    boost::container::flat_map<uint64_t , BeliefType> ogBeliefOgStateDistribution = (underApprox ? ogCheckingResult.beliefManagerUnder->getBelief(ogCheckingResult.mdpStateToBeliefIdMapUnder[ogBMDPState]) : ogCheckingResult.beliefManagerOver->getBelief(ogCheckingResult.mdpStateToBeliefIdMapOver[ogBMDPState])); // TODO getBelief is private, how to resolve?
                    // quick check to maybe avoid having to compare the maps in detail, no clue if that helps or just creates more overhead
                    if (ogBeliefOgStateDistribution.size() != unfBeliefOgStateDistribution.size()) {
                        continue;
                    }
                    bool match = true;
                    // for each (ogState, ogProb)
                    for (auto it = ogBeliefOgStateDistribution.begin(); it != ogBeliefOgStateDistribution.end(); it++){
                        // if the unfBelief does not have the state in its support or maps it to a different probability, the two beliefs dont match
                        // TODO handle potential floating point inaccuracies?
                        if (unfBeliefOgStateDistribution.find(it->first) == unfBeliefOgStateDistribution.end() || unfBeliefOgStateDistribution[it->first] != it->second) {
                            match = false;
                            break;
                        }
                    }

                    if (!match) {
                        // if they don't match, skip to next og BMDPstate
                        continue;
                    } else {
                        // if they do match, set same color and break loop (i am assuming that no two BMDP states have the exact same belief here) // TODO floating point inaccuracies?
                        found = true;
                        for (auto rgbIndex = 0; rgbIndex < 3; rgbIndex++){
                            unfStateColors[unfBMDPState][rgbIndex] = ogStateColors[ogBMDPState][rgbIndex];
                        }
                        break;
                    }
                }
                if (!found) {
                    // if it does not match any og BMDP, set color to white
                    for (auto rgbIndex = 0; rgbIndex < 3; rgbIndex++){
                        unfStateColors[unfBMDPState][rgbIndex] = 255;
                    }
                }
            }
        }

        template<class ValueType, typename BeliefType>
        void GEXFExporter<ValueType, BeliefType>::createGEFXOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo, std::string filename) {
            // TODO for now, this assumes both under and overapproximation have been applied
            // underapprox part
            uint64_t ogUnderStateNumber = ogCheckingResult.beliefMdpUnder->getNumberOfStates();
            auto ogUnderColors = std::vector<std::vector<uint64_t>>();
            ogUnderColors.reserve(ogUnderStateNumber);
            for (auto i = 0; i < ogUnderStateNumber; i++){
                ogUnderColors.push_back(std::vector<uint64_t>(3));
            }
            uint64_t unfUnderStateNumber = unfCheckingResult.beliefMdpUnder->getNumberOfStates();
            auto unfUnderColors = std::vector<std::vector<uint64_t>>();
            unfUnderColors.reserve(unfUnderStateNumber);
            for (auto i = 0; i < unfUnderStateNumber; i++){
                unfUnderColors.push_back(std::vector<uint64_t>(3));
            }

            determineOgColors(ogUnderColors);
            determineUnfColors(ogUnderColors, unfUnderColors, ogCheckingResult, unfCheckingResult, unfoldingInfo, true);

            auto unfUnderNumberOfEpochs = std::vector<std::string>(unfUnderStateNumber);
            auto unfUnderMaxSingleStateEpochNumbers = std::vector<std::string>(unfUnderStateNumber);
            determineNumberOfEpochs(unfUnderNumberOfEpochs, unfUnderMaxSingleStateEpochNumbers, unfCheckingResult, unfoldingInfo, true);
            auto unfUnderExtraAttr = std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>();
            unfUnderExtraAttr["numberOfEpochs"] = std::make_pair(GEXFAttributeType::GEXF_integer, unfUnderNumberOfEpochs);
            unfUnderExtraAttr["maxSingleStateEpochNumber"] = std::make_pair(GEXFAttributeType::GEXF_integer, unfUnderMaxSingleStateEpochNumbers);

            std::ofstream stream;
            storm::utility::openFile(filename + "_ogUnder.gexf", stream);
            exportGEFXToStream(ogCheckingResult.beliefMdpUnder, stream, ogUnderColors);
            storm::utility::closeFile(stream);
            stream.clear();
            storm::utility::openFile(filename + "_unfUnder.gexf", stream);
            exportGEFXToStream(unfCheckingResult.beliefMdpUnder, stream, unfUnderColors, unfUnderExtraAttr);
            storm::utility::closeFile(stream);
            stream.clear();

            // Overapprox part
            uint64_t ogOverStateNumber = ogCheckingResult.beliefMdpOver->getNumberOfStates();
            std::vector<std::vector<uint64_t>> ogOverColors = std::vector<std::vector<uint64_t>>();
            ogOverColors.reserve(ogOverStateNumber);
            for (auto i = 0; i < ogOverStateNumber; i++){
                ogOverColors.push_back(std::vector<uint64_t>(3));
            }
            uint64_t unfOverStateNumber = unfCheckingResult.beliefMdpOver->getNumberOfStates();
            auto unfOverColors = std::vector<std::vector<uint64_t>>();
            unfOverColors.reserve(unfOverStateNumber);
            for (auto i = 0; i < unfOverStateNumber; i++){
                unfOverColors.push_back(std::vector<uint64_t>(3));
            }

            determineOgColors(ogOverColors);
            determineUnfColors(ogOverColors, unfOverColors, ogCheckingResult, unfCheckingResult, unfoldingInfo, false);

            auto unfOverNumberOfEpochs = std::vector<std::string>(unfOverStateNumber);
            auto unfOverMaxSingleStateEpochNumbers = std::vector<std::string>(unfOverStateNumber);
            determineNumberOfEpochs(unfOverNumberOfEpochs, unfOverMaxSingleStateEpochNumbers, unfCheckingResult, unfoldingInfo, false);
            auto unfOverExtraAttr = std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>();
            unfOverExtraAttr["numberOfEpochs"] = std::make_pair(GEXFAttributeType::GEXF_integer, unfOverNumberOfEpochs);
            unfOverExtraAttr["maxSingleStateEpochNumber"] = std::make_pair(GEXFAttributeType::GEXF_integer, unfOverMaxSingleStateEpochNumbers);

            storm::utility::openFile(filename + "_ogOver.gexf", stream);
            exportGEFXToStream(ogCheckingResult.beliefMdpOver, stream, ogOverColors);
            storm::utility::closeFile(stream);
            stream.clear();
            storm::utility::openFile(filename + "_unfOver.gexf", stream);
            exportGEFXToStream(unfCheckingResult.beliefMdpOver, stream, unfOverColors, unfOverExtraAttr);
            storm::utility::closeFile(stream);
        }

        template class GEXFExporter<double, double>;
        template class GEXFExporter<double, storm::RationalNumber>;
        template class GEXFExporter<storm::RationalNumber, double>;
        template class GEXFExporter<storm::RationalNumber, storm::RationalNumber>;
    } // storm
} // exporter