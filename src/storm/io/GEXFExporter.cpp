//
// Created by spook on 11.01.24.
//

#include "GEXFExporter.h"
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
        void GEXFExporter<ValueType, BeliefType>::exportGEXFToStream(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp, std::ostream &outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes) {
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

        template class GEXFExporter<double, double>;
        template class GEXFExporter<double, storm::RationalNumber>;
        template class GEXFExporter<storm::RationalNumber, double>;
        template class GEXFExporter<storm::RationalNumber, storm::RationalNumber>;
    } // storm
} // exporter