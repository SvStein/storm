#include "storm/models/sparse/Mdp.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/constants.h"
#include "storm/utility/vector.h"

#include "storm/models/sparse/StandardRewardModel.h"

namespace storm {
namespace models {
namespace sparse {

template<typename ValueType, typename RewardModelType>
Mdp<ValueType, RewardModelType>::Mdp(storm::storage::SparseMatrix<ValueType> const& transitionMatrix, storm::models::sparse::StateLabeling const& stateLabeling,
                                     std::unordered_map<std::string, RewardModelType> const& rewardModels, ModelType type)
    : Mdp<ValueType, RewardModelType>(storm::storage::sparse::ModelComponents<ValueType, RewardModelType>(transitionMatrix, stateLabeling, rewardModels),
                                      type) {
    // Intentionally left empty
}

template<typename ValueType, typename RewardModelType>
Mdp<ValueType, RewardModelType>::Mdp(storm::storage::SparseMatrix<ValueType>&& transitionMatrix, storm::models::sparse::StateLabeling&& stateLabeling,
                                     std::unordered_map<std::string, RewardModelType>&& rewardModels, ModelType type)
    : Mdp<ValueType, RewardModelType>(
          storm::storage::sparse::ModelComponents<ValueType, RewardModelType>(std::move(transitionMatrix), std::move(stateLabeling), std::move(rewardModels)),
          type) {
    // Intentionally left empty
}

template<typename ValueType, typename RewardModelType>
Mdp<ValueType, RewardModelType>::Mdp(storm::storage::sparse::ModelComponents<ValueType, RewardModelType> const& components, ModelType type)
    : NondeterministicModel<ValueType, RewardModelType>(type, components) {
    assert(type == storm::models::ModelType::Mdp || type == storm::models::ModelType::Pomdp);
    // Intentionally left empty
}

template<typename ValueType, typename RewardModelType>
Mdp<ValueType, RewardModelType>::Mdp(storm::storage::sparse::ModelComponents<ValueType, RewardModelType>&& components, ModelType type)
    : NondeterministicModel<ValueType, RewardModelType>(type, std::move(components)) {
    assert(type == storm::models::ModelType::Mdp || type == storm::models::ModelType::Pomdp);
    // Intentionally left empty
}

template<class ValueType, typename RewardModelType>
void Mdp<ValueType, RewardModelType>::exportGEFXToStream(std::ostream &outStream,
                                                         std::vector<std::vector<uint64_t>> colors) {
    STORM_LOG_ASSERT(this->getNumberOfStates() == colors.size(), "The color vector's size does not equal the number of states.");

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
                 "    <description>Belief MDP Visualization</description>\n"
                 "  </meta>\n"
                 "  <graph defaultedgetype=\"directed\">\n"
                 "    <nodes>\n";
    // State nodes + intermediary action nodes
    for (auto state = 0; state < this->getNumberOfStates(); state++) {
        outStream << "      <node id=\"" << state << "\" label=\"" << state << "\"><!-- State node -->\n"
                     "        <viz:color r=\"" << colors[state][0] << "\" g=\"" << colors[state][1] << "\" b=\"" << colors[state][2] << "\" a=\"1\" />\n"
                     "      </node>\n";
        for (auto action = 0; action < this->getTransitionMatrix().getRowGroupSize(state); action++) {
            outStream << "<node id=\"" << state << "a" << action << "\" label=\"\"> <!-- Intermediate node -->\n"
                         "        <viz:color r=\"0\" g=\"0\" b=\"0\" a=\"1\" /> <!-- Black and small-->\n"
                         "        <viz:size value=\"3\"/>\n"
                         "      </node>\n";
        }
    }

    outStream << "    </nodes>\n"
                 "    <edges>\n";
    // Transitions
    for (auto state = 0; state < this->getNumberOfStates(); state++) {
        for (auto action = 0; action < this->getTransitionMatrix().getRowGroupSize(state); action++) {
            auto row = this->getTransitionMatrix().getRow(state, action);
            uint64_t rowIndex = this->getTransitionMatrix().getRowGroupIndices()[state] + action;
            std::string actionLabel;
            if (this->hasChoiceLabeling() && this->getChoiceLabeling().getLabelsOfChoice(rowIndex).size() == 1) {
                actionLabel = *this->getChoiceLabeling().getLabelsOfChoice(rowIndex).begin();
            } else {
                actionLabel = std::to_string(action);
            }
            std::string intermediate = std::to_string(state) + "a" + std::to_string(action);
            outStream << "      <edge source=\"" << state << "\" target=\"" << intermediate << "\" label=\"" << actionLabel << "\"> <!-- State to Intermediate -->\n";
            outStream << "          <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
            outStream << "      </edge>\n";
            for (auto entry : row) {
                outStream << "      <edge source=\"" << intermediate << "\" target=\"" << entry.getColumn() << "\" label=\"" << entry.getValue() << "\"> <!-- Intermediate to Successor State -->\n";
                outStream << "          <viz:color r=\"0\" g=\"0\" b=\"0\"/>\n";
                outStream << "      </edge>\n";
            }
        }
    }

    outStream << "    </edges>\n"
                 "  </graph>\n"
                 "</gexf>";

}

template class Mdp<double>;
template class Mdp<storm::RationalNumber>;

template class Mdp<double, storm::models::sparse::StandardRewardModel<storm::Interval>>;
template class Mdp<storm::RationalFunction>;
}  // namespace sparse
}  // namespace models
}  // namespace storm
