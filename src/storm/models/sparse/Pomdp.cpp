#include "storm/models/sparse/Pomdp.h"

#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace models {
namespace sparse {

template<typename ValueType, typename RewardModelType>
Pomdp<ValueType, RewardModelType>::Pomdp(storm::storage::SparseMatrix<ValueType> const &transitionMatrix,
                                         storm::models::sparse::StateLabeling const &stateLabeling,
                                         std::unordered_map<std::string, RewardModelType> const &rewardModels)
    : Mdp<ValueType, RewardModelType>(transitionMatrix, stateLabeling, rewardModels, storm::models::ModelType::Pomdp) {
    computeNrObservations();
}

template<typename ValueType, typename RewardModelType>
Pomdp<ValueType, RewardModelType>::Pomdp(storm::storage::SparseMatrix<ValueType> &&transitionMatrix, storm::models::sparse::StateLabeling &&stateLabeling,
                                         std::unordered_map<std::string, RewardModelType> &&rewardModels)
    : Mdp<ValueType, RewardModelType>(transitionMatrix, stateLabeling, rewardModels, storm::models::ModelType::Pomdp) {
    computeNrObservations();
}

template<typename ValueType, typename RewardModelType>
Pomdp<ValueType, RewardModelType>::Pomdp(storm::storage::sparse::ModelComponents<ValueType, RewardModelType> const &components, bool canonicFlag)
    : Mdp<ValueType, RewardModelType>(components, storm::models::ModelType::Pomdp),
      observations(components.observabilityClasses.value()),
      canonicFlag(canonicFlag),
      observationValuations(components.observationValuations) {
    computeNrObservations();
}

template<typename ValueType, typename RewardModelType>
Pomdp<ValueType, RewardModelType>::Pomdp(storm::storage::sparse::ModelComponents<ValueType, RewardModelType> &&components, bool canonicFlag)
    : Mdp<ValueType, RewardModelType>(components, storm::models::ModelType::Pomdp),
      observations(components.observabilityClasses.value()),
      canonicFlag(canonicFlag),
      observationValuations(components.observationValuations) {
    computeNrObservations();
}

template<typename ValueType, typename RewardModelType>
void Pomdp<ValueType, RewardModelType>::printModelInformationToStream(std::ostream &out) const {
    this->printModelInformationHeaderToStream(out);
    out << "Choices: \t" << this->getNumberOfChoices() << '\n';
    out << "Observations: \t" << this->nrObservations << '\n';
    this->printModelInformationFooterToStream(out);
}

template<typename ValueType, typename RewardModelType>
void Pomdp<ValueType, RewardModelType>::computeNrObservations() {
    uint64_t highestEntry = 0;
    for (uint32_t entry : observations) {
        if (entry > highestEntry) {
            highestEntry = entry;
        }
    }
    nrObservations = highestEntry + 1;  // Smallest entry should be zero.
    // In debug mode, ensure that every observability is used.
}

template<typename ValueType, typename RewardModelType>
uint32_t Pomdp<ValueType, RewardModelType>::getObservation(uint64_t state) const {
    return observations.at(state);
}

template<typename ValueType, typename RewardModelType>
uint64_t Pomdp<ValueType, RewardModelType>::getNrObservations() const {
    return nrObservations;
}

template<typename ValueType, typename RewardModelType>
uint64_t Pomdp<ValueType, RewardModelType>::getMaxNrStatesWithSameObservation() const {
    std::map<uint32_t, uint64_t> counts;
    for (auto const &obs : observations) {
        auto insertionRes = counts.emplace(obs, 1ull);
        if (!insertionRes.second) {
            ++insertionRes.first->second;
        }
    }
    uint64_t result = 0;
    for (auto const &count : counts) {
        result = std::max(result, count.second);
    }
    return result;
}

template<typename ValueType, typename RewardModelType>
std::vector<uint32_t> const &Pomdp<ValueType, RewardModelType>::getObservations() const {
    return observations;
}

template<typename ValueType, typename RewardModelType>
void Pomdp<ValueType, RewardModelType>::updateObservations(std::vector<uint32_t> &&newObservations, bool preservesCanonicity) {
    observations = std::move(newObservations);
    computeNrObservations();
    setIsCanonic(isCanonic() && preservesCanonicity);
}

template<typename ValueType, typename RewardModelType>
std::string Pomdp<ValueType, RewardModelType>::additionalDotStateInfo(uint64_t state) const {
    return "<" + std::to_string(getObservation(state)) + ">";
}

template<typename ValueType, typename RewardModelType>
std::vector<uint64_t> Pomdp<ValueType, RewardModelType>::getStatesWithObservation(uint32_t observation) const {
    std::vector<uint64_t> result;
    for (uint64_t state = 0; state < this->getNumberOfStates(); ++state) {
        if (this->getObservation(state) == observation) {
            result.push_back(state);
        }
    }
    return result;
}

template<typename ValueType, typename RewardModelType>
bool Pomdp<ValueType, RewardModelType>::hasObservationValuations() const {
    return static_cast<bool>(observationValuations);
}

template<typename ValueType, typename RewardModelType>
storm::storage::sparse::StateValuations const &Pomdp<ValueType, RewardModelType>::getObservationValuations() const {
    return observationValuations.value();
}

template<typename ValueType, typename RewardModelType>
std::optional<storm::storage::sparse::StateValuations> const &Pomdp<ValueType, RewardModelType>::getOptionalObservationValuations() const {
    return observationValuations;
}

template<typename ValueType, typename RewardModelType>
bool Pomdp<ValueType, RewardModelType>::isCanonic() const {
    return canonicFlag;
}

template<typename ValueType, typename RewardModelType>
void Pomdp<ValueType, RewardModelType>::setIsCanonic(bool newValue) {
    this->canonicFlag = newValue;
}

template<typename ValueType, typename RewardModelType>
bool Pomdp<ValueType, RewardModelType>::isPartiallyObservable() const {
    return true;
}

template<typename ValueType, typename RewardModelType>
std::size_t Pomdp<ValueType, RewardModelType>::hash() const {
    std::size_t seed = 0;
    boost::hash_combine(seed, sparse::Model<ValueType, RewardModelType>::hash());
    boost::hash_combine(seed, boost::hash_range(observations.begin(), observations.end()));
    return seed;
}

template<class ValueType, typename RewardModelType>
void Pomdp<ValueType, RewardModelType>::toJuliaOutput(std::ostream &outStream) {
    // Prep the action names (don't want to rely on action labels)
    // Actions are not necessarily labeled but states with the same observation always have the same action order
    // TODO maybe add the option to use action labels if existent
    std::map<std::pair<uint_fast64_t, uint_fast64_t>, std::string> actionNameMapping =
        std::map<std::pair<uint_fast64_t, uint_fast64_t>, std::string>();  // state + action to action name
    std::set<std::string> actionNames = std::set<std::string>();
    std::string actionNamesEnumeration;
    for (uint_fast64_t state = 0; state < this->getNumberOfStates(); state++) {
        auto obs = observations[state];
        for (uint_fast64_t action = 0; action < this->getTransitionMatrix().getRowGroupSize(state); action++) {
            std::string actionName = std::to_string(obs) + "_" + std::to_string(action);
            actionNameMapping[{state, action}] = actionName;
            actionNames.insert(actionName);
        }
    }
    bool first = true;
    for (const auto &name : actionNames) {
        if (first) {
            first = false;
            actionNamesEnumeration += "\"" + name + "\"";
        } else {
            actionNamesEnumeration += ", \"" + name + "\"";
        }
    }

    // Prep a vector of all states with a given observation
    std::map<uint32_t, std::set<uint_fast64_t>> obsToStates = std::map<uint32_t, std::set<uint_fast64_t>>();
    for (uint_fast64_t obs = 0; obs < nrObservations; obs++) {
        obsToStates[obs] = std::set<uint_fast64_t>();
    }
    for (uint_fast64_t state = 0; state < this->getNumberOfStates(); state++) {
        obsToStates[observations[state]].insert(state);
    }

    // First part
    outStream << "using QuickPOMDPs: QuickPOMDP\n"
                 "using POMDPTools: Deterministic, SparseCat\n\n"
                 "pomdp = QuickPOMDP(\n"
                 "    states = range(0, length="
              << std::to_string(this->getNumberOfStates())
              << "),\n"
                 "    actions = ["
              << actionNamesEnumeration
              << "],\n"
                 "    observations = range(0, length="
              << std::to_string(nrObservations) << "),\n\n";

    // Transition function
    outStream << "    transition = function (s, a)\n";
    uint_fast64_t numberOfStates = this->getNumberOfStates();
    for (uint_fast64_t state = 0; state < numberOfStates; state++) {
        if (state == 0) {
            outStream << "        if s == " << std::to_string(state) << "\n";
        } else if (state == numberOfStates - 1) {
            outStream << "        else\n";
        } else {
            outStream << "        elseif s == " << std::to_string(state) << "\n";
        }

        uint_fast64_t rowGroupSize = this->getTransitionMatrix().getRowGroupSize(state);
        for (uint_fast64_t action = 0; action < rowGroupSize; action++) {
            if (action == 0) {
                outStream << "            if a == \"" << actionNameMapping[{state, action}] << "\"\n";
            } else if (action == rowGroupSize - 1) {
                outStream << "            else\n";
            } else {
                outStream << "            elseif a == \"" << actionNameMapping[{state, action}] << "\"\n";
            }

            if (this->getTransitionMatrix().getRow(state, action).getNumberOfEntries() == 1) {
                uint_fast64_t succ = this->getTransitionMatrix().getRow(state, action).begin()->getColumn();
                STORM_LOG_ASSERT(storm::utility::isOne(this->getTransitionMatrix().getRow(state, action).begin()->getValue()),
                                 "There is only one entry. Why is it not equal to one?");
                outStream << "                return Deterministic(" << std::to_string(succ) << ")\n";
            } else {
                std::string succsVecString = "[";
                std::string probsVecString = "[";
                first = true;
                for (const auto &entry : this->getTransitionMatrix().getRow(state, action)) {
                    if (first) {
                        first = false;
                        succsVecString += std::to_string(entry.getColumn());
                        probsVecString += storm::utility::to_string(entry.getValue());
                    } else {
                        succsVecString += ", " + std::to_string(entry.getColumn());
                        probsVecString += ", " + storm::utility::to_string(entry.getValue());
                    }
                }
                succsVecString += "]";
                probsVecString += "]";
                outStream << "                return SparseCat(" << succsVecString << ", " << probsVecString << ")\n";
            }
            if (action == rowGroupSize - 1) {
                outStream << "            end\n";
            }
        }
        if (state == numberOfStates - 1) {
            outStream << "        end\n";
        }
    }
    outStream << "    end,\n\n";

    // Observation function
    outStream << "    observation = function (a, sp)\n";
    for (uint32_t obs = 0; obs < nrObservations; obs++) {
        std::string obsStatesVecString = "[";
        first = true;
        for (auto state : obsToStates[obs]) {
            if (first) {
                obsStatesVecString += std::to_string(state);
                first = false;
            } else {
                obsStatesVecString += ", " + std::to_string(state);
            }
        }
        obsStatesVecString += "]";
        if (obs == 0) {  // first obs
            outStream << "        if sp in " << obsStatesVecString << "\n"
                      << "            return Deterministic(" << std::to_string(obs) << ")\n";
        } else if (obs == nrObservations - 1) {  // last obs
            outStream << "        else\n"
                      << "            return Deterministic(" << std::to_string(obs) << ")\n"
                      << "        end\n"
                      << "    end,\n\n";
        } else {  // all other obs
            outStream << "        elseif sp in " << obsStatesVecString << "\n"
                      << "            return Deterministic(" << std::to_string(obs) << ")\n";
        }
    }

    // Reward function
    STORM_LOG_ASSERT(this->hasUniqueRewardModel(), "This output only supports a single reward model :(");
    auto rewModel = this->getUniqueRewardModel();

    bool sRews = rewModel.hasStateRewards();
    bool saRews = rewModel.hasStateActionRewards();
    bool tRews = rewModel.hasTransitionRewards();

    STORM_LOG_ASSERT(sRews || saRews || tRews, "Why does this have no rewards :o");

    if (tRews) {
        // Use the (s, a, sp) version of reward function
        outStream << "    reward = function (s, a, sp)\n";
        for (uint_fast64_t state = 0; state < this->getNumberOfStates(); state++) {
            if (state == 0) {
                outStream << "        if s == " << std::to_string(state) << "\n";
            } else if (state == this->getNumberOfStates() - 1) {
                outStream << "        else\n";
            } else {
                outStream << "        elseif s == " << std::to_string(state) << "\n";
            }
            auto rowGroupIndices = this->getTransitionMatrix().getRowGroupIndices(state);
            typename RewardModelType::ValueType stateRew = sRews ? rewModel.getStateReward(state) : storm::utility::zero<typename RewardModelType::ValueType>();
            auto rowGroupSize = this->getTransitionMatrix().getRowGroupSize(state);
            for (uint_fast64_t action = 0; action < rowGroupSize; action++) {
                typename RewardModelType::ValueType stateActionRew =
                    saRews ? rewModel.getStateActionReward(rowGroupIndices[action]) : storm::utility::zero<typename RewardModelType::ValueType>();
                if (action == 0) {
                    auto actionString = "\"" + actionNameMapping[{state, action}] + "\"";
                    outStream << "            if a == " << actionString << "\n";
                } else if (action == rowGroupSize - 1) {
                    outStream << "            else\n";
                } else {
                    auto actionString = "\"" + actionNameMapping[{state, action}] + "\"";
                    outStream << "            elseif a == " << actionString << "\n";
                }

                auto transRews = std::map<uint_fast64_t, typename RewardModelType::ValueType>();
                for (const auto &entry : rewModel.getTransitionRewardMatrix().getRow(state, action)) {
                    auto succState = entry.getColumn();
                    transRews[succState] = entry.getValue();
                }
                uint_fast64_t succNr = 0;
                auto rowSize = this->getTransitionMatrix().getRow(state, action).getNumberOfEntries();
                for (const auto &entry : this->getTransitionMatrix().getRow(state, action)) {
                    // we need this extra iteration through the transition matrix
                    // to account for transitions with transition reward zero but non-zero state or stateaction rewards,
                    // these would otherwise not get into the rewards map due to the sparse nature of the reward matrix
                    auto succState = entry.getColumn();
                    typename RewardModelType::ValueType finalReward = stateRew + stateActionRew;
                    if (transRews.contains(succState)) {
                        finalReward += transRews[succState];
                    }
                    if (succNr == 0) {
                        outStream << "                if sp == " << std::to_string(succState) << "\n"
                                  << "                    return Deterministic(" << storm::utility::to_string(finalReward) << ")\n";
                    } else if (succNr == rowSize - 1) {
                        outStream << "                else\n"
                                  << "                    return Deterministic(" << storm::utility::to_string(finalReward) << ")\n"
                                  << "                end\n";
                    } else {
                        outStream << "                elseif sp == " << std::to_string(succState) << "\n"
                                  << "                    return Deterministic(" << storm::utility::to_string(finalReward) << ")\n";
                    }
                }
            }
            outStream << "            end\n";
        }
        outStream << "        end\n";
    } else {
        // Use the (s, a) version of reward function
        outStream << "    reward = function (s, a)\n";
        if (rewModel.hasOnlyStateRewards()) {
            // only need to differentiate between states
            for (uint_fast64_t state = 0; state < this->getNumberOfStates(); state++) {
                if (state == 0) {
                    outStream << "        if s == " << std::to_string(state) << "\n"
                              << "            return Deterministic(" << storm::utility::to_string(rewModel.getStateReward(state)) << ")\n";
                } else if (state == this->getNumberOfStates() - 1) {
                    outStream << "        else\n"
                              << "            return Deterministic(" << storm::utility::to_string(rewModel.getStateReward(state)) << ")\n"
                              << "        end\n";
                } else {
                    outStream << "        elseif s == " << std::to_string(state) << "\n"
                              << "            return Deterministic(" << storm::utility::to_string(rewModel.getStateReward(state)) << ")\n";
                }
            }
        } else {
            // differentiate between states AND actions
            for (uint_fast64_t state = 0; state < this->getNumberOfStates(); state++) {
                if (state == 0) {
                    outStream << "        if s == " << std::to_string(state) << "\n";
                } else if (state == this->getNumberOfStates() - 1) {
                    outStream << "        else\n";
                } else {
                    outStream << "        elseif s == " << std::to_string(state) << "\n";
                }
                auto rowGroupIndices = this->getTransitionMatrix().getRowGroupIndices(state);
                auto rowGroupSize = this->getTransitionMatrix().getRowGroupSize(state);
                for (uint_fast64_t action = 0; action < rowGroupSize; action++) {
                    auto reward = storm::utility::zero<typename RewardModelType::ValueType>();
                    if (sRews) {
                        reward += rewModel.getStateReward(state);
                    }
                    reward += rewModel.getStateActionReward(rowGroupIndices[action]);
                    if (action == 0) {
                        auto actionString = "\"" + actionNameMapping[{state, action}] + "\"";
                        outStream << "            if a == " << actionString << "\n"
                                  << "                return Deterministic(" << storm::utility::to_string(reward) << ")\n";
                    } else if (action == rowGroupSize - 1) {
                        outStream << "            else\n"
                                  << "                return Deterministic(" << storm::utility::to_string(reward) << ")\n"
                                  << "            end\n";
                    } else {
                        auto actionString = "\"" + actionNameMapping[{state, action}] + "\"";
                        outStream << "            elseif a == " << actionString << "\n"
                                  << "                return Deterministic(" << storm::utility::to_string(reward) << ")\n";
                    }
                }
            }
            outStream << "        end\n";
        }
    }
    outStream << "    end,\n\n";

    // Initial State
    storm::storage::BitVector initStates = this->getInitialStates();
    STORM_LOG_ASSERT(initStates.getNumberOfSetBits() == 1, "I was told there is only ever one init state :(");
    auto initState = initStates.getNextSetIndex(0);
    outStream << "    initialstate = Deterministic(" << std::to_string(initState) << "),\n";

    // Finish up
    outStream << ");";
}

template class Pomdp<double>;
template class Pomdp<storm::RationalNumber>;
template class Pomdp<double, storm::models::sparse::StandardRewardModel<storm::Interval>>;
template class Pomdp<storm::RationalFunction>;
template class Pomdp<storm::Interval>;
}  // namespace sparse
}  // namespace models
}  // namespace storm
