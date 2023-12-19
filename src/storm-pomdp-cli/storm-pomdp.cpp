#include "storm/utility/initialize.h"

#include "storm-pomdp-cli/settings/modules/BeliefExplorationSettings.h"
#include "storm-pomdp-cli/settings/modules/POMDPSettings.h"
#include "storm-pomdp-cli/settings/modules/QualitativePOMDPAnalysisSettings.h"
#include "storm-pomdp-cli/settings/modules/ToParametricSettings.h"
#include "storm/settings/modules/DebugSettings.h"
#include "storm/settings/modules/GeneralSettings.h"

#include "storm-pomdp-cli/settings/PomdpSettings.h"
#include "storm/analysis/GraphConditions.h"

#include "storm-cli-utilities/cli.h"
#include "storm-cli-utilities/model-handling.h"

#include "storm-pomdp/analysis/FormulaInformation.h"
#include "storm-pomdp/analysis/IterativePolicySearch.h"
#include "storm-pomdp/analysis/JaniBeliefSupportMdpGenerator.h"
#include "storm-pomdp/analysis/OneShotPolicySearch.h"
#include "storm-pomdp/analysis/QualitativeAnalysisOnGraphs.h"
#include "storm-pomdp/analysis/UniqueObservationStates.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelChecker.h"
#include "storm-pomdp/transformer/ApplyFiniteSchedulerToPomdp.h"
#include "storm-pomdp/transformer/BinaryPomdpTransformer.h"
#include "storm-pomdp/transformer/BoundUnfolder.h"
#include "storm-pomdp/transformer/GlobalPOMDPSelfLoopEliminator.h"
#include "storm-pomdp/transformer/GlobalPomdpMecChoiceEliminator.h"
#include "storm-pomdp/transformer/KnownProbabilityTransformer.h"
#include "storm-pomdp/transformer/MakePOMDPCanonic.h"
#include "storm-pomdp/transformer/PomdpMemoryUnfolder.h"
#include "storm/api/storm.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/utility/NumberTraits.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/Stopwatch.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"

#include "storm-pars/transformer/ParametricTransformer.h"

#include <typeinfo>

namespace storm {
namespace pomdp {
namespace cli {

/// Perform preprocessings based on the graph structure (if requested or necessary). Return true, if some preprocessing has been done
template<typename ValueType, storm::dd::DdType DdType>
bool performPreprocessing(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>>& pomdp, storm::pomdp::analysis::FormulaInformation& formulaInfo,
                          storm::logic::Formula const& formula) {
    auto const& pomdpSettings = storm::settings::getModule<storm::settings::modules::POMDPSettings>();
    bool preprocessingPerformed = false;
    if (pomdpSettings.isSelfloopReductionSet()) {
        storm::transformer::GlobalPOMDPSelfLoopEliminator<ValueType> selfLoopEliminator(*pomdp);
        if (selfLoopEliminator.preservesFormula(formula)) {
            STORM_PRINT_AND_LOG("Eliminating self-loop choices ...");
            uint64_t oldChoiceCount = pomdp->getNumberOfChoices();
            pomdp = selfLoopEliminator.transform();
            STORM_PRINT_AND_LOG(oldChoiceCount - pomdp->getNumberOfChoices() << " choices eliminated through self-loop elimination.\n");
            preprocessingPerformed = true;
        } else {
            STORM_PRINT_AND_LOG("Not eliminating self-loop choices as it does not preserve the formula.\n");
        }
    }
    if (pomdpSettings.isQualitativeReductionSet() && formulaInfo.isNonNestedReachabilityProbability()) {
        storm::analysis::QualitativeAnalysisOnGraphs<ValueType> qualitativeAnalysis(*pomdp);
        STORM_PRINT_AND_LOG("Computing states with probability 0 ...");
        storm::storage::BitVector prob0States = qualitativeAnalysis.analyseProb0(formula.asProbabilityOperatorFormula());
        std::cout << prob0States << '\n';
        STORM_PRINT_AND_LOG(" done. " << prob0States.getNumberOfSetBits() << " states found.\n");
        STORM_PRINT_AND_LOG("Computing states with probability 1 ...");
        storm::storage::BitVector prob1States = qualitativeAnalysis.analyseProb1(formula.asProbabilityOperatorFormula());
        STORM_PRINT_AND_LOG(" done. " << prob1States.getNumberOfSetBits() << " states found.\n");
        storm::pomdp::transformer::KnownProbabilityTransformer<ValueType> kpt = storm::pomdp::transformer::KnownProbabilityTransformer<ValueType>();
        pomdp = kpt.transform(*pomdp, prob0States, prob1States);
        // Update formulaInfo to changes from Preprocessing
        formulaInfo.updateTargetStates(*pomdp, std::move(prob1States));
        formulaInfo.updateSinkStates(*pomdp, std::move(prob0States));
        preprocessingPerformed = true;
    }
    return preprocessingPerformed;
}

template<typename ValueType>
void printResult(ValueType const& lowerBound, ValueType const& upperBound) {
    if (lowerBound == upperBound) {
        if (storm::utility::isInfinity(lowerBound)) {
            STORM_PRINT_AND_LOG("inf");
        } else {
            STORM_PRINT_AND_LOG(lowerBound);
        }
    } else if (storm::utility::isInfinity<ValueType>(-lowerBound)) {
        if (storm::utility::isInfinity(upperBound)) {
            STORM_PRINT_AND_LOG("[-inf, inf] (width=inf)");
        } else {
            // Only upper bound is known
            STORM_PRINT_AND_LOG("≤ " << upperBound);
        }
    } else if (storm::utility::isInfinity(upperBound)) {
        STORM_PRINT_AND_LOG("≥ " << lowerBound);
    } else {
        STORM_PRINT_AND_LOG("[" << lowerBound << ", " << upperBound << "] (width=" << ValueType(upperBound - lowerBound) << ")");
    }
    if (storm::NumberTraits<ValueType>::IsExact) {
        STORM_PRINT_AND_LOG(" (approx. ");
        double roundedLowerBound =
            storm::utility::isInfinity<ValueType>(-lowerBound) ? -storm::utility::infinity<double>() : storm::utility::convertNumber<double>(lowerBound);
        double roundedUpperBound =
            storm::utility::isInfinity(upperBound) ? storm::utility::infinity<double>() : storm::utility::convertNumber<double>(upperBound);
        printResult(roundedLowerBound, roundedUpperBound);
        STORM_PRINT_AND_LOG(")");
    }
}

MemlessSearchOptions fillMemlessSearchOptionsFromSettings() {
    storm::pomdp::MemlessSearchOptions options;
    auto const& qualSettings = storm::settings::getModule<storm::settings::modules::QualitativePOMDPAnalysisSettings>();

    options.onlyDeterministicStrategies = qualSettings.isOnlyDeterministicSet();
    uint64_t loglevel = 0;
    // TODO a big ugly, but we have our own loglevels (for technical reasons)
    if (storm::utility::getLogLevel() == l3pp::LogLevel::INFO) {
        loglevel = 1;
    } else if (storm::utility::getLogLevel() == l3pp::LogLevel::DEBUG) {
        loglevel = 2;
    } else if (storm::utility::getLogLevel() == l3pp::LogLevel::TRACE) {
        loglevel = 3;
    }
    options.setDebugLevel(loglevel);
    options.validateEveryStep = qualSettings.validateIntermediateSteps();
    options.validateResult = qualSettings.validateFinalResult();

    options.pathVariableType = storm::pomdp::pathVariableTypeFromString(qualSettings.getLookaheadType());

    if (qualSettings.isExportSATCallsSet()) {
        options.setExportSATCalls(qualSettings.getExportSATCallsPath());
    }

    return options;
}

template<typename ValueType>
void performQualitativeAnalysis(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> const& origpomdp,
                                storm::pomdp::analysis::FormulaInformation const& formulaInfo, storm::logic::Formula const& formula) {
    auto const& qualSettings = storm::settings::getModule<storm::settings::modules::QualitativePOMDPAnalysisSettings>();
    auto const& coreSettings = storm::settings::getModule<storm::settings::modules::CoreSettings>();
    std::stringstream sstr;
    origpomdp->printModelInformationToStream(sstr);
    STORM_LOG_INFO(sstr.str());
    STORM_LOG_THROW(formulaInfo.isNonNestedReachabilityProbability(), storm::exceptions::NotSupportedException,
                    "Qualitative memoryless scheduler search is not implemented for this property type.");
    STORM_LOG_TRACE("Run qualitative preprocessing...");
    storm::models::sparse::Pomdp<ValueType> pomdp(*origpomdp);
    storm::analysis::QualitativeAnalysisOnGraphs<ValueType> qualitativeAnalysis(pomdp);
    // After preprocessing, this might be done cheaper.
    storm::storage::BitVector surelyNotAlmostSurelyReachTarget = qualitativeAnalysis.analyseProbSmaller1(formula.asProbabilityOperatorFormula());
    pomdp.getTransitionMatrix().makeRowGroupsAbsorbing(surelyNotAlmostSurelyReachTarget);
    storm::storage::BitVector targetStates = qualitativeAnalysis.analyseProb1(formula.asProbabilityOperatorFormula());
    bool computedSomething = false;
    if (qualSettings.isMemlessSearchSet()) {
        computedSomething = true;
        std::shared_ptr<storm::utility::solver::SmtSolverFactory> smtSolverFactory = std::make_shared<storm::utility::solver::Z3SmtSolverFactory>();
        uint64_t lookahead = qualSettings.getLookahead();
        if (lookahead == 0) {
            lookahead = pomdp.getNumberOfStates();
        }
        if (qualSettings.getMemlessSearchMethod() == "one-shot") {
            storm::pomdp::OneShotPolicySearch<ValueType> memlessSearch(pomdp, targetStates, surelyNotAlmostSurelyReachTarget, smtSolverFactory);
            if (qualSettings.isWinningRegionSet()) {
                STORM_LOG_ERROR("Computing winning regions is not supported by the one-shot method.");
            } else {
                bool result = memlessSearch.analyzeForInitialStates(lookahead);
                if (result) {
                    STORM_PRINT_AND_LOG("From initial state, one can almost-surely reach the target.\n");
                } else {
                    STORM_PRINT_AND_LOG("From initial state, one may not almost-surely reach the target .\n");
                }
            }
        } else if (qualSettings.getMemlessSearchMethod() == "iterative") {
            storm::pomdp::MemlessSearchOptions options = fillMemlessSearchOptionsFromSettings();
            storm::pomdp::IterativePolicySearch<ValueType> search(pomdp, targetStates, surelyNotAlmostSurelyReachTarget, smtSolverFactory, options);
            if (qualSettings.isWinningRegionSet()) {
                search.computeWinningRegion(lookahead);
            } else {
                bool result = search.analyzeForInitialStates(lookahead);
                if (result) {
                    STORM_PRINT_AND_LOG("From initial state, one can almost-surely reach the target.");
                } else {
                    // TODO consider adding check for end components to improve this message.
                    STORM_PRINT_AND_LOG("From initial state, one may not almost-surely reach the target.");
                }
            }

            if (qualSettings.isPrintWinningRegionSet()) {
                search.getLastWinningRegion().print();
                std::cout << '\n';
            }
            if (qualSettings.isExportWinningRegionSet()) {
                std::size_t hash = pomdp.hash();
                search.getLastWinningRegion().storeToFile(qualSettings.exportWinningRegionPath(), "model hash: " + std::to_string(hash));
            }

            search.finalizeStatistics();
            if (pomdp.getInitialStates().getNumberOfSetBits() == 1) {
                uint64_t initialState = pomdp.getInitialStates().getNextSetIndex(0);
                uint64_t initialObservation = pomdp.getObservation(initialState);
                // TODO this is inefficient.
                uint64_t offset = 0;
                for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                    if (state == initialState) {
                        break;
                    }
                    if (pomdp.getObservation(state) == initialObservation) {
                        ++offset;
                    }
                }

                if (search.getLastWinningRegion().isWinning(initialObservation, offset)) {
                    STORM_PRINT_AND_LOG("Initial state is safe!\n");
                } else {
                    STORM_PRINT_AND_LOG("Initial state may not be safe.\n");
                }
            } else {
                STORM_LOG_WARN("Output for multiple initial states is incomplete");
            }

            if (coreSettings.isShowStatisticsSet()) {
                STORM_PRINT_AND_LOG("#STATS Number of belief support states: " << search.getLastWinningRegion().beliefSupportStates() << '\n');
                if (qualSettings.computeExpensiveStats()) {
                    auto wbss = search.getLastWinningRegion().computeNrWinningBeliefs();
                    STORM_PRINT_AND_LOG("#STATS Number of winning belief support states: [" << wbss.first << "," << wbss.second << "]");
                }
                search.getStatistics().print();
            }

        } else {
            STORM_LOG_ERROR("This method is not implemented.");
        }
    }
    if (qualSettings.isComputeOnBeliefSupportSet()) {
        computedSomething = true;
        storm::pomdp::qualitative::JaniBeliefSupportMdpGenerator<ValueType> janicreator(pomdp);
        janicreator.generate(targetStates, surelyNotAlmostSurelyReachTarget);
        bool initialOnly = !qualSettings.isWinningRegionSet();
        janicreator.verifySymbolic(initialOnly);
        STORM_PRINT_AND_LOG("Initial state is safe: " << janicreator.isInitialWinning() << "\n");
    }
    STORM_LOG_THROW(computedSomething, storm::exceptions::InvalidSettingsException, "Nothing to be done, did you forget to set a method?");
}

template<typename ValueType, storm::dd::DdType DdType, typename BeliefType>
std::pair<bool, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result> performSpecialAnalysis(
        const shared_ptr<storm::models::sparse::Pomdp<ValueType>> &pomdp, const analysis::FormulaInformation &formulaInfo,
        const logic::Formula &formula) {
    // need this bc the regular analysis does not return the result. using only the belief exploring part here bc that is what im here for
    auto const& pomdpSettings = storm::settings::getModule<storm::settings::modules::POMDPSettings>();
    bool analysisPerformed = false;
    STORM_PRINT_AND_LOG("Exploring the belief MDP... \n");
    auto options = storm::pomdp::modelchecker::BeliefExplorationPomdpModelCheckerOptions<ValueType>(pomdpSettings.isBeliefExplorationDiscretizeSet(),
                                                                                                    pomdpSettings.isBeliefExplorationUnfoldSet());
    auto const& beliefExplorationSettings = storm::settings::getModule<storm::settings::modules::BeliefExplorationSettings>();
    beliefExplorationSettings.setValuesInOptionsStruct(options);
    storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType> checker(pomdp, options);
    auto result = checker.check(formula);
    checker.printStatisticsToStream(std::cout);
    if (storm::utility::resources::isTerminate()) {
        STORM_PRINT_AND_LOG("\nResult till abort: ")
    } else {
        STORM_PRINT_AND_LOG("\nResult: ")
    }
    printResult(result.lowerBound, result.upperBound);
    STORM_PRINT_AND_LOG('\n');
    analysisPerformed = true;
    return std::make_pair(analysisPerformed, result);
}

void determineOgColors(std::vector<std::vector<uint64_t>> &stateColors) {
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

template<typename ValueType, storm::dd::DdType DdType, typename BeliefType> // TODO idk is this the right thing to put here?
void determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox) {
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

template<typename ValueType, storm::dd::DdType DdType, typename BeliefType> // TODO idk is this the right thing to put here?
void createGEFXOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo) {
    // TODO for now, this assumes both under and overapproximation have been applied
    // underapprox part
    uint64_t ogUnderStateNumber = ogCheckingResult.beliefMdpUnder->getNumberOfStates();
    std::vector<std::vector<uint64_t>> ogUnderColors = std::vector<std::vector<uint64_t>>();
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
    determineUnfColors<ValueType, DdType, BeliefType>(ogUnderColors, unfUnderColors, ogCheckingResult, unfCheckingResult, unfoldingInfo, true);

    std::string filename = storm::settings::getModule<storm::settings::modules::POMDPSettings>().getBMDPCompareFilename();

    std::ofstream stream;
    storm::utility::openFile(filename + "_ogUnder.gexf", stream);
    ogCheckingResult.beliefMdpUnder->exportGEFXToStream(stream, ogUnderColors);
    storm::utility::closeFile(stream);
    stream.clear();
    storm::utility::openFile(filename + "_unfUnder.gexf", stream);
    unfCheckingResult.beliefMdpUnder->exportGEFXToStream(stream, unfUnderColors);
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
    determineUnfColors<ValueType, DdType, BeliefType>(ogOverColors, unfOverColors, ogCheckingResult, unfCheckingResult, unfoldingInfo, false);

    storm::utility::openFile("ogOver.gexf", stream);
    ogCheckingResult.beliefMdpOver->exportGEFXToStream(stream, ogOverColors);
    storm::utility::closeFile(stream);
    stream.clear();
    storm::utility::openFile("unfOver.gexf", stream);
    unfCheckingResult.beliefMdpOver->exportGEFXToStream(stream, unfOverColors);
    storm::utility::closeFile(stream);
}

template<typename ValueType, storm::dd::DdType DdType, typename BeliefType>
bool performAnalysis(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> const& pomdp, storm::pomdp::analysis::FormulaInformation const& formulaInfo,
                     storm::logic::Formula const& formula) {
    auto const& pomdpSettings = storm::settings::getModule<storm::settings::modules::POMDPSettings>();
    bool analysisPerformed = false;
    if (pomdpSettings.isBeliefExplorationSet()) {
        STORM_PRINT_AND_LOG("Exploring the belief MDP... \n");
        auto options = storm::pomdp::modelchecker::BeliefExplorationPomdpModelCheckerOptions<ValueType>(pomdpSettings.isBeliefExplorationDiscretizeSet(),
                                                                                                        pomdpSettings.isBeliefExplorationUnfoldSet());
        auto const& beliefExplorationSettings = storm::settings::getModule<storm::settings::modules::BeliefExplorationSettings>();
        beliefExplorationSettings.setValuesInOptionsStruct(options);
        storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType> checker(pomdp, options);
        auto result = checker.check(formula);
        checker.printStatisticsToStream(std::cout);
        if (storm::utility::resources::isTerminate()) {
            STORM_PRINT_AND_LOG("\nResult till abort: ")
        } else {
            STORM_PRINT_AND_LOG("\nResult: ")
        }
        printResult(result.lowerBound, result.upperBound);
        STORM_PRINT_AND_LOG('\n');
        analysisPerformed = true;
    }
    if (pomdpSettings.isQualitativeAnalysisSet()) {
        performQualitativeAnalysis(pomdp, formulaInfo, formula);
        analysisPerformed = true;
    }
    if (pomdpSettings.isCheckFullyObservableSet()) {
        STORM_PRINT_AND_LOG("Analyzing the formula on the fully observable MDP ... ");
        auto resultPtr = storm::api::verifyWithSparseEngine<ValueType>(pomdp->template as<storm::models::sparse::Mdp<ValueType>>(),
                                                                       storm::api::createTask<ValueType>(formula.asSharedPointer(), true));
        if (resultPtr) {
            auto result = resultPtr->template asExplicitQuantitativeCheckResult<ValueType>();
            result.filter(storm::modelchecker::ExplicitQualitativeCheckResult(pomdp->getInitialStates()));
            if (storm::utility::resources::isTerminate()) {
                STORM_PRINT_AND_LOG("\nResult till abort: ")
            } else {
                STORM_PRINT_AND_LOG("\nResult: ")
            }
            printResult(result.getMin(), result.getMax());
            STORM_PRINT_AND_LOG('\n');
        } else {
            STORM_PRINT_AND_LOG("\nResult: Not available.\n");
        }
        analysisPerformed = true;
    }
    return analysisPerformed;
}

template<typename ValueType, storm::dd::DdType DdType>
bool performTransformation(std::shared_ptr<storm::models::sparse::Pomdp<ValueType>>& pomdp, storm::logic::Formula const& formula) {
    auto const& pomdpSettings = storm::settings::getModule<storm::settings::modules::POMDPSettings>();
    auto const& ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    auto const& transformSettings = storm::settings::getModule<storm::settings::modules::ToParametricSettings>();
    bool transformationPerformed = false;
    bool memoryUnfolded = false;
    if (pomdpSettings.getMemoryBound() > 1) {
        STORM_PRINT_AND_LOG("Computing the unfolding for memory bound " << pomdpSettings.getMemoryBound() << " and memory pattern '"
                                                                        << storm::storage::toString(pomdpSettings.getMemoryPattern()) << "' ...");
        storm::storage::PomdpMemory memory = storm::storage::PomdpMemoryBuilder().build(pomdpSettings.getMemoryPattern(), pomdpSettings.getMemoryBound());
        std::cout << memory.toString() << '\n';
        storm::transformer::PomdpMemoryUnfolder<ValueType> memoryUnfolder(*pomdp, memory);
        pomdp = memoryUnfolder.transform();
        STORM_PRINT_AND_LOG(" done.\n");
        pomdp->printModelInformationToStream(std::cout);
        transformationPerformed = true;
        memoryUnfolded = true;
    }

    // From now on the POMDP is considered memoryless

    if (transformSettings.isMecReductionSet()) {
        STORM_PRINT_AND_LOG("Eliminating mec choices ...");
        // Note: Elimination of mec choices only preserves memoryless schedulers.
        uint64_t oldChoiceCount = pomdp->getNumberOfChoices();
        storm::transformer::GlobalPomdpMecChoiceEliminator<ValueType> mecChoiceEliminator(*pomdp);
        pomdp = mecChoiceEliminator.transform(formula);
        STORM_PRINT_AND_LOG(" done.\n");
        STORM_PRINT_AND_LOG(oldChoiceCount - pomdp->getNumberOfChoices() << " choices eliminated through MEC choice elimination.\n");
        pomdp->printModelInformationToStream(std::cout);
        transformationPerformed = true;
    }

    if (transformSettings.isTransformBinarySet() || transformSettings.isTransformSimpleSet()) {
        if (transformSettings.isTransformSimpleSet()) {
            STORM_PRINT_AND_LOG("Transforming the POMDP to a simple POMDP.");
            pomdp = storm::transformer::BinaryPomdpTransformer<ValueType>().transform(*pomdp, true).transformedPomdp;
        } else {
            STORM_PRINT_AND_LOG("Transforming the POMDP to a binary POMDP.");
            pomdp = storm::transformer::BinaryPomdpTransformer<ValueType>().transform(*pomdp, false).transformedPomdp;
        }
        pomdp->printModelInformationToStream(std::cout);
        STORM_PRINT_AND_LOG(" done.\n");
        transformationPerformed = true;
    }

    if (pomdpSettings.isExportToParametricSet()) {
        STORM_PRINT_AND_LOG("Transforming memoryless POMDP to pMC...");
        storm::transformer::ApplyFiniteSchedulerToPomdp<ValueType> toPMCTransformer(*pomdp);
        std::string transformMode = transformSettings.getFscApplicationTypeString();
        auto pmc = toPMCTransformer.transform(storm::transformer::parsePomdpFscApplicationMode(transformMode));
        STORM_PRINT_AND_LOG(" done.\n");
        if (transformSettings.allowPostSimplifications()) {
            STORM_PRINT_AND_LOG("Simplifying pMC...");
            pmc = storm::api::performBisimulationMinimization<storm::RationalFunction>(pmc->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>(),
                                                                                       {formula.asSharedPointer()}, storm::storage::BisimulationType::Strong)
                      ->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
            STORM_PRINT_AND_LOG(" done.\n");
            pmc->printModelInformationToStream(std::cout);
        }
        if (pmc->hasRewardModel() && transformSettings.isConstantRewardsSet()) {
            STORM_PRINT_AND_LOG("Ensuring constant rewards...");
            pmc = storm::transformer::makeRewardsConstant(*(pmc->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>()));
            STORM_PRINT_AND_LOG(" done.\n");
            pmc->printModelInformationToStream(std::cout);
        }
        STORM_PRINT_AND_LOG("Exporting pMC...");
        storm::analysis::ConstraintCollector<storm::RationalFunction> constraints(*pmc);
        auto const& parameterSet = constraints.getVariables();
        std::vector<storm::RationalFunctionVariable> parameters(parameterSet.begin(), parameterSet.end());
        std::vector<std::string> parameterNames;
        for (auto const& parameter : parameters) {
            parameterNames.push_back(parameter.name());
        }
        storm::api::exportSparseModelAsDrn(pmc, pomdpSettings.getExportToParametricFilename(), parameterNames,
                                           !ioSettings.isExplicitExportPlaceholdersDisabled());
        STORM_PRINT_AND_LOG(" done.\n");
        transformationPerformed = true;
    }
    if (transformationPerformed && !memoryUnfolded) {
        STORM_PRINT_AND_LOG("Implicitly assumed restriction to memoryless schedulers for at least one transformation.\n");
    }
    return transformationPerformed;
}

template<typename ValueType, storm::dd::DdType DdType>
void processOptionsWithValueTypeAndDdLib(storm::cli::SymbolicInput const& symbolicInput, storm::cli::ModelProcessingInformation const& mpi) {
    auto const& pomdpSettings = storm::settings::getModule<storm::settings::modules::POMDPSettings>();

    auto model = storm::cli::buildPreprocessExportModelWithValueTypeAndDdlib<DdType, ValueType>(symbolicInput, mpi);
    if (!model) {
        STORM_PRINT_AND_LOG("No input model given.\n");
        return;
    }
    STORM_LOG_THROW(model->getType() == storm::models::ModelType::Pomdp && model->isSparseModel(), storm::exceptions::WrongFormatException,
                    "Expected a POMDP in sparse representation.");

    std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> pomdp = model->template as<storm::models::sparse::Pomdp<ValueType>>();

    std::shared_ptr<storm::logic::Formula const> formula;
    if (!symbolicInput.properties.empty()) {
        formula = symbolicInput.properties.front().getRawFormula();
        STORM_PRINT_AND_LOG("Analyzing property '" << *formula << "'\n");
        STORM_LOG_WARN_COND(symbolicInput.properties.size() == 1,
                            "There is currently no support for multiple properties. All other properties will be ignored.");
    }

    if (!pomdpSettings.isNoCanonicSet()) {
        storm::transformer::MakePOMDPCanonic<ValueType> makeCanonic(*pomdp);
        pomdp = makeCanonic.transform();
    }

    if (pomdpSettings.isAnalyzeUniqueObservationsSet()) {
        STORM_PRINT_AND_LOG("Analyzing states with unique observation ...\n");
        storm::analysis::UniqueObservationStates<ValueType> uniqueAnalysis(*pomdp);
        std::cout << uniqueAnalysis.analyse() << '\n';
    }

    // TODO set/unset flag to create gefx comparison outputs here
    bool compareBMDPs = storm::settings::getModule<storm::settings::modules::POMDPSettings>().isCompareBMDPsSet();
    if (formula && !compareBMDPs) {
        if (formula->asOperatorFormula().getSubformula().isBoundedUntilFormula()) {
            auto unfolder = storm::transformer::BoundUnfolder<ValueType>();
            auto unfoldedStuff = unfolder.unfold(pomdp, *formula.get());
            pomdp = unfoldedStuff.pomdp;
            formula = std::make_shared<storm::logic::ProbabilityOperatorFormula const>(unfoldedStuff.formula);
        }
        auto formulaInfo = storm::pomdp::analysis::getFormulaInformation(*pomdp, *formula);
        STORM_LOG_THROW(!formulaInfo.isUnsupported(), storm::exceptions::InvalidPropertyException,
                        "The formula '" << *formula << "' is not supported by storm-pomdp.");

        storm::utility::Stopwatch sw(true);
        // Note that formulaInfo contains state-based information which potentially needs to be updated during preprocessing
        if (performPreprocessing<ValueType, DdType>(pomdp, formulaInfo, *formula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for graph-based POMDP (pre-)processing: " << sw << ".\n");
            pomdp->printModelInformationToStream(std::cout);
        }

        sw.restart();
        if (performTransformation<ValueType, DdType>(pomdp, *formula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP transformation(s): " << sw << "s.\n");
        }

        sw.restart();
        if (performAnalysis<ValueType, DdType, ValueType>(pomdp, formulaInfo, *formula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP analysis: " << sw << "s.\n");
        }
    } else if (formula && compareBMDPs){
        STORM_LOG_ASSERT(formula->asOperatorFormula().getSubformula().isBoundedUntilFormula(), "Not bounded Until formula");
        auto unfolder = storm::transformer::BoundUnfolder<ValueType>();
        auto unfoldedStuff = unfolder.unfold(pomdp, *formula.get());
        auto unfpomdp = unfoldedStuff.pomdp;
        shared_ptr<const logic::Formula> unfformula = std::make_shared<storm::logic::ProbabilityOperatorFormula const>(unfoldedStuff.formula);
        // drop reward bound from formula ((currently assumes Pmin/max = ? operator bounded until formula))
        std::string droppedBoundString = "";
        STORM_LOG_ASSERT(formula->isProbabilityOperatorFormula(), "Dropping of Reward bounds currently only implemented for ProbabilityOperator Formulas");

        if (formula->asProbabilityOperatorFormula().getOptimalityType() == storm::solver::OptimizationDirection::Maximize){
            droppedBoundString += "Pmax=? ";
        } else {
            droppedBoundString += "Pmin=? ";
        }
        droppedBoundString += "[" + formula->asOperatorFormula().getSubformula().asBoundedUntilFormula().getLeftSubformula().toString() + " U " + formula->asOperatorFormula().getSubformula().asBoundedUntilFormula().getRightSubformula().toString() + "]";
        std::vector<storm::jani::Property> propertyVector = storm::api::parseProperties(droppedBoundString);
        storm::logic::ProbabilityOperatorFormula newFormula = storm::api::extractFormulasFromProperties(propertyVector).front()->asProbabilityOperatorFormula();
        formula = std::make_shared<storm::logic::ProbabilityOperatorFormula const>(newFormula);

        std::cout << std::endl << formula->toString() << std::endl;

        auto formulaInfo = storm::pomdp::analysis::getFormulaInformation(*pomdp, *formula);
        STORM_LOG_THROW(!formulaInfo.isUnsupported(), storm::exceptions::InvalidPropertyException,
                        "The formula '" << *formula << "' is not supported by storm-pomdp.");
        auto unfformulaInfo = storm::pomdp::analysis::getFormulaInformation(*pomdp, *unfformula);
        STORM_LOG_THROW(!unfformulaInfo.isUnsupported(), storm::exceptions::InvalidPropertyException,
                        "The formula '" << *unfformula << "' is not supported by storm-pomdp.");

        // TODO run both here
        storm::utility::Stopwatch sw(true);
        // Note that formulaInfo contains state-based information which potentially needs to be updated during preprocessing
        if (performPreprocessing<ValueType, DdType>(pomdp, formulaInfo, *formula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for graph-based POMDP (pre-)processing: " << sw << ".\n");
            pomdp->printModelInformationToStream(std::cout);
        }

        sw.restart();
        if (performTransformation<ValueType, DdType>(pomdp, *formula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP transformation(s): " << sw << "s.\n");
        }

        sw.restart();
        auto performedResult = performSpecialAnalysis<ValueType, DdType, ValueType>(pomdp, formulaInfo, *formula);
        if (performedResult.first) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP analysis: " << sw << "s.\n") ;
        }

        if (performPreprocessing<ValueType, DdType>(unfpomdp, unfformulaInfo, *unfformula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for graph-based POMDP (pre-)processing: " << sw << ".\n");
            pomdp->printModelInformationToStream(std::cout);
        }

        sw.restart();
        if (performTransformation<ValueType, DdType>(unfpomdp, *unfformula)) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP transformation(s): " << sw << "s.\n");
        }

        auto unfperformedResult = performSpecialAnalysis<ValueType, DdType, ValueType>(unfpomdp, unfformulaInfo, *unfformula);
        sw.restart();
        if (unfperformedResult.first) {
            sw.stop();
            STORM_PRINT_AND_LOG("Time for POMDP analysis: " << sw << "s.\n");
        }

        createGEFXOutputs<ValueType, DdType, ValueType>(performedResult.second, unfperformedResult.second, unfoldedStuff);
    } else {
        STORM_LOG_WARN("Nothing to be done. Did you forget to specify a formula?");
    }
}

template<storm::dd::DdType DdType>
void processOptionsWithDdLib(storm::cli::SymbolicInput const& symbolicInput, storm::cli::ModelProcessingInformation const& mpi) {
    STORM_LOG_ERROR_COND(mpi.buildValueType == mpi.verificationValueType,
                         "Build value type differs from verification value type. Will ignore Verification value type.");
    switch (mpi.buildValueType) {
        case storm::cli::ModelProcessingInformation::ValueType::FinitePrecision:
            processOptionsWithValueTypeAndDdLib<double, DdType>(symbolicInput, mpi);
            break;
        case storm::cli::ModelProcessingInformation::ValueType::Exact:
            STORM_LOG_THROW(DdType == storm::dd::DdType::Sylvan, storm::exceptions::UnexpectedException,
                            "Exact arithmetic is only supported with Dd library Sylvan.");
            processOptionsWithValueTypeAndDdLib<storm::RationalNumber, storm::dd::DdType::Sylvan>(symbolicInput, mpi);
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected ValueType for model building.");
    }
}

void processOptions() {
    auto symbolicInput = storm::cli::parseSymbolicInput();
    storm::cli::ModelProcessingInformation mpi;
    std::tie(symbolicInput, mpi) = storm::cli::preprocessSymbolicInput(symbolicInput);
    switch (mpi.ddType) {
        case storm::dd::DdType::CUDD:
            processOptionsWithDdLib<storm::dd::DdType::CUDD>(symbolicInput, mpi);
            break;
        case storm::dd::DdType::Sylvan:
            processOptionsWithDdLib<storm::dd::DdType::Sylvan>(symbolicInput, mpi);
            break;
        default:
            STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Unexpected Dd Type.");
    }
}
}  // namespace cli
}  // namespace pomdp
}  // namespace storm

/*!
 * Entry point for the pomdp backend.
 *
 * @param argc The argc argument of main().
 * @param argv The argv argument of main().
 * @return Return code, 0 if successfull, not 0 otherwise.
 */
int main(const int argc, const char** argv) {
    // try {
    storm::utility::setUp();
    storm::cli::printHeader("Storm-pomdp", argc, argv);
    storm::settings::initializePomdpSettings("Storm-POMDP", "storm-pomdp");

    bool optionsCorrect = storm::cli::parseOptions(argc, argv);
    if (!optionsCorrect) {
        return -1;
    }
    storm::utility::Stopwatch totalTimer(true);
    storm::cli::setUrgentOptions();

    // Invoke storm-pomdp with obtained settings
    storm::pomdp::cli::processOptions();

    totalTimer.stop();
    if (storm::settings::getModule<storm::settings::modules::ResourceSettings>().isPrintTimeAndMemorySet()) {
        storm::cli::printTimeAndMemoryStatistics(totalTimer.getTimeInMilliseconds());
    }

    // All operations have now been performed, so we clean up everything and terminate.
    storm::utility::cleanUp();
    return 0;
    // } catch (storm::exceptions::BaseException const &exception) {
    //    STORM_LOG_ERROR("An exception caused Storm-pomdp to terminate. The message of the exception is: " << exception.what());
    //    return 1;
    //} catch (std::exception const &exception) {
    //    STORM_LOG_ERROR("An unexpected exception occurred and caused Storm-pomdp to terminate. The message of this exception is: " << exception.what());
    //    return 2;
    //}
}
