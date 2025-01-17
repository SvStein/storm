#include "storm-pomdp-cli/settings/modules/BeliefExplorationSettings.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/ArgumentBuilder.h"

#include "storm/utility/NumberTraits.h"
#include "storm/adapters/RationalNumberAdapter.h"
#include "storm-pomdp/modelchecker/BeliefExplorationPomdpModelCheckerOptions.h"

#include "storm/exceptions/InvalidArgumentException.h"


namespace storm {
    namespace settings {
        namespace modules {
            
            const std::string BeliefExplorationSettings::moduleName = "belexpl";

            const std::string beliefTypeOption = "belieftype";
            const std::string refineOption = "refine";
            const std::string explorationTimeLimitOption = "exploration-time";
            const std::string resolutionOption = "resolution";
            const std::string clipGridResolutionOption = "clip-resolution";
            const std::string sizeThresholdOption = "size-threshold";
            const std::string gapThresholdOption = "gap-threshold";
            const std::string schedulerThresholdOption = "scheduler-threshold";
            const std::string observationThresholdOption = "obs-threshold";
            const std::string numericPrecisionOption = "numeric-precision";
            const std::string triangulationModeOption = "triangulationmode";
            const std::string explHeuristicOption = "expl-heuristic";
            const std::string clippingOption = "use-clipping";
            const std::string cutZeroGapOption = "cut-zero-gap";
            const std::string parametricPreprocessingOption = "par-preprocessing";
            const std::string stateEliminationCutoffOption = "state-elimination-cutoff";
            const std::string alphaVectorOption = "import-alphavec";
            const std::string preProcMinMaxMethodOption = "preproc-minmax";

            BeliefExplorationSettings::BeliefExplorationSettings() : ModuleSettings(moduleName) {
                
                this->addOption(storm::settings::OptionBuilder(moduleName, refineOption, false,"Refines the result bounds until reaching either the goal precision or the refinement step limit").addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("prec","The goal precision.").setDefaultValueDouble(1e-4).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(0.0)).build()).addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("steps","The number of allowed refinement steps (0 means no limit).").setDefaultValueUnsignedInteger(0).makeOptional().build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, explorationTimeLimitOption, false, "Sets after which time no further states shall be explored.").addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("time","In seconds.").build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, resolutionOption, false,"Sets the resolution of the discretization and how it is increased in case of refinement").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("init","the initial resolution (higher means more precise)").setDefaultValueUnsignedInteger(3).addValidatorUnsignedInteger(storm::settings::ArgumentValidatorFactory::createUnsignedGreaterValidator(0)).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("factor","Multiplied to the resolution of refined observations (higher means more precise).").setDefaultValueDouble(2).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterValidator(1)).build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, clipGridResolutionOption, false, "Sets the resolution of the clipping grid").addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("resolution", "the resolution (higher means more precise)").setDefaultValueUnsignedInteger(2).addValidatorUnsignedInteger(storm::settings::ArgumentValidatorFactory::createUnsignedGreaterValidator(0)).build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, observationThresholdOption, false,"Only observations whose score is below this threshold will be refined.").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("init","initial threshold (higher means more precise").setDefaultValueDouble(0.1).addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleRangeValidatorIncluding(0,1)).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("factor","Controlls how fast the threshold is increased in each refinement step (higher means more precise).").setDefaultValueDouble(0.1).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleRangeValidatorIncluding(0,1)).build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, sizeThresholdOption, false,"Sets how many new states are explored or rewired in a refinement step and how this value is increased in case of refinement.").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("init","initial limit (higher means more precise, 0 means automatic choice)").setDefaultValueUnsignedInteger(0).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("factor","Before each step the new threshold is set to the current state count times this number (higher means more precise).").setDefaultValueDouble(4).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(1)).build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, gapThresholdOption, false,"Sets how large the gap between known lower- and upper bounds at a beliefstate needs to be in order to explore").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("init","initial threshold (higher means less precise").setDefaultValueDouble(0.1).addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(0)).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("factor","Multiplied to the gap in each refinement step (higher means less precise).").setDefaultValueDouble(0.25).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleRangeValidatorIncluding(0,1)).build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, schedulerThresholdOption, false,"Sets how much worse a sub-optimal choice can be in order to be included in the relevant explored fragment").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("init","initial threshold (higher means more precise").setDefaultValueDouble(1e-3).addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(0)).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("factor","Multiplied to the threshold in each refinement step (higher means more precise).").setDefaultValueDouble(1).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(1)).build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, numericPrecisionOption, false,"Sets the precision used to determine whether two belief-states are equal.").setIsAdvanced().addArgument(
                        storm::settings::ArgumentBuilder::createDoubleArgument("value","the precision").setDefaultValueDouble(1e-9).makeOptional().addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleRangeValidatorIncluding(0, 1)).build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, triangulationModeOption, false,"Sets how to triangulate beliefs when discretizing.").setIsAdvanced().addArgument(
                        storm::settings::ArgumentBuilder::createStringArgument("value","the triangulation mode").setDefaultValueString("dynamic").addValidatorString(storm::settings::ArgumentValidatorFactory::createMultipleChoiceValidator({"dynamic", "static"})).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, clippingOption, false, "If this is set, unfolding will use  (grid) clipping instead of cut-offs only.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, explHeuristicOption, false,"Sets how to sort the states into the exploration queue.").setIsAdvanced().addArgument(
                        storm::settings::ArgumentBuilder::createStringArgument("value","the exploration heuristic").setDefaultValueString("bfs").addValidatorString(storm::settings::ArgumentValidatorFactory::createMultipleChoiceValidator({"bfs", "lowerBound", "upperBound", "gap", "prob"})).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, beliefTypeOption, false,"Sets number type used to handle probabilities in beliefs").setIsAdvanced().addArgument(
                        storm::settings::ArgumentBuilder::createStringArgument("value","the number type. 'default' is the POMDP datatype").setDefaultValueString("default").addValidatorString(storm::settings::ArgumentValidatorFactory::createMultipleChoiceValidator({"default", "float", "rational"})).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, cutZeroGapOption, false,"Cut beliefs where the gap between over- and underapproximation is 0.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, parametricPreprocessingOption, false, "If this is set, the POMDP will be transformed to a pMC for preprocessing steps.").addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("memoryBound", "number of memory states").setDefaultValueUnsignedInteger(0).addValidatorUnsignedInteger(storm::settings::ArgumentValidatorFactory::createUnsignedGreaterEqualValidator(0)).build()).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("gd-eps", "epsilon for gradient descent").setDefaultValueDouble(1e-6).addValidatorDouble(storm::settings::ArgumentValidatorFactory::createDoubleGreaterEqualValidator(0)).build()).addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("maxInstantiations", "max. number of initial instantiations to use for gradient descent").setDefaultValueUnsignedInteger(1).addValidatorUnsignedInteger(storm::settings::ArgumentValidatorFactory::createUnsignedGreaterEqualValidator(1)).build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, stateEliminationCutoffOption, false, "If this is set, an additional unfolding step for cut-off beliefs is performed.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, alphaVectorOption, false, "Loads a set of alpha vectors that is used for preprocessing.").addArgument(storm::settings::ArgumentBuilder::createStringArgument("filename", "The name of the file containing the alpha vectors").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, preProcMinMaxMethodOption, false,"Sets the method to be used for model checking during pre-processing.").setIsAdvanced().addArgument(storm::settings::ArgumentBuilder::createStringArgument("method","the method to use").setDefaultValueString("svi").addValidatorString(storm::settings::ArgumentValidatorFactory::createMultipleChoiceValidator({"svi", "pi"})).build()).build());
            }

            bool BeliefExplorationSettings::isRefineSet() const {
                return this->getOption(refineOption).getHasOptionBeenSet();
            }

            bool BeliefExplorationSettings::isStateEliminationCutoffSet() const {
                return this->getOption(stateEliminationCutoffOption).getHasOptionBeenSet();
            }

            double BeliefExplorationSettings::getRefinePrecision() const {
                return this->getOption(refineOption).getArgumentByName("prec").getValueAsDouble();
            }
            
            bool BeliefExplorationSettings::isRefineStepLimitSet() const {
                return this->getOption(refineOption).getArgumentByName("steps").getValueAsUnsignedInteger() != 0;
            }
            
            uint64_t BeliefExplorationSettings::getRefineStepLimit() const {
                assert(isRefineStepLimitSet());
                return this->getOption(refineOption).getArgumentByName("steps").getValueAsUnsignedInteger();
            }
            
            bool BeliefExplorationSettings::isExplorationTimeLimitSet() const {
                return this->getOption(explorationTimeLimitOption).getHasOptionBeenSet();
            }
            
            uint64_t BeliefExplorationSettings::getExplorationTimeLimit() const {
                return this->getOption(explorationTimeLimitOption).getArgumentByName("time").getValueAsUnsignedInteger();
            }
            
            uint64_t BeliefExplorationSettings::getResolutionInit() const {
                return this->getOption(resolutionOption).getArgumentByName("init").getValueAsUnsignedInteger();
            }

            uint64_t BeliefExplorationSettings::getClippingGridResolution() const {
                return this->getOption(clipGridResolutionOption).getArgumentByName("resolution").getValueAsUnsignedInteger();
            }

            double BeliefExplorationSettings::getResolutionFactor() const {
                return this->getOption(resolutionOption).getArgumentByName("factor").getValueAsDouble();
            }
            
            uint64_t BeliefExplorationSettings::getSizeThresholdInit() const {
                return this->getOption(sizeThresholdOption).getArgumentByName("init").getValueAsUnsignedInteger();
            }
            
            double BeliefExplorationSettings::getSizeThresholdFactor() const {
                return this->getOption(sizeThresholdOption).getArgumentByName("factor").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getGapThresholdInit() const {
                return this->getOption(gapThresholdOption).getArgumentByName("init").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getGapThresholdFactor() const {
                return this->getOption(gapThresholdOption).getArgumentByName("factor").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getOptimalChoiceValueThresholdInit() const {
                return this->getOption(schedulerThresholdOption).getArgumentByName("init").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getOptimalChoiceValueThresholdFactor() const {
                return this->getOption(schedulerThresholdOption).getArgumentByName("factor").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getObservationScoreThresholdInit() const {
                return this->getOption(observationThresholdOption).getArgumentByName("init").getValueAsDouble();
            }
            
            double BeliefExplorationSettings::getObservationScoreThresholdFactor() const {
                return this->getOption(observationThresholdOption).getArgumentByName("factor").getValueAsDouble();
            }
            
            bool BeliefExplorationSettings::isNumericPrecisionSetFromDefault() const {
                return !this->getOption(numericPrecisionOption).getHasOptionBeenSet() || this->getOption(numericPrecisionOption).getArgumentByName("value").wasSetFromDefaultValue();
            }
            
            double BeliefExplorationSettings::getNumericPrecision() const {
                return this->getOption(numericPrecisionOption).getArgumentByName("value").getValueAsDouble();
            }
            
            bool BeliefExplorationSettings::isDynamicTriangulationModeSet() const {
                return this->getOption(triangulationModeOption).getArgumentByName("value").getValueAsString() == "dynamic";
                
            }
            bool BeliefExplorationSettings::isStaticTriangulationModeSet() const {
                return this->getOption(triangulationModeOption).getArgumentByName("value").getValueAsString() == "static";
            }

            bool BeliefExplorationSettings::isUseClippingSet() const {
                return this->getOption(clippingOption).getHasOptionBeenSet();
            }

            storm::builder::ExplorationHeuristic BeliefExplorationSettings::getExplorationHeuristic() const {
                if(this->getOption(explHeuristicOption).getArgumentByName("value").getValueAsString() == "bfs") {
                    return storm::builder::ExplorationHeuristic::BreadthFirst;
                }
                if(this->getOption(explHeuristicOption).getArgumentByName("value").getValueAsString() == "lowerBound") {
                    return storm::builder::ExplorationHeuristic::LowerBoundPrio;
                }
                if(this->getOption(explHeuristicOption).getArgumentByName("value").getValueAsString() == "upperBound") {
                    return storm::builder::ExplorationHeuristic::UpperBoundPrio;
                }
                if(this->getOption(explHeuristicOption).getArgumentByName("value").getValueAsString() == "gap") {
                    return storm::builder::ExplorationHeuristic::GapPrio;
                }
                if(this->getOption(explHeuristicOption).getArgumentByName("value").getValueAsString() == "prob") {
                    return storm::builder::ExplorationHeuristic::ProbabilityPrio;
                }
                return storm::builder::ExplorationHeuristic::BreadthFirst;
            }

            bool BeliefExplorationSettings::isCutZeroGapSet() const {
                return this->getOption(cutZeroGapOption).getHasOptionBeenSet();
            }

            bool BeliefExplorationSettings::isParametricPreprocessingSet() const {
                return this->getOption(parametricPreprocessingOption).getArgumentByName("memoryBound").getValueAsUnsignedInteger() > 0;
            }

            uint64_t BeliefExplorationSettings::getParametricPreprocessingMemoryBound() const {
                return this->getOption(parametricPreprocessingOption).getArgumentByName("memoryBound").getValueAsUnsignedInteger();
            }

            uint64_t BeliefExplorationSettings::getParametricGDMaxInstantiations() const {
                return this->getOption(parametricPreprocessingOption).getArgumentByName("maxInstantiations").getValueAsUnsignedInteger();
            }

            double BeliefExplorationSettings::getParametricGDEpsilon() const {
                return this->getOption(parametricPreprocessingOption).getArgumentByName("gd-eps").getValueAsDouble();
            }

            storm::solver::MinMaxMethod BeliefExplorationSettings::getPreProcMinMaxMethod() const {
                if(this->getOption(preProcMinMaxMethodOption).getArgumentByName("method").getValueAsString() == "svi") {
                    return storm::solver::MinMaxMethod::SoundValueIteration;
                }
                if(this->getOption(preProcMinMaxMethodOption).getArgumentByName("method").getValueAsString() == "pi") {
                    return storm::solver::MinMaxMethod::PolicyIteration;
                }
            }

            template<typename ValueType>
            void BeliefExplorationSettings::setValuesInOptionsStruct(storm::pomdp::modelchecker::BeliefExplorationPomdpModelCheckerOptions<ValueType>& options) const {
                options.refine = isRefineSet();
                options.refinePrecision = storm::utility::convertNumber<ValueType>(getRefinePrecision());
                if (isRefineStepLimitSet()) {
                    options.refineStepLimit = getRefineStepLimit();
                } else {
                    options.refineStepLimit = boost::none;
                }
                if (isExplorationTimeLimitSet()) {
                    options.explorationTimeLimit = getExplorationTimeLimit();
                } else {
                    options.explorationTimeLimit = boost::none;
                }
                options.clippingGridRes = getClippingGridResolution();
                options.resolutionInit = getResolutionInit();
                options.resolutionFactor = storm::utility::convertNumber<ValueType>(getResolutionFactor());
                options.sizeThresholdInit = getSizeThresholdInit();
                options.sizeThresholdFactor = storm::utility::convertNumber<ValueType>(getSizeThresholdFactor());
                options.gapThresholdInit = storm::utility::convertNumber<ValueType>(getGapThresholdInit());
                options.gapThresholdFactor = storm::utility::convertNumber<ValueType>(getGapThresholdFactor());
                options.optimalChoiceValueThresholdInit = storm::utility::convertNumber<ValueType>(getOptimalChoiceValueThresholdInit());
                options.optimalChoiceValueThresholdFactor = storm::utility::convertNumber<ValueType>(getOptimalChoiceValueThresholdFactor());
                options.obsThresholdInit = storm::utility::convertNumber<ValueType>(getObservationScoreThresholdInit());
                options.obsThresholdIncrementFactor = storm::utility::convertNumber<ValueType>(getObservationScoreThresholdFactor());
                options.useGridClipping = isUseClippingSet();
                options.useStateEliminationCutoff = isStateEliminationCutoffSet();

                options.useParametricPreprocessing = isParametricPreprocessingSet();
                options.paramMemBound = getParametricPreprocessingMemoryBound();
                options.paramGDEps = getParametricGDEpsilon();
                options.paramGDMaxInstantiations = getParametricGDMaxInstantiations();
                
                options.numericPrecision = storm::utility::convertNumber<ValueType>(getNumericPrecision());
                if (storm::NumberTraits<ValueType>::IsExact) {
                    if (isNumericPrecisionSetFromDefault()) {
                        STORM_LOG_WARN_COND(storm::utility::isZero(options.numericPrecision), "Setting numeric precision to zero because exact arithmethic is used.");
                        options.numericPrecision = storm::utility::zero<ValueType>();
                    } else {
                        STORM_LOG_WARN_COND(storm::utility::isZero(options.numericPrecision), "A non-zero numeric precision was set although exact arithmethic is used. Results might be inexact.");
                    }
                }
                options.dynamicTriangulation = isDynamicTriangulationModeSet();

                options.explorationHeuristic = getExplorationHeuristic();

                options.cutZeroGap = isCutZeroGapSet();

                options.preProcMinMaxMethod = getPreProcMinMaxMethod();
            }

            bool BeliefExplorationSettings::isBeliefTypeSetFromDefault() const {
                return this->getOption(beliefTypeOption).getArgumentByName("value").getValueAsString() != "default";
            }

            storm::pomdp::BeliefNumberType BeliefExplorationSettings::getBeliefType() const {
                if(this->getOption(beliefTypeOption).getArgumentByName("value").getValueAsString() == "default") {
                    return storm::pomdp::Default;
                }
                if(this->getOption(beliefTypeOption).getArgumentByName("value").getValueAsString() == "float") {
                    return storm::pomdp::Float;
                }
                if(this->getOption(beliefTypeOption).getArgumentByName("value").getValueAsString() == "rational") {
                    return storm::pomdp::Rational;
                }
                STORM_LOG_WARN("Number Type for belief unknown, use default.");
                return storm::pomdp::Default;
            }

            bool BeliefExplorationSettings::isAlphaVectorProcessingSet() const {
                return this->getOption(alphaVectorOption).getHasOptionBeenSet();
            }
            std::string BeliefExplorationSettings::getAlphaVectorFileName() const {
                return this->getOption(alphaVectorOption).getArgumentByName("filename").getValueAsString();
            }

            template void BeliefExplorationSettings::setValuesInOptionsStruct<double>(storm::pomdp::modelchecker::BeliefExplorationPomdpModelCheckerOptions<double>& options) const;
            template void BeliefExplorationSettings::setValuesInOptionsStruct<storm::RationalNumber>(storm::pomdp::modelchecker::BeliefExplorationPomdpModelCheckerOptions<storm::RationalNumber>& options) const;

            
            
        } // namespace modules
    } // namespace settings
} // namespace storm
