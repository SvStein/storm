#include "FaultTreeSettings.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/SettingMemento.h"
#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Argument.h"
#include "storm/exceptions/IllegalArgumentValueException.h"
#include "storm/exceptions/InvalidSettingsException.h"
#include "storm/parser/CSVParser.h"

namespace storm {
    namespace settings {
        namespace modules {

            const std::string FaultTreeSettings::moduleName = "dft";
            const std::string FaultTreeSettings::symmetryReductionOptionName = "symmetryreduction";
            const std::string FaultTreeSettings::symmetryReductionOptionShortName = "symred";
            const std::string FaultTreeSettings::modularisationOptionName = "modularisation";
            const std::string FaultTreeSettings::disableDCOptionName = "disabledc";
            const std::string FaultTreeSettings::relevantEventsOptionName = "relevantevents";
            const std::string FaultTreeSettings::approximationErrorOptionName = "approximation";
            const std::string FaultTreeSettings::approximationErrorOptionShortName = "approx";
            const std::string FaultTreeSettings::approximationHeuristicOptionName = "approximationheuristic";
            const std::string FaultTreeSettings::firstDependencyOptionName = "firstdep";
#ifdef STORM_HAVE_Z3
            const std::string FaultTreeSettings::solveWithSmtOptionName = "smt";
#endif

            FaultTreeSettings::FaultTreeSettings() : ModuleSettings(moduleName) {
                this->addOption(storm::settings::OptionBuilder(moduleName, symmetryReductionOptionName, false, "Exploit symmetric structure of model.").setShortName(symmetryReductionOptionShortName).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, modularisationOptionName, false, "Use modularisation (not applicable for expected time).").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, disableDCOptionName, false, "Disable Dont Care propagation.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, firstDependencyOptionName, false, "Avoid non-determinism by always taking the first possible dependency.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, relevantEventsOptionName, false, "Specifies the relevant events from the DFT.")
                    .addArgument(storm::settings::ArgumentBuilder::createStringArgument("values", "A comma separated list of relevant events. 'all' marks all events as relevant, The default '' or 'none' mark only the top level event as relevant.").setDefaultValueString("").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, approximationErrorOptionName, false, "Approximation error allowed.").setShortName(approximationErrorOptionShortName).addArgument(storm::settings::ArgumentBuilder::createDoubleArgument("error", "The relative approximation error to use.").addValidatorDouble(ArgumentValidatorFactory::createDoubleGreaterEqualValidator(0.0)).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, approximationHeuristicOptionName, false, "Set the heuristic used for approximation.")
                    .addArgument(storm::settings::ArgumentBuilder::createStringArgument("heuristic", "The name of the heuristic used for approximation.")
                    .setDefaultValueString("depth")
                    .addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator({"depth", "probability", "bounddifference"})).build()).build());
#ifdef STORM_HAVE_Z3
                this->addOption(storm::settings::OptionBuilder(moduleName, solveWithSmtOptionName, true, "Solve the DFT with SMT.").build());
#endif
            }

            bool FaultTreeSettings::useSymmetryReduction() const {
                return this->getOption(symmetryReductionOptionName).getHasOptionBeenSet();
            }

            bool FaultTreeSettings::useModularisation() const {
                return this->getOption(modularisationOptionName).getHasOptionBeenSet();
            }

            bool FaultTreeSettings::isDisableDC() const {
                return this->getOption(disableDCOptionName).getHasOptionBeenSet();
            }

            bool FaultTreeSettings::areRelevantEventsSet() const {
                return this->getOption(relevantEventsOptionName).getHasOptionBeenSet() && (this->getOption(relevantEventsOptionName).getArgumentByName("values").getValueAsString() != "");
            }

            std::vector<std::string> FaultTreeSettings::getRelevantEvents() const {
                return storm::parser::parseCommaSeperatedValues(this->getOption(relevantEventsOptionName).getArgumentByName("values").getValueAsString());
            }

            bool FaultTreeSettings::isApproximationErrorSet() const {
                return this->getOption(approximationErrorOptionName).getHasOptionBeenSet();
            }

            double FaultTreeSettings::getApproximationError() const {
                return this->getOption(approximationErrorOptionName).getArgumentByName("error").getValueAsDouble();
            }

            storm::builder::ApproximationHeuristic FaultTreeSettings::getApproximationHeuristic() const {
                std::string heuristicAsString = this->getOption(approximationHeuristicOptionName).getArgumentByName("heuristic").getValueAsString();
                if (heuristicAsString == "depth") {
                    return storm::builder::ApproximationHeuristic::DEPTH;
                } else if (heuristicAsString == "probability") {
                    return storm::builder::ApproximationHeuristic::PROBABILITY;
                } else if (heuristicAsString == "bounddifference") {
                    return storm::builder::ApproximationHeuristic::BOUNDDIFFERENCE;
                }
                STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Illegal value '" << heuristicAsString << "' set as heuristic for approximation.");
            }

            bool FaultTreeSettings::isTakeFirstDependency() const {
                return this->getOption(firstDependencyOptionName).getHasOptionBeenSet();
            }

#ifdef STORM_HAVE_Z3
            bool FaultTreeSettings::solveWithSMT() const {
                return this->getOption(solveWithSmtOptionName).getHasOptionBeenSet();
            }
#endif

            void FaultTreeSettings::finalize() {
            }

            bool FaultTreeSettings::check() const {
                // Ensure that disableDC and relevantEvents are not set at the same time
                STORM_LOG_THROW(!isDisableDC() || !areRelevantEventsSet(), storm::exceptions::InvalidSettingsException, "DisableDC and relevantSets can not both be set.");
                return true;
            }

        } // namespace modules
    } // namespace settings
} // namespace storm
