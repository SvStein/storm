// Include other headers.
#include <chrono>
#include "src/exceptions/BaseException.h"
#include "src/utility/macros.h"
#include "src/cli/cli.h"
#include "src/utility/initialize.h"

#include "src/settings/SettingsManager.h"
#include "src/settings/modules/GeneralSettings.h"

/*!
 * Main entry point of the executable storm.
 */
int main(const int argc, const char** argv) {

    try {
        auto start = std::chrono::high_resolution_clock::now();
        storm::utility::setUp();
        storm::cli::printHeader("Storm", argc, argv);
        storm::settings::initializeAll("Storm", "storm");
        bool optionsCorrect = storm::cli::parseOptions(argc, argv);
        if (!optionsCorrect) {
            return -1;
        }
        
        // From this point on we are ready to carry out the actual computations.
        storm::cli::processOptions();
        
        // All operations have now been performed, so we clean up everything and terminate.
        storm::utility::cleanUp();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        if(storm::settings::getModule<storm::settings::modules::GeneralSettings>().isPrintTimingsSet()) {
            std::cout << "Overal runtime: " << duration.count() << " ms. (approximately " << durationSec.count() << " seconds)." << std::endl;
        }
        return 0;
    } catch (storm::exceptions::BaseException const& exception) {
        STORM_LOG_ERROR("An exception caused Storm to terminate. The message of the exception is: " << exception.what());
    } catch (std::exception const& exception) {
        STORM_LOG_ERROR("An unexpected exception occurred and caused Storm to terminate. The message of this exception is: " << exception.what());
    }
}
