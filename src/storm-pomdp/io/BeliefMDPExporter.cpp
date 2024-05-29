#include "BeliefMDPExporter.h"
#include <boost/container/flat_map.hpp>
#include "storm/io/file.h"
#include "storm/adapters/RationalFunctionAdapter.h"


namespace storm{
namespace exporter{
template<class ValueType, typename BeliefType>
void BeliefMDPExporter<ValueType, BeliefType>::determineOgColors(std::vector<std::vector<uint64_t>> &stateColors) {
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
    if (numberOfColors > 1530) {
        determineOgColorsExtended(stateColors);
        return;
    }
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
void BeliefMDPExporter<ValueType, BeliefType>::determineOgColorsExtended(std::vector<std::vector<uint64_t>> &stateColors) {
    auto numberOfColors = stateColors.size();
    // We wanna divide the color "cube" with:
    // x "cuts" orthogonal to r axis,
    // y "cuts" orthogonal to g axis,
    // z "cuts" orthogonal to b axis,
    // and take the corner points of the resulting sections (minus black and white bc they are reserved) as our different color values
    // we check value triples with x >= y >= z >= x-1
    // with x,y,z close to each other for a relatively even distribution across the color space
    // but allowing that y and z may be 1 smaller than x to provide intermediate sizes for the color space
    if (256 * 256 * 256 - 2 < numberOfColors) {
        STORM_LOG_ERROR("Color Space too small for required number of colors");
    }
    uint64_t x = 1;
    uint64_t y;
    uint64_t z;
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
        stateColors[state][0] = std::get<0>(rgb) * (255/(x + 1));
        stateColors[state][1] = std::get<1>(rgb) * (255/(y + 1));
        stateColors[state][2] = std::get<2>(rgb) * (255/(z + 1));
    }
}

template<typename ValueType, typename BeliefType>
std::tuple<uint64_t, uint64_t, uint64_t> BeliefMDPExporter<ValueType, BeliefType>::adaptedEuclid(uint64_t i, uint64_t x, uint64_t y) {
    uint64_t b = i / ((x + 2) * (y + 2));
    i = i - b * (x + 2) * (y + 2);
    uint64_t g = i / (x + 2);
    i = i - g * (x + 2);
    uint64_t r = i;
    return std::make_tuple(r, g, b);
}

template<class ValueType, typename BeliefType>
void BeliefMDPExporter<ValueType, BeliefType>::determineNumberOfEpochs(std::vector<std::string> &numbersOfEpochs, std::vector<std::string> &maxSingleStateEpochNumbers, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox) {
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
void BeliefMDPExporter<ValueType, BeliefType>::determineUnfColors(std::vector<std::vector<uint64_t>> &ogStateColors, std::vector<std::vector<uint64_t>> &unfStateColors, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo,  bool underApprox) {
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

template<typename ValueType, typename BeliefType>
void BeliefMDPExporter<ValueType, BeliefType>::outputSummary(std::vector<std::string> numEpochs, std::vector<std::string> maxSingleStateEpochs, bool under, std::string filename) {
    auto stringIntComp = [] (const std::string &v1, const std::string &v2) -> bool {
        return std::stoi(v1) < std::stoi(v2);
    };
    std::string epochMax = *(std::max_element(numEpochs.begin(), numEpochs.end(), stringIntComp));
    std::string maxSingleStateEpochsMax = *(std::max_element(maxSingleStateEpochs.begin(), maxSingleStateEpochs.end(), stringIntComp));
    std::string fileNameWithoutPath = filename.substr(filename.find_last_of('/'));
    std::ofstream stream;
    storm::utility::openFile(summaryFile, stream, true);
    std::string lineToWrite = epochMax + "\t" + maxSingleStateEpochsMax + "\t" + fileNameWithoutPath + (under? "_under" : "_over") + "\n";
    stream << lineToWrite;
    storm::utility::closeFile(stream);
}

template<class ValueType, typename BeliefType>
void BeliefMDPExporter<ValueType, BeliefType>::createGEXFOutputs(typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &ogCheckingResult, typename storm::pomdp::modelchecker::BeliefExplorationPomdpModelChecker<storm::models::sparse::Pomdp<ValueType>, BeliefType>::Result &unfCheckingResult, typename storm::transformer::BoundUnfolder<ValueType>::UnfoldingResult unfoldingInfo, std::string filename) {
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
    auto unfUnderExtraAttr = std::map<std::string, std::pair<typename GEXFExporter<ValueType, BeliefType>::GEXFAttributeType, std::vector<std::string>>>();
    unfUnderExtraAttr["numberOfEpochs"] = std::make_pair(GEXFExporter<ValueType, BeliefType>::GEXFAttributeType::GEXF_integer, unfUnderNumberOfEpochs);
    unfUnderExtraAttr["maxSingleStateEpochNumber"] = std::make_pair(GEXFExporter<ValueType, BeliefType>::GEXFAttributeType::GEXF_integer, unfUnderMaxSingleStateEpochNumbers);

    std::ofstream stream;
    storm::utility::openFile(filename + "_ogUnder.gexf", stream);
    this->exportGEXFToStream(ogCheckingResult.beliefMdpUnder, stream, ogUnderColors);
    storm::utility::closeFile(stream);
    stream.clear();
    storm::utility::openFile(filename + "_unfUnder.gexf", stream);
    this->exportGEXFToStream(unfCheckingResult.beliefMdpUnder, stream, unfUnderColors, unfUnderExtraAttr);
    storm::utility::closeFile(stream);
    stream.clear();

    outputSummary(unfUnderNumberOfEpochs, unfUnderMaxSingleStateEpochNumbers, true, filename);

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
    auto unfOverExtraAttr = std::map<std::string, std::pair<typename GEXFExporter<ValueType, BeliefType>::GEXFAttributeType, std::vector<std::string>>>();
    unfOverExtraAttr["numberOfEpochs"] = std::make_pair(GEXFExporter<ValueType, BeliefType>::GEXFAttributeType::GEXF_integer, unfOverNumberOfEpochs);
    unfOverExtraAttr["maxSingleStateEpochNumber"] = std::make_pair(GEXFExporter<ValueType, BeliefType>::GEXFAttributeType::GEXF_integer, unfOverMaxSingleStateEpochNumbers);

    storm::utility::openFile(filename + "_ogOver.gexf", stream);
    this->exportGEXFToStream(ogCheckingResult.beliefMdpOver, stream, ogOverColors);
    storm::utility::closeFile(stream);
    stream.clear();
    storm::utility::openFile(filename + "_unfOver.gexf", stream);
    this->exportGEXFToStream(unfCheckingResult.beliefMdpOver, stream, unfOverColors, unfOverExtraAttr);
    storm::utility::closeFile(stream);

    outputSummary(unfOverNumberOfEpochs, unfOverMaxSingleStateEpochNumbers, false, filename);
}

template class BeliefMDPExporter<double, double>;
template class BeliefMDPExporter<double, storm::RationalNumber>;
template class BeliefMDPExporter<storm::RationalNumber, double>;
template class BeliefMDPExporter<storm::RationalNumber, storm::RationalNumber>;

}
}