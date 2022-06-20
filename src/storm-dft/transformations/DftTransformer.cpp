#include "DftTransformer.h"

#include "storm/exceptions/UnexpectedException.h"
#include "storm/utility/macros.h"

#include "storm-dft/builder/DFTBuilder.h"

namespace storm::dft {
namespace transformations {

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftTransformer<ValueType>::transformUniqueFailedBE(storm::dft::storage::DFT<ValueType> const &dft) {
    STORM_LOG_DEBUG("Start transformation UniqueFailedBE");
    storm::dft::builder::DFTBuilder<ValueType> builder;
    // NOTE: if probabilities for constant BEs are introduced, change this to vector of tuples (name, prob)
    std::vector<std::string> failedBEs;

    for (size_t i = 0; i < dft.nrElements(); ++i) {
        std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> element = dft.getElement(i);
        switch (element->type()) {
            case storm::dft::storage::elements::DFTElementType::BE: {
                auto be = std::static_pointer_cast<storm::dft::storage::elements::DFTBE<ValueType> const>(element);
                switch (be->beType()) {
                    case storm::dft::storage::elements::BEType::CONSTANT: {
                        // Remember constant failed BEs for later
                        auto beConst = std::static_pointer_cast<storm::dft::storage::elements::BEConst<ValueType> const>(element);
                        if (beConst->canFail()) {
                            STORM_LOG_TRACE("Transform " << beConst);
                            failedBEs.push_back(beConst->name());
                        }
                        // All original constant BEs are set to failsafe, failed BEs are later triggered by a new element
                        builder.addBasicElementConst(beConst->name(), false);
                        break;
                    }
                    case storm::dft::storage::elements::BEType::PROBABILITY: {
                        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                                        "BE with constant probability distribution are not supported and need to be transformed before.");
                        break;
                    }
                    default:
                        // Clone other types of BEs
                        builder.cloneElement(element);
                        break;
                }
                break;
            }
            default:
                // Clone other elements
                builder.cloneElement(element);
                break;
        }
    }
    // At this point the DFT is an exact copy of the original, except for all constant failure probabilities being 0

    // Introduce new constantly failed BE and FDEPs to trigger all failures
    if (!failedBEs.empty()) {
        STORM_LOG_TRACE("Add unique constant failed BE 'Unique_Constant_Failure'");
        builder.addBasicElementConst("Unique_Constant_Failure", true);
        failedBEs.insert(failedBEs.begin(), "Unique_Constant_Failure");
        STORM_LOG_TRACE("Add FDEP 'Failure_Trigger'");
        builder.addPdep("Failure_Trigger", failedBEs, storm::utility::one<ValueType>());
    }

    builder.setTopLevel(dft.getTopLevelElement()->name());

    STORM_LOG_DEBUG("Transformation UniqueFailedBE complete");
    return std::make_shared<storm::dft::storage::DFT<ValueType>>(builder.build());
}

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftTransformer<ValueType>::transformBinaryDependencies(storm::dft::storage::DFT<ValueType> const &dft) {
    STORM_LOG_DEBUG("Start transformation BinaryDependencies");
    storm::dft::builder::DFTBuilder<ValueType> builder;

    for (size_t i = 0; i < dft.nrElements(); ++i) {
        std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> element = dft.getElement(i);
        switch (element->type()) {
            case storm::dft::storage::elements::DFTElementType::PDEP: {
                auto dep = std::static_pointer_cast<storm::dft::storage::elements::DFTDependency<ValueType> const>(element);
                if (dep->dependentEvents().size() == 1) {
                    // Already binary dependency -> simply clone element
                    builder.cloneElement(dep);
                } else {
                    if (!storm::utility::isOne(dep->probability())) {
                        // PDEP with probability < 1
                        STORM_LOG_TRACE("Transform " << element);
                        // Introduce additional element to first capture the probabilistic dependency
                        std::string nameAdditional = dep->name() + "_additional";
                        STORM_LOG_TRACE("Add auxiliary BE " << nameAdditional);
                        builder.addBasicElementConst(nameAdditional, false);
                        STORM_LOG_TRACE("Add PDEP " << dep->name() << "_pdep");
                        // First consider probabilistic dependency
                        builder.addPdep(dep->name() + "_pdep", {dep->triggerEvent()->name(), nameAdditional}, dep->probability());
                        // Then consider dependencies to the children if probabilistic dependency failed
                        for (size_t j = 0; j < dep->dependentEvents().size(); ++j) {
                            std::string nameDep = dep->name() + "_" + std::to_string(j);
                            std::string dependentName = dep->dependentEvents()[j]->name();
                            STORM_LOG_TRACE("Add FDEP " << nameDep << " for " << dependentName);
                            builder.addPdep(nameDep, {nameAdditional, dependentName}, storm::utility::one<ValueType>());
                        }
                    } else {
                        // FDEP -> add explicit dependencies for each dependent event
                        STORM_LOG_TRACE("Transform " << element);
                        for (size_t j = 0; j < dep->dependentEvents().size(); ++j) {
                            std::string nameDep = dep->name() + "_" + std::to_string(j);
                            std::string dependentName = dep->dependentEvents()[j]->name();
                            STORM_LOG_TRACE("Add FDEP " << nameDep << " for " << dependentName);
                            builder.addPdep(nameDep, {dep->triggerEvent()->name(), dependentName}, storm::utility::one<ValueType>());
                        }
                    }
                }
                break;
            }
            default:
                // Clone other elements
                builder.cloneElement(element);
                break;
        }
    }

    builder.setTopLevel(dft.getTopLevelElement()->name());

    STORM_LOG_DEBUG("Transformation BinaryDependencies complete");
    return std::make_shared<storm::dft::storage::DFT<ValueType>>(builder.build());
}

// Explicitly instantiate the class.
template class DftTransformer<double>;
template class DftTransformer<RationalFunction>;

}  // namespace transformations
}  // namespace storm::dft