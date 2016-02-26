#pragma once


#include "DFTElement.h"
namespace storm {
    namespace storage {
        
        template<typename ValueType>
        class DFTDependency : public DFTElement<ValueType> {

            using DFTGatePointer = std::shared_ptr<DFTGate<ValueType>>;
            using DFTBEPointer = std::shared_ptr<DFTBE<ValueType>>;
            
        protected:
            std::string mNameTrigger;
            std::string mNameDependent;
            ValueType mProbability;
            DFTGatePointer mTriggerEvent;
            DFTBEPointer mDependentEvent;

        public:
            DFTDependency(size_t id, std::string const& name, std::string const& trigger, std::string const& dependent, ValueType probability) :
                DFTElement<ValueType>(id, name), mNameTrigger(trigger), mNameDependent(dependent), mProbability(probability)
            {
            }

            virtual ~DFTDependency() {}

            void initialize(DFTGatePointer triggerEvent, DFTBEPointer dependentEvent) {
                assert(triggerEvent->name() == mNameTrigger);
                assert(dependentEvent->name() == mNameDependent);
                mTriggerEvent = triggerEvent;
                mDependentEvent = dependentEvent;
            }

            std::string nameTrigger() const {
                return mNameTrigger;
            }

            std::string nameDependent() const {
                return mNameDependent;
            }

            ValueType const& probability() const {
                return mProbability;
            }

            DFTGatePointer const& triggerEvent() const {
                assert(mTriggerEvent);
                return mTriggerEvent;
            }

            DFTBEPointer const& dependentEvent() const {
                assert(mDependentEvent);
                return mDependentEvent;
            }

            DFTElementType type() const override {
                return DFTElementType::PDEP;
            }

            virtual size_t nrChildren() const override {
                return 1;
            }

            virtual bool isDependency() const override {
                return true;
            }
            
            virtual bool isTypeEqualTo(DFTElement<ValueType> const& other) const override {
                if(!DFTElement<ValueType>::isTypeEqualTo(other)) return false;
                DFTDependency<ValueType> const& otherDEP= static_cast<DFTDependency<ValueType> const&>(other);
                return (mProbability == otherDEP.mProbability);
            }
            

            virtual std::vector<size_t> independentUnit() const override {
                std::set<size_t> unit = {this->mId};
                mDependentEvent->extendUnit(unit);
                if(unit.count(mTriggerEvent->id()) != 0) {
                    return {};
                }
                return std::vector<size_t>(unit.begin(), unit.end());
            }

            virtual void extendSubDft(std::set<size_t>& elemsInSubtree, std::vector<size_t> const& parentsOfSubRoot, bool blockParents) const override {
                 if(elemsInSubtree.count(this->id())) return;
                DFTElement<ValueType>::extendSubDft(elemsInSubtree, parentsOfSubRoot, blockParents);
                if(elemsInSubtree.empty()) {
                    // Parent in the subdft, ie it is *not* a subdft
                    return;
                }
                mDependentEvent->extendSubDft(elemsInSubtree, parentsOfSubRoot, blockParents);
                if(elemsInSubtree.empty()) {
                    // Parent in the subdft, ie it is *not* a subdft
                    return;
                }
                mTriggerEvent->extendSubDft(elemsInSubtree, parentsOfSubRoot, blockParents);
                
                
            }
            
            virtual std::string toString() const override {
                std::stringstream stream;
                bool fdep = storm::utility::isOne(mProbability);
                stream << "{" << this->name() << "} " << (fdep ? "FDEP" : "PDEP") << "(" << mTriggerEvent->name() << " => " << mDependentEvent->name() << ")";
                if (!fdep) {
                    stream << " with probability " << mProbability;
                }
                return stream.str();
            }

        protected:

        };

    }
}