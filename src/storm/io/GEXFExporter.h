//
// Created by spook on 11.01.24.
//

#pragma once

#include "storm/models/sparse/Mdp.h"

namespace storm {
    namespace exporter {

        template<typename ValueType, typename BeliefType> // TODO are these sufficient / right?
        class GEXFExporter {
        public:
            enum GEXFAttributeType { GEXF_string, GEXF_integer, GEXF_long, GEXF_float, GEXF_double, GEXF_boolean, GEXF_short, GEXF_byte, GEXF_date, GEXF_anyURI };

            GEXFExporter() = default;

        protected:
            uint64_t encodeColor(uint64_t r, uint64_t g, uint64_t b);

            std::string attributeTypeToString(GEXFAttributeType attributeType);

            void exportGEXFToStream(std::shared_ptr<storm::models::sparse::Mdp<ValueType>> mdp, std::ostream& outStream, std::vector<std::vector<uint64_t>> colors, std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>> additionalAttributes = std::map<std::string, std::pair<GEXFAttributeType, std::vector<std::string>>>());

        };

    } // storm
} // exporter

