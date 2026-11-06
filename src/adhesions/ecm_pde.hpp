#pragma once

#include "vec2.hpp"
#include "ecm_boundary_state.hpp"
#include "cell_ecm_interactions.hpp"

#include <unordered_map>

// Particle ID to its position in the ECM.
struct ECMParticlePDEValues {
    std::vector<ParId> par_ids;
    std::vector<std::vector<double>> values;
};
