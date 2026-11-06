#pragma once

#include "vec2.hpp"
#include "ecm_boundary_state.hpp"
#include "cell_ecm_interactions.hpp"

#include <unordered_map>

// Particle ID to its position in the ECM.
struct ECMParticlePositions {
    std::unordered_map<ParId, ParPos> positions;
    
    void update_particle_locations(CellECMInteractions const &interactions);
};
