#include "ecm_particle_positions.hpp"
#include "cell_ecm_interactions.hpp"

void ECMParticlePositions::update_particle_locations(CellECMInteractions const &interactions) {
    // Access the adhesion particles movements
    const auto &moves = interactions.move_adhesion_particles;
    
    // Loop through all the movements
    for (std::size_t i = 0; i < moves.par_id.size(); i++) {
        const ParId id = moves.par_id[i];
        const ParPos new_pos = moves.new_pos[i];
        
        auto it = positions.find(id);
        it->second = new_pos;
        
        std::cout << "Particle has been moved";
    }
}
