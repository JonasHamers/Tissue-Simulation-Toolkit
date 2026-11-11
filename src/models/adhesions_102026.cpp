/*

Copyright 1996-2006 Roeland Merks
Copyright 2023 Netherlands eScience Center
Copyright 2006-2026 Tissue Simulation Toolkit contributors

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/

// Standard library headers for input-output operations, 
// memory management, and mathematical functions
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// Platform-specific header for memory allocation functions
#ifndef __APPLE__
#include <malloc.h>
#endif

// Tissue Simulation Toolkit components
#include "cell.hpp"              // cellular potts: cell-related definitions and operations
#include "dish.hpp"              // cellular potts: manages the simulation environment

#include "adhesion_creation.hpp" // adhesions:
#include "cpm_ecm/io.hpp"        // cpm_ecm:
#include "ecm_particle_positions.hpp" // Contains the ECMParticlePositions struct
#include "ecm_pde.hpp"           // Contains the ECMParticlePDEValues struct

#include "graph.hpp"             // graphics: graphics-related functions and definitions
#include "parameter.hpp"         // parameter: manages simulation parameters
#include "plotter.hpp"           // plotter: handles plotting and visualization of the simulation

#include "info.hpp"              // util: handles simulation information and user interactions
#include "profiler.hpp"          // util: profiler functions for performance measurement
#include "random.hpp"            // util: random number generation functions

// Muscle 3 components
#include "util/muscle3/settings.hpp"
#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>

// Chemicals
const int VEGF = 0;
const int LOX = 1;

// Layers
int L_VEGF = -1;
int L_LOX = -1;

// Plotting
int L_plotgrad = -1;
int L_plotcont = -1;


using namespace std;

using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using libmuscle::PortsDescription;
using libmuscle::StorageOrder;
using ymmsl::Operator;

extern Parameter par;

std::unique_ptr<Instance> instance;

/* Implements the function to initialize the simulation state
   INIT is defined as a macro in dish.hpp 
   the code will be added to the Dish::Init method
   dish.hpp declares the pointer CPM to an instance of the 
   CellularPotts class, which is itself defined in ca.hpp.
*/
INIT {
    try {
    	// Define initial distribution of cells
      	CPM -> GrowInCells(par.n_init_cells, par.size_init_cells, par.subfield);
		CPM -> ConstructInitCells(*this);

		// Start with only one big cell and divide it a few times
		// to generate a clump of cells.
		// The behavior can be changed in the parameter file using
		// parameters n_init_cells, size_init_cells, and divisions.
		for (int i = 0; i < par.divisions; i++) {
    		CPM -> DivideCells();
      	}
      
        // Assign layer indices to specific layer variables
        for (int l=0; l<PDEfield->Layers(); l++) {
            switch ((int)par.chem_type[l]) {
                case VEGF: L_VEGF=l; break;
                case LOX: L_LOX=l; break;
            }
        }
        
        // Assign layer indices to plot variables
        for (int l=0; l<PDEfield->Layers(); l++) {
            if (par.cont_color[l]) L_plotcont = l;
            if (par.grad_color[l]) L_plotgrad = l;
        }

		// The edgelist keeps track of pairs of lattice points  
		// that are eligible to change the CPM configuration.
		CPM -> InitialiseEdgeList();

    } catch (const char *error) {
    	cerr << "Caught exception\n";
    	std::cerr << error << "\n";
    	exit(1);
    }
}

TIMESTEP {
	try {

		// Initialize static variables and pointers
		static int i = 0;
		static Dish *dish = new Dish();
		static Info *info = new Info(*dish, *this);
		static Plotter plotter = Plotter(dish, this);
        
        static ECMParticlePositions ecm_particle_positions;
        static ECMParticlePDEValues ecm_particle_pde_values;

		/* cell-ECM model */
		if (i == 0) {
			// request creation of initial adhesions
			CellECMInteractions interactions;

			auto adh_zone = adhesion_zone(*(dish->CPM));

			interactions.change_type_in_area.change_area = adh_zone;
			interactions.change_type_in_area.num_particles = par.num_initial_adhesions;
			interactions.change_type_in_area.from_type = ParticleType::free;
			interactions.change_type_in_area.to_type = ParticleType::adhesion;

			auto data_mem = encode_cell_ecm_interactions(interactions);
            
            auto message_data = Data::dict("interactions", data_mem.first);

			instance -> send("cell_ecm_interactions_out", Message(i, message_data));

		} else {
			// get any adhesion particle movements from CPM and send them out
			auto interactions = dish->CPM->GetCellECMInteractions();
			auto data_mem = encode_cell_ecm_interactions(interactions);
            
            // Use the interactions object to update the particle locations
            ecm_particle_positions.update_particle_locations(interactions);
            
            ecm_particle_pde_values.par_ids.clear();
            ecm_particle_pde_values.values.clear();
            
            // The following block is temporarily here, needs relocation
            // Find the PDE field values and save them
            for (const auto &[par_id, position] : ecm_particle_positions.positions) {
                ecm_particle_pde_values.par_ids.push_back(par_id);
                
                int x = static_cast<int>(std::floor(position.x));
                int y = static_cast<int>(std::floor(position.y));
                
                std::vector<double> row_values;
                for (int layer = 0; layer < dish->PDEfield->Layers(); layer++) {
                    std::cerr << "x = " << x << ", y" << y << std::endl;
                    row_values.push_back(dish->PDEfield->get_PDEvars(layer, x, y));
                }
                
                ecm_particle_pde_values.values.push_back(row_values);
            }
            
            auto data_pde_values_mem = encode_ecm_pde_values(ecm_particle_pde_values);
            
            // Combine the data into one message
            auto message_data = Data::dict("interactions", data_mem.first,
                                           "pde_values", data_pde_values_mem.first);

            // Send the message
			instance -> send("cell_ecm_interactions_out", Message(i, message_data));
    	}

        dish->CPM->ResetCellECMInteractions();

        // Check if the relaxation has finished
        if (i >= par.relaxation) {
            if (par.useopencl) {
                PROFILE(opencl_diff, dish->PDEfield->SecreteAndDiffuseCL(dish->CPM, par.pde_its);)
            } else if (i == par.relaxation) {
                // Initialise the PDE fields
                dish->PDEfield->InitialisePDE(dish->CPM);
                dish->PDEfield->InitialiseDiffusionCoefficients(dish->CPM);

                #ifdef CUDA_ENABLED
                if (par.usecuda) dish->PDEfield->InitialiseCuda();
                #endif

            } else {
                for (int r = 0; r < par.pde_its; r++) {
                    // Run the PDE equations for the preset iterations
                    if (!par.usecuda) {
                        dish->PDEfield->ReactionDiffusion(dish->CPM);
                    }
                    #ifdef CUDA_ENABLED
                    if (par.usecuda) dish->PDEfield->cuPDEsteps(dish->CPM, par.pde_its);
                    #endif
                }
            }
        }

        // Receive the boundary state message from Muscle3
        // Note: after sending a message earlier in TIMESTEP, the code waits here for a
        // message in. This means that the HOOMD code runs between the instance->send() and
        // the instance->receive()
        auto ecm_boundary_state_msg = instance->receive("ecm_boundary_state_in");
        auto ecm_boundary_state = decode_ecm_boundary_state(ecm_boundary_state_msg.data());
        ecm_particle_positions = decode_ecm_particle_positions(ecm_boundary_state_msg.data());
        dish->CPM->SetECMBoundaryState(ecm_boundary_state);
        
        // Here we print some of the particle positions from the message to check that they work
        std::size_t printed = 0;
        for (const auto &[par_id, pos] : ecm_particle_positions.positions) {
            std::cout << "particle " << par_id
                      << " at (" << pos.x << ", " << pos.y << ")\n";
            if (++printed == 5) {
                break;
            }
        }

        // Call AmoebaeMove to perform the CPM timestep
        // Note that dish->PDEfield is a 3-dimensional object, with [layers, sizex, sizey] structure
        PROFILE(amoebamove, dish->CPM->AmoebaeMove(dish->PDEfield);)

        // Send a state_out message if the interval is right
        if (instance->is_connected("state_out")) {
            if (i % instance->get_setting_as<int64_t>("state_output_interval") == 0) {
                std::cerr << "i = " << i << ", sending on state_out" << std::endl;
        
                auto *cpm_sigma = dish->CPM->getSigma()[0];
                Data cpm_state = Data::grid(cpm_sigma, {static_cast<std::size_t>(dish->CPM->SizeX()), static_cast<std::size_t>(dish->CPM->SizeY())}, {"x", "y"}, StorageOrder::last_adjacent);
        
                auto const &pde = dish->PDEfield;
                auto *pde_sigma = pde->getPDEvars()[0][0];
                Data pde_state = Data::grid(pde_sigma, {static_cast<std::size_t>(pde->Layers()), static_cast<std::size_t>(pde->SizeX()), static_cast<std::size_t>(pde->SizeY())}, {"layer", "x", "y"}, StorageOrder::first_adjacent);

                ACT::ActField act_field = dish->CPM->act_field;
                Data ActField2DataDict(ACT::ActField act_field);
                Data act_field_dic = ActField2DataDict(act_field);
                Data state = Data::dict("cpm", cpm_state, "pde", pde_state, "act_field", act_field_dic);
                instance->send("state_out", Message(i, state));
            }
        }

        if (par.graphics && !(i % par.storage_stride)) {
            PROFILE(all_plots, plotter.Plot();)
            char title[400];
            snprintf(title, 399, "CellularPotts: %.2f hr", dish->PDEfield->TheTime() / 3600);
            info->Menu();
        }

        i++;
        
    } catch (const char *error) {
        cerr << "Caught exception\n";
        std::cerr << error << "\n";
        exit(1);
        
    } catch (std::exception const &e) {
        // ensure we crash if there's a problem, Qt swallows exceptions
        std::terminate();
    }
    
    PROFILE_PRINT
}

// Converting Actfield dictionary into a MUSCLE data dictionary   
Data ActField2DataDict(ACT::ActField AF){
  Data AF_dic = Data::dict();
  for (int i = 0; i < par.sizex; i++){
    for (int j = 0; j < par.sizey; j++ ){
      PixelPos pos(i,j);
      // Writing only non-zero values for speed
      if (AF.Value(pos) !=0){
        // assuming a comma seperated format for keys
        // keys should be strings as required by MUSCLE
        std::string str;
        str.append(std::to_string(i));
        str.append(",");
        str.append(std::to_string(j));
        AF_dic[str]= std::int32_t (AF.Value(pos));
        }
    }
  }
  return AF_dic;
}

void PDE::InitialisePDE(CellularPotts *cpm) {
    for (int x = 0; x < sizex; x++) {
        for (int y = 0; y < sizey; y++) {
            for (int l=0; l<Layers(); l++) {
                setValue(l,x,y,par.init_conc[l]);
                //PDEvars[l][x][y] = par.init_conc[l]; // Old method
            }
            //PDEvars[0][x][y] = 0; // Old method for a single layer
        }
    }
    PROFILE_PRINT
}

void PDE::InitialiseDiffusionCoefficients(CellularPotts *cpm) {
    for (int x = 0; x < sizex; x++) {
        for (int y = 0; y < sizey; y++) {
            for (int l = 0; l < par.n_chem; l++) {
                DiffCoeffs[l][x][y] = par.diff_coeff[l];
            }
        }
    }
    PROFILE_PRINT
}

void PDE::DerivativesPDE(CellularPotts *cpm, PDEFIELD_TYPE *derivs, int x, int y) {
    // inside cells
    if (cpm->Sigma(x, y)>0) {
        for (int l=0; l < par.n_chem; l++) {
            derivs[l] = par.secr_rate[l];
        }
        //derivs[0] = par.secr_rate[0]; // Old code for 1 layer
    } else { // outside cells
        for (int l=0; l < par.n_chem; l++) {
            derivs[l] = -par.decay_rate[l] * PDEvars[l][x][y];
        }
        //derivs[0] = -par.decay_rate[0] * PDEvars[0][x][y]; // Old code for 1 layer
    }
    PROFILE_PRINT
}

void PDE::Secrete(CellularPotts *cpm) {
    const double dt = par.dt;
    for (int x = 0; x < sizex; x++) {
        for (int y = 0; y < sizey; y++) {
            // inside cells
            if (cpm->Sigma(x, y)>0) {
                for (int l=0; l < par.n_chem; l++) {
                    PDEvars[l][x][y] = alt_PDEvars[l][x][y] + par.secr_rate[l] * dt;
                }
                //PDEvars[0][x][y] = alt_PDEvars[0][x][y] + par.secr_rate[0] * dt;
            } else {
                // outside cells
                for (int l=0; l < par.n_chem; l++) {
                    PDEvars[l][x][y] = alt_PDEvars[l][x][y] - par.decay_rate[l] * dt * alt_PDEvars[l][x][y];
                }
                //PDEvars[0][x][y] = alt_PDEvars[0][x][y] - par.decay_rate[0] * dt * alt_PDEvars[0][x][y];
            }
        }
    }
    PROFILE_PRINT
}

int PDE::MapColour(double val) {
  return (((int)((val / ((val) + 1.)) * 100)) % 100) + 155;
}

void Plotter::Plot() {
    graphics->BeginScene();
    graphics->ClearImage();

    if (L_plotgrad > -1) {
        SetPDEDensityLayer(L_plotgrad);
        plotPDEDensity();
    }
    if (L_plotgrad > -1) {
        SetPDEContourLayer(L_plotcont);
        plotPDEContourLines();
    }
    
    //plotPDEDensity();
    plotCPMCellTypes();
    plotCPMLines();
    //plotPDEContourLines();

    graphics->EndScene();
}

int main(int argc, char *argv[]) {

  /* muscle3 */
  PortsDescription ports(
      {{Operator::O_I, {"cell_ecm_interactions_out", "state_out"}},
       {Operator::S, {"ecm_boundary_state_in"}}});
  instance = std::make_unique<Instance>(0, nullptr, ports);

  instance->reuse_instance();
  set_parameters_from_settings(*instance);
  Seed(par.rseed);

  try {
    start_graphics(argc, argv);
  } catch (const char *error) {
    std::cerr << error << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "An unknown exception was caught" << std::endl;
    return 1;
  }

  // This is a hack, the whole model is really supposed to be inside a while
  // loop guarded by this statement. The architecture here won't allow that
  // and fortunately we don't need to actually run more than once, but
  // the function still needs to be called here for MUSCLE3 to shut down
  // correctly.
  instance->reuse_instance();

  return 0;
}
