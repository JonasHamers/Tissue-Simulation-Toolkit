from tissue_simulation_toolkit.ecm.muscle3 import (decode_cell_ecm_interactions, decode_mdstate, decode_pde_values, encode_mdstate, encode_ecm_boundary_state, encode_particle_positions, from_settings)
from tissue_simulation_toolkit.ecm.muscle3_mpi_wrapper import Instance
from tissue_simulation_toolkit.ecm.parameters import EvolutionParameters
from tissue_simulation_toolkit.ecm.simulation import Simulation

import hoomd
import hoomd.md
from libmuscle import Message
import numpy as np
import os
from ymmsl import Operator

import logging
from typing import Any, Dict, List


_logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.DEBUG)
    instance = Instance({
        Operator.F_INIT: ['ecm_in'],
        Operator.O_I: ['ecm_boundary_state_out', 'state_out'],
        Operator.S: ['cell_ecm_interactions_in'],
        Operator.O_F: ['ecm_out']})

    while instance.reuse_instance():
        # F_INIT
        par = from_settings(EvolutionParameters, instance)
        
        mcs = instance.get_setting('mcs', 'int')
        try:
            state_output_interval = instance.get_setting('state_output_interval', 'int')
        except KeyError:
            state_output_interval = mcs + 1

        msg = instance.receive('ecm_in')
        ecm = decode_mdstate(msg.data)
        sim = Simulation(par, ecm)

        for i in range(mcs):
            # O_I
            if instance.is_connected('state_out'):
                if i % state_output_interval == 0:
                    state = sim.get_state()
                    instance.send('state_out', Message(i, data=encode_mdstate(state)))
            
            # Old method
            #boundary = encode_ecm_boundary_state(sim.get_boundary_state())
            
            # New method
            ecm_boundary_state = sim.get_boundary_state()
            particle_positions = sim.get_particle_positions()
            
            if ecm_boundary_state is not None:
                data = encode_ecm_boundary_state(ecm_boundary_state)
                
                assert particle_positions is not None
                data["particle_positions"] = encode_particle_positions(particle_positions)
            
                msg = Message(msg.timestamp, data=data)
                instance.send('ecm_boundary_state_out', msg)

            # S
            msg = instance.receive('cell_ecm_interactions_in', default=Message(0.0))
            if msg.data is not None:
                # extract interactions data and pde_values separately
                if 'interactions' in msg.data and msg.data['interactions'] is not None:
                    sim.apply_interactions(decode_cell_ecm_interactions(msg.data['interactions']))
                if 'pde_values' in msg.data and msg.data['pde_values'] is not None:
                    sim._pde_values = decode_pde_values(msg.data['pde_values'])
                
                    # Update the crosslinks
                    sim.create_new_crosslinks(par) # Perhaps set r_max automatically via the parameter file?
                
            sim.run()

        # O_F
        message = Message(msg.timestamp, data=encode_mdstate(sim.get_state()))
        instance.send('ecm_out', message)


if __name__ == '__main__':
    main()
