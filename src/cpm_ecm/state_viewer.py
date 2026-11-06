"""Model state viewer

Based on artist.py in TST-MD
"""
import logging
from pathlib import Path
from typing import List

from libmuscle import Instance
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from ymmsl import Operator

from tissue_simulation_toolkit.ecm.ecm import ParticleType
from tissue_simulation_toolkit.cpm_ecm.state_plotter import StatePlotter


_logger = logging.getLogger(__name__)


def main() -> None:
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger('PIL').setLevel(logging.INFO)
    logging.getLogger('matplotlib').setLevel(logging.INFO)

    instance = Instance({
        Operator.S: ['cpm_state_in', 'ecm_state_in']})
        
    # Set the output directory
    run_dir = Path.cwd()
    while run_dir != run_dir.parent and not run_dir.name.startswith("run_adhesions_"):
        run_dir = run_dir.parent
    
    if not run_dir.name.startswith("run_adhesions_"):
        raise RuntimeError("Coudl not find the MUSCLE run directory")
        
    out_dir = run_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    while instance.reuse_instance():
        sizex = instance.get_setting('sizex', 'int')
        sizey = instance.get_setting('sizey', 'int')
        image_height = instance.get_setting('image_height', 'int')
        plotter = StatePlotter(sizex, sizey, image_height)

        state_output_interval = instance.get_setting('state_output_interval', 'int')
        mcs = instance.get_setting('mcs', 'int')

        for i in range(mcs):
            if i % state_output_interval == 0:
                cpm_state_msg = instance.receive('cpm_state_in')
                ecm_state_msg = instance.receive('ecm_state_in')

                particles = ecm_state_msg.data['particles']
                
                # Get bond data
                bonds = ecm_state_msg.data['bonds']

                plotter.draw(
                        i,
                        particles['positions'].array,
                        particles['types'].array,
                        ecm_state_msg.data['bonds']['groups'].array,
                        bonds['types'].array,
                        cpm_state_msg.data['pde'].array,
                        cpm_state_msg.data['cpm'].array,
                        cpm_state_msg.data['act_field'],
                        save=True,
                        out_dir=out_dir)
                        
        print(cpm_state_msg.data['pde'].array)

if __name__ == '__main__':
    main()
