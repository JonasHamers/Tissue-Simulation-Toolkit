from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from tissue_simulation_toolkit.ecm.util import empty_f64, empty_i32
from typing import Optional

@dataclass
class PDEValues:
    """PDE values of the different PDE layers at the positions of the beads

    This data is sent from the Cellular Potts Model to the MD-based ECM model.

    Attributes:
        par_ids: Shape (N, ), particle ID for each row of values.
        values: Shape (N, n_layers), PDE values for each particle.
    """
    par_ids: npt.NDArray[np.int32] = field(default_factory=empty_i32)
    values: npt.NDArray[np.float64] = field(default_factory=empty_f64)
    
    def clear(self) -> None:
        self.par_ids = empty_i32()
        self.values = empty_f64()
        
    def find(self, bead_id: np.int32) -> Optional[npt.NDArray[np.float64]]:
        matches = np.flatnonzero(self.par_ids == bead_id)
        
        if matches.size == 0:
            return None
        
        return self.values[matches[0]]
        
        
        
