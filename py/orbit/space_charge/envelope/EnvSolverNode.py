from orbit.utils import orbitFinalize
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class EnvSolverNode(SC_Base_AccNode):
    
    def __init__(self, sc_calculator, name='no name'):
        """Class implementation of envelope solver node.
    
        sc_calculator : EnvSolverKV or EnvSolverRotating object
            The solver object used to track the bunch.
        """
        SC_Base_AccNode.__init__(self, sc_calculator, name)
        self.setType('envsolver')

    def track(self, paramsDict):
        if not self.switcher:
            return
        bunch = paramsDict['bunch']
        self.sc_calculator.trackBunch(bunch, self.sc_length)
