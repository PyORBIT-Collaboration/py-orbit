
from orbit.space_charge.envelope import EnvSolverNode
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes

def setEnvAccNodes(lattice, path_length_min, sc_calculator):
    """Place set of envelope solver nodes into the lattice.
    
    The method will place the set into the lattice as child nodes of
    the first level accelerator nodes. The nodes will be inserted at
    the beginning of a particular part of the first level AccNode
    element.
    
    Parameters
    ----------
    lattice : AccLattice object
        The lattice in which to insert the nodes.
    path_length_min : float
        The minimum distance between the nodes.
    sc_calculator : EnvSolverKV or EnvSolverRotating object
        The solver object used to track the bunch.
        
    Returns
    -------
    list[EnvSolverNode]
        The list of inserted envelope solver nodes.
    """
    sc_nodes = setSC_General_AccNodes(lattice, path_length_min, sc_calculator, EnvSolverNode)
    for sc_node in sc_nodes:
        sc_node.setName(sc_node.getName() + 'envsolver')
    lattice.initialize()
    return sc_nodes
