## \namespace orbit::utils::fitting
##
## Classes:
## - PolynomialFit - fitting a Function or SplineCH instances with the plynomial


from orbit.utils.fitting.PolynomialFit import PolynomialFit

from orbit.utils.fitting.general_minimization.SimplexSearch import SimplexSearchAlgorithm
from orbit.utils.fitting.general_minimization.GoldenSectionSearch1D import GoldenSectionSearchAlgorithm
from orbit.utils.fitting.general_minimization.BisectionSearch1D import BisectionSearchAlgorithm

from orbit.utils.fitting.general_minimization.Solver import Solver
from orbit.utils.fitting.general_minimization.Solver import TrialPoint
from orbit.utils.fitting.general_minimization.Solver import SolveStopper
from orbit.utils.fitting.general_minimization.Solver import SolveStopperFactory
from orbit.utils.fitting.general_minimization.Solver import ScoreboardActionListener
from orbit.utils.fitting.general_minimization.Solver import VariableProxy
from orbit.utils.fitting.general_minimization.Solver import Scorer
from orbit.utils.fitting.general_minimization.Solver import SearchAgorithm 

__all__ = []
__all__.append("PolynomialFit")
__all__.append("SimplexSearchAlgorithm")
__all__.append("GoldenSectionSearchAlgorithm")
__all__.append("BisectionSearchAlgorithm")
__all__.append("Solver")
__all__.append("TrialPoint")
__all__.append("SolveStopper")
__all__.append("SolveStopperFactory")
__all__.append("ScoreboardActionListener")
__all__.append("VariableProxy")
__all__.append("Scorer")
__all__.append("SearchAgorithm")


