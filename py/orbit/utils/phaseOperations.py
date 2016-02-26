"""
The collection of functions and classes for operations with phases.
They utilized the 360 deg or 2*math.pi periodicity of the phases.
"""

import math

def phaseNearTargetPhase(phase,phase_trgt):
	"""
	Adds or subtracts 2*math.pi to get the phase near the target phase.
	"""
	pi2 = 2*math.pi
	delta = pi2*int((phase_trgt - phase)/(pi2))
	phase += delta 
	if(phase_trgt - phase > math.pi): 
		phase += pi2
		return phase
	if(phase_trgt - phase < -math.pi):
		phase -= pi2
	return phase
	
def phaseNearTargetPhaseDeg(phase,phase_trgt):
	"""
	Adds or subtracts 360 to get the phase near the target phase.
	"""
	delta = 360.*int((phase_trgt - phase)/(360.))
	phase += delta 
	if(phase_trgt - phase > 180.): 
		phase += 360.
		return phase
	if(phase_trgt - phase < -180.):
		phase -= 360.
	return phase