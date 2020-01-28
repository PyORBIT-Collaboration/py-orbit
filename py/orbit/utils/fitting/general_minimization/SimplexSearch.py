#!/usr/bin/env python
#
# $Id: Simplex.py,v 1.2 2004/05/31 14:01:06 vivake Exp $
# 
# Copyright (c) 2002-2004 Vivake Gupta (vivakeATlab49.com).  All rights reserved.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
#
# This software is maintained by Vivake (vivakeATlab49.com) and is available at:
#     http://shell.lab49.com/~vivake/python/Simplex.py
#
# 1.2  ( 5/2004) - Fixed a bug found by Noboru Yamamoto <noboru.yamamotoATkek.jp>
#                  which caused minimize() not to converge, and reach the maxiter
#                  limit under some conditions.

import math
import sys
import copy

from Solver import SearchAgorithm

#====================================================================
#       class SimplexSearchAlgorithm
#====================================================================

class SimplexSearchAlgorithm(SearchAgorithm):
	""" 
	Simplex - a regression method for arbitrary nonlinear function minimization
	
	Simplex minimizes an arbitrary nonlinear function of N variables by the
	Nedler-Mead Simplex method as described in:
	
	Nedler, J.A. and Mead, R. "A Simplex Method for Function Minimization." 
			Computer Journal 7 (1965): 308-313.
	
	It makes no assumptions about the smoothness of the function being minimized.
	It converges to a local minimum which may or may not be the global minimum
	depending on the initial guess used as a starting point.
	"""
	
	def __init__(self):
		SearchAgorithm.__init__(self)
		"""
		Initializes the simplex.
		_testFunc     the function to minimize
		guess[]       an list containing initial guesses
		increments[]  an list containing increments, perturbation size
		kR            reflection constant
		kE            expansion constant
		kC            contraction constant
		"""
		self.setName("Simplex Search")
		self.guess = None
		self.increments = None
		self.kR = -1.
		self.kE = 2.0
		self.kC = 2.0
		self.numvars = 0
		self.simplex = []
		
		self.lowest = -1
		self.highest = -1
		self.secondhighest = -1
		
		self.errors = []
		self.currenterror = 0.
		
		self.initTrialPoint = None
			
	def setTrialPoint(self,initTrialPoint):
		res = SearchAgorithm.setTrialPoint(self,initTrialPoint)
		if(not res): return res
		
		self.guess = self.initTrialPoint.getVariablesUsedInOptArr()
		self.increments = self.initTrialPoint.getStepsUsedInOptArr()
		
		self.numvars = len(self.guess)
		
		self.errors = []
		self.lowest = -1
		self.highest = -1
		self.secondhighest = -1
		self.currenterror = 0.
		
		# Initialize vertices
		for vertex in range(0, self.numvars + 3): # Two extras to store centroid and reflected point
			self.simplex.append(copy.copy(self.guess))
		# Use initial increments
		for vertex in range(0, self.numvars + 1):
			for x in range(0, self.numvars):
				if x == (vertex - 1):
					self.simplex[vertex][x] = self.guess[x] + self.increments[x]
			self.errors.append(0)
		res = self._calculate_errors_at_vertices()
		if(not res): 
			return False
		return True
		
	def _testFunc(self,guess):
		trialPoint = self.initTrialPoint.getCopy()
		trialPoint.setVariablesUsedInOptArr(guess)	
		trialPoint.setStepsUsedInOptArr(self.increments)
		if(not trialPoint.isAcceptable()):
			return None
		score = self.solver.getScorer().getScore(trialPoint)
		scoreBoard = self.solver.getScoreboard()
		scoreBoard.addScoreTrialPoint(score,trialPoint)
		self.solver.getStopper().checkStopConditions(self.solver)
		if(self.solver.getStopper().getShouldStop()):
			return None
		return score

	def _calculate_errors_at_vertices(self):
		for vertex in range(0, self.numvars + 1):
			if vertex == self.lowest:
				continue
			for x in range(0, self.numvars):
				self.guess[x] = self.simplex[vertex][x]
			val = self._testFunc(self.guess)
			if(val == None): return False
			self.currenterror = val
			self.errors[vertex] = self.currenterror
		return True
		
	def makeStep(self):			
		# Identify highest, second highest, and lowest vertices
		self.highest = 0
		self.lowest = 0
		for vertex in range(0, self.numvars + 1):
			if(self.errors[vertex] > self.errors[self.highest]):
				self.highest = vertex
			if self.errors[vertex] < self.errors[self.lowest]:
				self.lowest = vertex
		self.secondhighest = 0
		for vertex in range(0, self.numvars + 1):
			if vertex == self.highest:
				continue
			if self.errors[vertex] > self.errors[self.secondhighest]:
				self.secondhighest = vertex

		# Calculate centroid of simplex, excluding highest vertex
		for x in range(0, self.numvars):
			S = 0.0
			for vertex in range(0, self.numvars + 1):
				if vertex == self.highest:
					continue
				S = S + self.simplex[vertex][x]
			self.simplex[self.numvars + 1][x] = S / self.numvars
			
		#---- calculate the size of the simplex - that our new increments
		for x in range(0, self.numvars):
			max_deviation = 0.
			for vertex in range(0, self.numvars + 1):
				deviation = abs(self.simplex[self.numvars + 1][x] - self.simplex[vertex][x])
				if(max_deviation < deviation): max_deviation = deviation		
			self.increments[x] = max_deviation
		#-------------------------------------------------------------------
		
		self.reflect_simplex()
		
		#----------------------------------------
		val = self._testFunc(self.guess)
		if(val == None):
			self.solver.getStopper().setShouldStop(True)
			return
		#-----------------------------------------
		self.currenterror = val
		
		if self.currenterror < self.errors[self.lowest]:
			tmp = self.currenterror
			self.expand_simplex()
			#----------------------------------------
			val = self._testFunc(self.guess)
			if(val == None):
				self.solver.getStopper().setShouldStop(True)
				return
			#----------------------------------------
			self.currenterror = val
			if self.currenterror < tmp:
				self.accept_expanded_point()
			else:
				self.currenterror = tmp
				self.accept_reflected_point()

		elif(self.currenterror <= self.errors[self.secondhighest]):
			self.accept_reflected_point()

		elif(self.currenterror <= self.errors[self.highest]):
			self.accept_reflected_point()

			self.contract_simplex()
			#----------------------------------------
			val = self._testFunc(self.guess)
			if(val == None):
				self.solver.getStopper().setShouldStop(True)
				return
			#----------------------------------------			
			self.currenterror = val
			if self.currenterror < self.errors[self.highest]:
				self.accept_contracted_point()
			else:
				res = self.multiple_contract_simplex()
				if(not res):
					self.solver.getStopper().setShouldStop(True)
					return

		elif(self.currenterror >= self.errors[self.highest]):
			self.contract_simplex()
			#----------------------------------------
			val = self._testFunc(self.guess)
			if(val == None):
				self.solver.getStopper().setShouldStop(True)
				return
			#----------------------------------------	
			self.currenterror = val
			if self.currenterror < self.errors[self.highest]:
				self.accept_contracted_point()
			else:
				res = self.multiple_contract_simplex()
				if(not res):
					self.solver.getStopper().setShouldStop(True)
					return				

	def contract_simplex(self):
		for x in range(0, self.numvars):
			self.guess[x] = self.kC * self.simplex[self.highest][x] + (1 - self.kC) * self.simplex[self.numvars + 1][x]
		return

	def expand_simplex(self):
		for x in range(0, self.numvars):
			self.guess[x] = self.kE * self.guess[x]         + (1 - self.kE) * self.simplex[self.numvars + 1][x]
		return

	def reflect_simplex(self):
		for x in range(0, self.numvars):
			self.guess[x] = self.kR * self.simplex[self.highest][x] + (1 - self.kR) * self.simplex[self.numvars + 1][x]
			self.simplex[self.numvars + 2][x] = self.guess[x] # REMEMBER THE REFLECTED POINT
		return

	def multiple_contract_simplex(self):
		for vertex in range(0, self.numvars + 1):
			if vertex == self.lowest:
				continue
			for x in range(0, self.numvars):
				self.simplex[vertex][x] = 0.5 * (self.simplex[vertex][x] + self.simplex[self.lowest][x])
		res = self._calculate_errors_at_vertices()
		if(not res): 
			return False		
		return True

	def accept_contracted_point(self):
		self.errors[self.highest] = self.currenterror
		for x in range(0, self.numvars):
			self.simplex[self.highest][x] = self.guess[x]
		return

	def accept_expanded_point(self):
		self.errors[self.highest] = self.currenterror
		for x in range(0, self.numvars):
			self.simplex[self.highest][x] = self.guess[x]
		return

	def accept_reflected_point(self):
		self.errors[self.highest] = self.currenterror
		for x in range(0, self.numvars):
			self.simplex[self.highest][x] = self.simplex[self.numvars + 2][x]
		return
	