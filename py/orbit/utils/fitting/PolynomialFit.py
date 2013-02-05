#-------------------------------------------------------------------------
# This is a polynom fitting class. 
# As input you can use Function or SplineCH instances
# order of the polynomial a0+...+aN*x^N is N
#--------------------------------------------------------------------------
import sys
import math

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import Matrix

from orbit_utils import Polynomial

class PolynomialFit:
	def __init__(self,order):
		self.order = order
		self.polynomial = Polynomial()
		#self.coef_err_arr is a final array with coef_arr and err_arr for polinomial coefficients
		self.coef_err_arr = []
		#self.x_y_err_arr is initial data with (x,y,y_err) points
		self.x_y_err_arr = []
		
	def getPolynomial(self):
		polynomial = Polynomial()
		self.polynomial.copyTo(polynomial)
		return polynomial
	
	def getCoefficientsAndErr(self):
		return self.coef_err_arr
		
	def fitFunction(self, function):
		self.x_y_err_arr = []
		for i in range(function.getSize()):
			x = function.x(i)
			y = function.y(i)
			err = function.err(i)
			self.x_y_err_arr.append([x,y,err])
		self._makePolynomial()
		
	def fitSpline(self, spline):
		self.x_y_err_arr = []
		for i in range(spline.getSize()):
			x = spline.x(i)
			y = spline.y(i)
			err = 0.
			self.x_y_err_arr.append([x,y,err])
		self._makePolynomial()
		
	def _makePolynomial(self):
		nPoints = len(self.x_y_err_arr)
		if(nPoints < (self.order+1)):
			self.order = nPoints - 1
		#check if just one of errors is zero
		infoZeroErr = 1.0
		for [x,y,err] in self.x_y_err_arr:
			infoZeroErr *= err
		for i in range(nPoints):
			[x,y,err] = self.x_y_err_arr[i]
			sigma = 1.0
			if(infoZeroErr != 0.):
				sigma = 1.0/(err*err)
			self.x_y_err_arr[i][2] = sigma
		#now make A matrix
		aMatr = Matrix(nPoints,self.order+1)
		for i in range(nPoints):
			for j in range(self.order+1):
				x = self.x_y_err_arr[i][0]
				aMatr.set(i,j,math.pow(x,j))	
		aTCa = Matrix(self.order+1,self.order+1)
		for i in range(self.order+1):
			for j in range(self.order+1):
				a = 0.
				for k in range(nPoints):
					sigma = self.x_y_err_arr[k][2]
					a += aMatr.get(k,i)*sigma*aMatr.get(k,j)
				aTCa.set(i,j,a)
		#now the resuting coefficients and errors		
		aTCaI = aTCa.invert()
		e = aTCaI.mult(aTCa)		
		if(aTCa == None):
			print "python PolynomialFit: Problem with data."
			for i in range(nPoints):
				x = self.x_y_err_arr[i][0]
				y = self.x_y_err_arr[i][1]
				err = self.x_y_err_arr[i][2]
				print " x,y,err = %12.5g %12.5g %12.5g "%(x,y,err)
			print "Stop."
			sys.exit(1)
		coef_arr = [0.]*(self.order+1)			
		err_arr = [0.]*(self.order+1)
		for i in range(self.order+1):
			err_arr[i] = math.sqrt(math.fabs(aTCaI.get(i,i)))
		for i in range(self.order+1):
			coef_arr[i] = 0.
			for j in range(self.order+1):
				for k in range(nPoints):
					sigma = self.x_y_err_arr[k][2]
					y = self.x_y_err_arr[k][1]
					coef_arr[i] += aTCaI.get(i,j)*aMatr.get(k,j)*sigma*y
		# polinimial coefficients are found
		self.polynomial.order(self.order)
		for i in range(len(coef_arr)):
			self.polynomial.coefficient(i,coef_arr[i])	
		# now let's calculate errors
		if(infoZeroErr == 0.):
			total_sigma = 0.
			for k in range(nPoints):
				x = self.x_y_err_arr[k][0]
				y = self.x_y_err_arr[k][1]
				total_sigma += (self.polynomial.value(x)-y)**2
			total_sigma = math.sqrt(total_sigma/(nPoints-2))
			for i in range(len(err_arr)):
				err_arr[i] *= total_sigma
		# set the resulting coefficients and errors array
		self.coef_err_arr = [coef_arr,err_arr]	
