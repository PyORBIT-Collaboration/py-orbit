import sys
import math

import orbit_mpi
from orbit_mpi import mpi_comm

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator
from orbit.utils.fitting import PolynomialFit
from orbit_utils import Polynomial

class RF_AxisFieldAnalysis:
	"""
	This class analyzes the RF electric field on the axis of a whole cavity.
	The result of analysis are Time Transit Factors T,Tp,S,Sp for all gaps
	found in the cavity.
	"""
	def __init__(self,splineFiled):
		self.splineFiled = splineFiled
		#----------------------------------------------------
		self.eps_root = 1.0e-6
		self.rf_freq = 0.0
		#self.e0_normalized_arr - normilized amplitudes of the gaps
		self.e0_normalized_arr = []
		self.e0l_normalized_arr = []
		# self.beta_arr - relativistic beta, cappa = 2*math.pi*rf_freq/(c_light*beta)
		self.beta_arr = []
		self.cappa_arr = []	
		# self.ttp_ssp_gap_arr array with [T,Tp,S,Sp] for all gaps (T,Tp,S,Sp - Functions of cappa)
		self.ttp_ssp_gap_arr = []
		# self.gap_polynoms_coef_arr array with polynomial coefficients for [T,Tp,S,Sp] for each gap
		self.gap_polynoms_coef_arr = []
		# self.gap_polynoms_arr array with Polynomial instances for [T,Tp,S,Sp] for each gap
		self.gap_polynoms_arr = []		
		# self.gap_polynoms_t_tp_s_sp_err_arr - maximal relative errors for polynomial fitting
		self.gap_polynoms_t_tp_s_sp_err_arr = []
		#-----------------------------------------------------
		#calculate the roots
		self.roots_arr = self.rootAnalysis()
		#find the roots of derivative - yp = y' - RFgap center positions
		self.yp_roots_arr = self.gapCentersAnalysis()
		#print "debug yp roots=",self.yp_roots_arr
		if(len(self.roots_arr) - 1 != len(self.yp_roots_arr)):
			rank = orbit_mpi.MPI_Comm_rank(mpi_comm.MPI_COMM_WORLD)
			if(rank == 0):
				print "Class RF_AxisFieldAnalysis."
				print "The structure of the input rf field spline is wrong!"
				print "roots of the filed =",self.roots_arr 
				print "extrema positions =",self.yp_roots_arr
			sys.exit(1)
		# caluclate the position of the center of the cavity 
		rf_center = 0
		for i in range(1,len(self.yp_roots_arr)-1):
			rf_center += self.yp_roots_arr[i]
		rf_center /= (len(self.yp_roots_arr)-2)
		self.rf_center = rf_center
		# make spline for each RF gap
		self.gap_slpline_arr = []
		#print "debug roots_arr=",self.roots_arr
		#---make splineGap with x in the [m] instead of [cm]
		for i in range(len(self.roots_arr)-1):			
			x_center = self.yp_roots_arr[i]
			x0 = self.roots_arr[i]
			x1 = self.roots_arr[i+1]
			f = Function()
			f.add((x0-x_center),math.fabs(splineFiled.getY(x0)))
			for ix in range(splineFiled.getSize()-1):
				x = splineFiled.x(ix)
				if(x > x0 and x < x1):
					f.add((x-x_center),math.fabs(splineFiled.y(ix)))
			f.add((x1-x_center),math.fabs(splineFiled.getY(x1)))
			splineGap = SplineCH()					
			splineGap.compile(f)
			n = splineGap.getSize()
			x_min = splineGap.x(0)
			x_max = splineGap.x(n-1)
			gap_length = x_max - x_min
			self.gap_slpline_arr.append([gap_length,(x_center - self.rf_center),splineGap])

	def rootAnalysis(self):
		""" 
		The method will find the roots of the field distribution, and
		it will return the root array including the left and right edge points
		"""
		#calculate the pairs of roots
		roots_arr = []
		roots_arr.append(self.splineFiled.x(0))
		iz = 1
		while( iz < (self.splineFiled.getSize()-2)):
			x0 = self.splineFiled.x(iz)
			x1 = self.splineFiled.x(iz+1)
			y0 = self.splineFiled.y(iz)
			y1 = self.splineFiled.y(iz+1)
			if(y1*y0 <= 0.):
				if(y0*y1 == 0.):
					if(y0 == 0.):
						roots_arr.append(x0)
					else:
						roots_arr.append(x1)
						iz = iz + 1
				else:
					while(math.fabs(x0 - x1) > self.eps_root):
						x = (x0 + x1)/2.0
						if(self.splineFiled.getY(x) == 0.):
							roots_arr.append(x)
							break
						else:
							if(self.splineFiled.getY(x)*self.splineFiled.getY(x0) < 0.):
								x1 = x
							else:
								x0 = x
					roots_arr.append((x1+x0)/2.0)
			iz = iz + 1
		roots_arr.append(self.splineFiled.x(self.splineFiled.getSize()-1))
		for i in range(self.splineFiled.getSize()):
			x = self.splineFiled.x(i)
			y = self.splineFiled.y(i)
			yp = self.splineFiled.getYP(x)
			#print "debug i=",i," x= %8.4f y = %12.5e   yp= %12.5e "%(x,y,yp)
		return roots_arr
		
	def getNormilizedSpline(self):
		"""
		Returns the spline  normilized by the integral of the absolute value.
		"""
		n = self.splineFiled.getSize()
		f = Function()
		for i in range(n):
			f.add(self.splineFiled.x(i),math.fabs(self.splineFiled.y(i)))
		integral = GaussLegendreIntegrator(500)
		integral.setLimits(self.splineFiled.x(0),self.splineFiled.x(n-1))
		spline = SplineCH()					
		spline.compile(f)
		res = integral.integral(spline)	
		f = Function()
		for i in range(n):
			f.add(self.splineFiled.x(i),self.splineFiled.y(i)/res)		
		spline = SplineCH()					
		spline.compile(f)
		return spline		
		
	def gapCentersAnalysis(self):
		"""
		The method will calculate the positions of d(Ez)/dz = 0 as centers of the gaps
		"""
		yp_roots_arr = []
		iz = 0
		while( iz < (self.splineFiled.getSize()-1)):
			x0 = self.splineFiled.x(iz)
			x1 = self.splineFiled.x(iz+1)
			y0 = self.splineFiled.getYP(x0)
			y1 = self.splineFiled.getYP(x1)
			if(y1*y0 <= 0.):
				if(y0*y1 == 0.):
					if(y0 == 0.):
						yp_roots_arr.append(x0)
					else:
						yp_roots_arr.append(x1)
						iz = iz + 1
				else:
					while(math.fabs(x0 - x1) > self.eps_root):
						x = (x0 + x1)/2.0
						if(self.splineFiled.getYP(x) == 0.):
							yp_roots_arr.append(x)
							break
						else:
							if(self.splineFiled.getYP(x)*self.splineFiled.getYP(x0) < 0.):
								x1 = x
							else:
								x0 = x
					yp_roots_arr.append((x1+x0)/2.0)
			iz = iz + 1
		return yp_roots_arr
		
	def makeTransitTimeTables(self,beta_min,beta_max,n_table_points,rf_freq):
		"""
		It will calculate transit time factor tables for all RF gaps
		TTFs (T,S,Tp,Sp) are funcftions of the cappa variable = 2*pi*f/(c*beta)
		"""
		self.rf_freq = rf_freq
		c_light = 2.99792458e+8
		self.beta_arr = []
		self.cappa_arr = []
		for i_beta in range(n_table_points):
			beta = beta_min + i_beta*(beta_max-beta_min)/(n_table_points-1)
			cappa = 2*math.pi*rf_freq/(c_light*beta)
			self.beta_arr.append(beta)
			self.cappa_arr.append(cappa)
		self.beta_arr.reverse()
		self.cappa_arr.reverse()
		#--calculate realtive gap amplitudes
		integral = GaussLegendreIntegrator(500)
		e0l_arr = []
		e0l_sum = 0.
		for i in range(len(self.gap_slpline_arr)):
			[gap_length,x_center,splineGap] = self.gap_slpline_arr[i]
			n = splineGap.getSize()
			x_min = splineGap.x(0)
			x_max = splineGap.x(n-1)
			integral.setLimits(x_min,x_max)
			e0l = integral.integral(splineGap)
			e0l_sum += e0l
			e0l_arr.append(e0l)
		self.e0_normalized_arr = []
		self.e0l_normalized_arr = []
		e0_norm = e0l_arr[0]/self.gap_slpline_arr[0][0]
		e0l_norm = e0l_arr[0]
		for i in range(len(e0l_arr)):
			self.e0_normalized_arr.append((e0l_arr[i]/self.gap_slpline_arr[i][0])/e0_norm)
			self.e0l_normalized_arr.append((e0l_arr[i]/e0l_norm))
		#--- calculate transit time factors
		self.ttp_ssp_gap_arr = []
		for i in range(len(self.gap_slpline_arr)):
			func_T  = Function()
			func_TP = Function()
			func_S  = Function()
			func_SP = Function()
			self.ttp_ssp_gap_arr.append([func_T,func_TP,func_S,func_SP])
		for i_gap in range(len(self.gap_slpline_arr)):
			[func_T,func_TP,func_S,func_SP] = self.ttp_ssp_gap_arr[i_gap]
			[gap_length,x0,spline] = self.gap_slpline_arr[i_gap]
			x_min = spline.x(0)
			x_max = spline.x(spline.getSize()-1)
			integral.setLimits(x_min,x_max)		
			for i_beta in range(n_table_points):
				beta = self.beta_arr[i_beta]
				cappa = self.cappa_arr[i_beta]
				f_cos = Function()
				f_sin = Function()	
				for isp in range(spline.getSize()):
					x = spline.x(isp)
					y = spline.y(isp)
					phase = cappa*x
					s = math.sin(phase)
					c = math.cos(phase)
					f_cos.add(x,c*y)
					f_sin.add(x,s*y)
				f_sp_cos = SplineCH()	
				f_sp_sin = SplineCH()	
				f_sp_cos.compile(f_cos)	
				f_sp_sin.compile(f_sin)	
				T = integral.integral(f_sp_cos)	
				S = integral.integral(f_sp_sin)
				func_T.add(cappa,T/e0l_arr[i_gap])
				func_S.add(cappa,S/e0l_arr[i_gap])	
			spline_T = SplineCH()
			spline_S = SplineCH()
			spline_T.compile(func_T)
			spline_S.compile(func_S)
			for i_beta in range(spline_T.getSize()):
				cappa = spline_T.x(i_beta)
				TP = spline_T.getYP(cappa)
				SP = spline_S.getYP(cappa)
				func_TP.add(cappa,TP)
				func_SP.add(cappa,SP)
		return self.ttp_ssp_gap_arr

	def makePlynomialFittings(self,n_order):
		"""
		The method will prepare the polynomial fitting for the T,Tp,S,Sp 
		as functions of cappa for all RF gaps
		TTFs (T,S,Tp,Sp) are funcftions of the cappa variable = 2*pi*f/(c*beta)
		"""
		if(len(self.ttp_ssp_gap_arr) == 0):
			print "Please, call makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq) first!"
			print "Stop."
			sys.exit(1)
			return		
		polynomialFit = PolynomialFit(n_order)
		self.gap_polynoms_arr = []	
		self.gap_polynoms_coef_arr = []
		self.gap_polynoms_t_tp_s_sp_err_arr = []
		for i_gap in range(len(self.ttp_ssp_gap_arr)):
			[func_T,func_TP,func_S,func_SP] = self.ttp_ssp_gap_arr[i_gap]
			polynomialFit.fitFunction(func_T)
			t_coef_err_arr = polynomialFit.getCoefficientsAndErr()
			t_polynom = polynomialFit.getPolynomial()
			#-------------------------------------------
			polynomialFit.fitFunction(func_TP)
			tp_coef_err_arr = polynomialFit.getCoefficientsAndErr()
			tp_polynom = polynomialFit.getPolynomial()
			#-------------------------------------------
			polynomialFit.fitFunction(func_S)
			s_coef_err_arr = polynomialFit.getCoefficientsAndErr()
			s_polynom = polynomialFit.getPolynomial()
			#-------------------------------------------
			polynomialFit.fitFunction(func_SP)
			sp_coef_err_arr = polynomialFit.getCoefficientsAndErr()
			sp_polynom = polynomialFit.getPolynomial()
			#-------------------------------------------
			self.gap_polynoms_coef_arr.append([t_coef_err_arr,tp_coef_err_arr,s_coef_err_arr,sp_coef_err_arr])
			self.gap_polynoms_arr.append([t_polynom,tp_polynom,s_polynom,sp_polynom])
			#============ calculate max relative errors ===============
			err_t = 0.
			err_tp = 0.
			err_s = 0.
			err_sp = 0.
			n_points = func_T.getSize()
			for i in range(n_points):
				x = func_T.x(i)
				y_t = func_T.y(i)
				y_poly_t = t_polynom.value(x)
				y_tp = func_TP.y(i)
				y_poly_tp = tp_polynom.value(x)
				y_s = func_S.y(i)
				y_poly_s = s_polynom.value(x)
				y_sp = func_SP.y(i)
				y_poly_sp = sp_polynom.value(x)
				if(math.fabs((y_t-y_poly_t)) > err_t): err_t = math.fabs((y_t-y_poly_t))
				if(math.fabs((y_tp-y_poly_tp)) > err_tp): err_tp = math.fabs((y_tp-y_poly_tp))
				if(math.fabs((y_s-y_poly_s)) > err_s): err_s = math.fabs((y_s-y_poly_s))
				if(math.fabs((y_sp-y_poly_sp)) > err_sp): err_sp = math.fabs((y_sp-y_poly_sp))
			self.gap_polynoms_t_tp_s_sp_err_arr.append([err_t,err_tp,err_s,err_sp])
		return self.gap_polynoms_coef_arr

	def dumpTTFandFitting(self,file_ttf_and_fitting_out):
		""" 
		This method prints out the TTF T,Tp,S,Sp and polynomials fitting for all RF gaps.
		"""
		if(len(self.ttp_ssp_gap_arr) == 0 or len(self.gap_polynoms_coef_arr) == 0):
			print "RF_AxisFieldAnalysis::dumpTTFandFitting method:"
			print "Please, call makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq) first!"
			print "Please, call makePlynomialFittings(n_order) first!"
			print "Stop."
			sys.exit(1)
			return
		fl_out = open(file_ttf_and_fitting_out,"w")
		st = "beta_index beta  cappa "
		for i in range(len(self.gap_slpline_arr)):
			st = st + " " + str(i+1) + "T  Tfit  " + str(i+1) + "TP  TPfit  " + str(i+1) + "S  Sfit  " + str(i+1) + "SP    SPfit  "
		fl_out.write(st + "\n")
		for i_beta in range(len(self.beta_arr)):
			beta = self.beta_arr[i_beta]
			cappa = self.cappa_arr[i_beta]
			st = str(i_beta) + " %9.6f "%beta + " %14.7g "%cappa
			for i_gap in range(len(self.gap_slpline_arr)):
				[func_T,func_TP,func_S,func_SP] = self.ttp_ssp_gap_arr[i_gap]
				[t_polynom,tp_polynom,s_polynom,sp_polynom] = self.gap_polynoms_arr[i_gap]
				T = func_T.y(i_beta)
				S = func_S.y(i_beta)
				TP = func_TP.y(i_beta)
				SP = func_SP.y(i_beta)
				Tfit = t_polynom.value(cappa)
				TPfit = tp_polynom.value(cappa)
				Sfit = s_polynom.value(cappa)
				SPfit = sp_polynom.value(cappa)
				st = st + "  %14.7g %14.7g   %14.7g %14.7g   %14.7g %14.7g   %14.7g %14.7g   "%(T,Tfit,TP,TPfit,S,Sfit,SP,SPfit)
			fl_out.write(st + "\n")			
		fl_out.write(st + "\n")		
		fl_out.close()

	def dumpTransitTimeTables(self,file_ttf_out):
		""" 
		This method prints out the final parameters of the transit time tables for all RF gaps.
		"""
		if(len(self.ttp_ssp_gap_arr) == 0 or len(self.gap_polynoms_coef_arr) == 0):
			print "RF_AxisFieldAnalysis::dumpTransitTimeTables(self,file_ttf_out) method:"
			print "Please, call makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq) first!"
			print "Please, call makePlynomialFittings(n_order) first!"
			print "Stop."
			sys.exit(1)
			return
		fl_out = open(file_ttf_out,"w")
		st = "rf_frequency "+str(self.rf_freq)
		fl_out.write(st + "\n")
		st = "beta_min_max %9.6f %9.6f "%(self.beta_arr[len(self.beta_arr)-1],self.beta_arr[0])
		fl_out.write(st + "\n")
		st = "ngaps "+str(len(self.gap_slpline_arr))
		fl_out.write(st + "\n")
		st = "gap_border_points "
		for i in range(len(self.roots_arr)):
			st = st + " %9.6f "%(self.roots_arr[i]-self.rf_center)
		fl_out.write(st + "\n")
		st = "gap_positions "
		for i in range(len(self.gap_slpline_arr)):
			st = st + " %9.6f "%self.gap_slpline_arr[i][1]
		fl_out.write(st + "\n")
		st = "gap_lengths "
		for i in range(len(self.gap_slpline_arr)):
			st = st + " %9.6f "%self.gap_slpline_arr[i][0]
		fl_out.write(st + "\n")
		st = "gap_E0_amplitudes "
		for i in range(len(self.e0_normalized_arr)):
			st = st + " %8.6f "%(self.e0_normalized_arr[i])
		fl_out.write(st + "\n")
		st = "gap_E0L_amplitudes "
		for i in range(len(self.e0l_normalized_arr)):
			st = st + " %8.6f "%(self.e0l_normalized_arr[i])
		fl_out.write(st + "\n")
		for i in range(len(self.gap_polynoms_coef_arr)):
			[t_coef_err_arr,tp_coef_err_arr,s_coef_err_arr,sp_coef_err_arr] = self.gap_polynoms_coef_arr[i]
			st = "gap "+str(i+1)+" polynom Tcoef= "
			[coef_arr,err_arr] = t_coef_err_arr
			for j in range(len(coef_arr)):
				st = st + " %13.6e +- %9.2e "%(coef_arr[j],err_arr[j])
			fl_out.write(st + "\n")
			st = "gap "+str(i+1)+" polynom TPcoef= "
			[coef_arr,err_arr] = tp_coef_err_arr
			for j in range(len(coef_arr)):
				st = st + " %13.6e +- %9.2e "%(coef_arr[j],err_arr[j])
			fl_out.write(st + "\n")
			st = "gap "+str(i+1)+" polynom Scoef= "
			[coef_arr,err_arr] = s_coef_err_arr
			for j in range(len(coef_arr)):
				st = st + " %13.6e +- %9.2e "%(coef_arr[j],err_arr[j])
			fl_out.write(st + "\n")
			st = "gap "+str(i+1)+" polynom SPcoef= "
			[coef_arr,err_arr] = sp_coef_err_arr
			for j in range(len(coef_arr)):
				st = st + " %13.6e +- %9.2e "%(coef_arr[j],err_arr[j])
			fl_out.write(st + "\n")
		"""
		#--- dump of transit time factors -----------------
		st = "beta_index beta  cappa "
		for i in range(len(self.gap_slpline_arr)):
			st = st + " " + str(i+1) + "T(cappa)  " + str(i+1) + "TP(cappa)  " + str(i+1) + "S(cappa)  " + str(i+1) + "SP(cappa)  "
		fl_out.write(st + "\n")
		for i_beta in range(len(self.beta_arr)):
			beta = self.beta_arr[i_beta]
			cappa = self.cappa_arr[i_beta]
			st = str(i_beta) + " %9.6f "%beta + " %14.7g "%cappa
			for i_gap in range(len(self.gap_slpline_arr)):
				[func_T,func_TP,func_S,func_SP] = self.ttp_ssp_gap_arr[i_gap]
				T = func_T.y(i_beta)
				S = func_S.y(i_beta)
				TP = func_TP.y(i_beta)
				SP = func_SP.y(i_beta)
				st = st + "  %14.7g %14.7g  %14.7g %14.7g  "%(T,TP,S,SP)
			fl_out.write(st + "\n")
		"""
		fl_out.close()

