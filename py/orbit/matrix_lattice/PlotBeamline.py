import string
import sys
import os
from numpy import *


import matplotlib.patches as patches

class PlotBeamline:
	
	def __init__(self, lattice):
		self.lattice = lattice
		
	def get_bline_length(self):
		return self.lattice.getLength()
	

	# plot beamline
	def add_quad(self,ax,s,w,k):
		if k > 0 or k == 0:
			ax.add_patch(
				patches.Rectangle(
					(s, 0.5+0.05),	# (x,y)
					w,					# width
					0.05,				# height
					facecolor="red",	# color
					edgecolor="none"	# No border
					))
		elif k < 0:
			ax.add_patch(
				patches.Rectangle(
					(s, 0.5),			# (x,y)
					w,					# width
					0.05,				# height
					facecolor="red",	# color
					edgecolor="none"	# No border
					))

	
	def add_bend(self,ax,s,w):
		ax.add_patch(
			patches.Rectangle(
				(s, 0.505),				# (x,y)
				w,						# width
				0.08,					# height
				facecolor="blue",		# color
				edgecolor="none"		# No border
			)
		)
		
	def add_kicker(self,ax,s,w):
		ax.add_patch(
			patches.Rectangle(
				(s, 0.505),				# (x,y)
				w,						# width
				0.08,					# height
				facecolor="black",		# color
				edgecolor="none"		# No border
			)
		)
		
	def add_septum(self,ax,s,w):
		ax.add_patch(
			patches.Rectangle(
				(s-w, 0.505),				# (x,y)
				w,						# width
				0.08,					# height
				facecolor="green",		# color
				edgecolor="none"		# No border
			)
		)
		
	
	def plot_beamline(self,ax, l1=0.0, l2=0, title='', septum_name=''):
		for node in self.lattice.getNodes():
			if node.getType() == "quad teapot" or node.getType() == "quad teapot with aperture":
				s = self.lattice.getNodePositionsDict()[node][0]
				length = node.getLength()
				kq = node.getParam("kq")
				self.add_quad(ax,s,length,kq)
			if node.getType() == "bend teapot":
				s = self.lattice.getNodePositionsDict()[node][0]
				length = node.getLength()
				self.add_bend(ax,s,length)
			if node.getType() == "kick teapot":
				s = self.lattice.getNodePositionsDict()[node][0]
				length = node.getLength()
				self.add_kicker(ax,s,0.25)
			if node.getName() == septum_name:
				s = self.lattice.getNodePositionsDict()[node][0]
				self.add_septum(ax,s,1.0)
		if l2==0.0:
			l2 = self.get_bline_length()
		ax.text(0, 1,title, horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
		ax.set_xlim(l1,l2)
		ax.set_ylim(0.45,0.65)
		ax.axhline(y=0.5+0.05, color='black')
		ax.axis('off')

