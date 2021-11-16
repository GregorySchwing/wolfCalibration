#!/usr/bin/python

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.optimize import minimize, brute
from scipy import interpolate, optimize
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle, Ellipse

def add_point(ax, x, y, z, fc = None, ec = None, radius = 0.005, labelArg = None):
	xy_len, z_len = ax.get_figure().get_size_inches()
	axis_length = [x[1] - x[0] for x in [ax.get_xbound(), ax.get_ybound(), ax.get_zbound()]]
	axis_rotation =  {'z': ((x, y, z), axis_length[1]/axis_length[0]),
			'y': ((x, z, y), axis_length[2]/axis_length[0]*xy_len/z_len),
			'x': ((y, z, x), axis_length[2]/axis_length[1]*xy_len/z_len)}
	i = 0
	for a, ((x0, y0, z0), ratio) in axis_rotation.items():
		p = Ellipse((x0, y0), width = radius, height = radius*ratio, fc=fc, ec=ec, label = labelArg if i == 0 else "")
		ax.add_patch(p)
		art3d.pathpatch_2d_to_3d(p, z=z0, zdir=a)
		i = i + 1

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print('findOptimumRCutWolfAlpha.py -i <intputfile.p> -o <outputfile>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('findOptimumRCutWolfAlpha.py -i <intputfile.p> -o <outputfile>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
			parsingInputs = False
	print('Input file is', inputfile)
	print('Output file is ', outputfile)

	df = pd.read_pickle(inputfile)
	method_pots = []
	methods = [ "Gross" , "Vlugt", "Cassandra" ]
	potentials = [ "DSP" , "DSF"]
	tuples = list()
	for method in methods:
		for potential in potentials:
			colName = "{}_{}".format(method, potential)
			method_pots.append(colName)
	for metpot in method_pots:
		df2 = df.iloc[:,df.columns.get_level_values(1)==(metpot)]
		alphas = df.iloc[:,df.columns.get_level_values(1)==('Alpha')]
		#print(alphas)

		df2 = pd.MultiIndex.from_frame(df2)
		Z = np.vstack(df2)

		x = []
		y = []
		z = []
		pointsbottom = []
		for name in df2.names:
			for alpha, value in zip(alphas[(name[0], "Alpha")], df2.get_level_values(name)):
				x.append(alpha) 
				y.append(name[0]) 
				z.append(value)

		
		x2 = np.linspace(0.00, 0.30, 65)
		y2 = np.linspace(10, 15, 65)
		rranges = slice(0.00, 0.30, 0.3/65), slice(10, 15, 5/65)
		X2, Y2 = np.meshgrid(x2, y2)


		F2 = interpolate.interp2d(x, y, z, kind='quintic')
		Z2 = F2(x2, y2)
		f = lambda x: np.abs(F2(*x))
		
		x0 = (0.1, 12)
		bounds = [(0.0,0.3),(10,15)]
		gd = minimize(f, x0, method='SLSQP', bounds=bounds)
		bf = brute(f, rranges, full_output=True, finish=optimize.fmin)

		print(gd)
		gdXY = np.array(gd.x)
		print(gdXY[0])
		print(gdXY[1])
		bfXY = np.array(bf[0])
		print(bfXY[0])
		print(bfXY[1])
		ax = plt.axes(projection='3d')
		ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none');
		ax.set_title(metpot, fontsize=20)
		ax.set_xlabel('Alpha', fontsize=20, labelpad=20)
		ax.set_ylabel('RCut', fontsize=20, labelpad=20)
		ax.set_zlabel('Relative Error', fontsize=20, labelpad=20)
		xbf,ybf = bf[0]
		add_point(ax, gdXY[0], gdXY[1], gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01, labelArg = "Gradient Descent")
		add_point(ax, gdXY[0], gdXY[1], 10*gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01)
		add_point(ax, gdXY[0], gdXY[1], 20*gd.fun[0], fc = 'orange', ec = 'orange', radius=0.01)

		add_point(ax, bfXY[0], bfXY[1], bf[1], fc = 'r', ec = 'r', radius=0.01, labelArg = "Brute Force")
		add_point(ax, bfXY[0], bfXY[1], 10*bf[1], fc = 'r', ec = 'r', radius=0.01)
		add_point(ax, bfXY[0], bfXY[1], 20*bf[1], fc = 'r', ec = 'r', radius=0.01)
		ax.legend(loc='best')
		plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
