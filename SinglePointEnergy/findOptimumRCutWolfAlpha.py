#!/usr/bin/python

import pandas as pd
import sys, getopt
import matplotlib.pyplot as plt
import numpy as np

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
	
	df2 = df.iloc[:,df.columns.get_level_values(1)==('Vlugt_DSP')]
	alphas = df.iloc[:,df.columns.get_level_values(1)==('Alpha')]
	print(alphas)

	df2 = pd.MultiIndex.from_frame(df2)
	print(df2)	
	print(df2.names)
	x = []
	y = []
	z = []
	for name in df2.names:
		for alpha, value in zip(alphas[(name[0], "Alpha")], df2.get_level_values(name)):
			x.append(alpha) 
			y.append(name[0]) 
			z.append(value)


	ax = plt.axes(projection='3d')
	ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none');
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
