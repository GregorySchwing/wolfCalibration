Text_File_Import = open('GCMC_water_O2_liq.psf', 'r')

Text_lines = Text_File_Import.readlines()
counter = 00
segID = "W00"
lastResid = 0
with open('GCMC_water_O2_liq_seg_fix.psf', 'w') as the_file:
	for line in Text_lines:
		User_Inputs = line.split()
		resName = 0
		try:
			resName = User_Inputs.index("TIP3")
		except ValueError:
			"print no tip3"		
		try:
			resName = User_Inputs.index("DIOX")
		except ValueError:
			"print no diox"	

		if (resName > 0):
			if (User_Inputs[resName-1] == str(1) and lastResid == str(9999)):
				counter = counter + 1
				segID = "W{0:0=2d}".format(counter)
			newLine = line.replace('SYS', segID)
			the_file.write(newLine)
			lastResid = User_Inputs[resName-1]
		else:
			the_file.write(line)




