Text_File_Import = open('GCMC_water_O2_res_seg_fix.psf', 'r')

Text_lines = Text_File_Import.readlines()
with open('GCMC_water_O2_res_seg_fix_atom_type_fix.psf', 'w') as the_file:
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
			if (User_Inputs[resName+2] == "A"):
				newLine = line.replace('A ', "HT")
				the_file.write(newLine)
			elif (User_Inputs[resName+2] == "B"):
				newLine = line.replace('B ', "ON")
				the_file.write(newLine)
			elif (User_Inputs[resName+2] == "C"):
				newLine = line.replace('C ', "OP")
				the_file.write(newLine)
			elif (User_Inputs[resName+2] == "D"):
				newLine = line.replace('D ', "OT")
				the_file.write(newLine)
		else:
			the_file.write(line)




