#  Creates a NAMD control file, from a template file, which cools and equilibrates the system using a find/replace.
import os

def get_sphere_builder_path():
    return str(os.path.dirname(os.path.abspath(__file__))) + "/makeSphere_template.tcl"
def get_water_box_builder_path():
    return str(os.path.dirname(os.path.abspath(__file__))) + "/makeWaterBox_template.tcl"
def get_pdb2bincoords_path():
    return str(os.path.dirname(os.path.abspath(__file__))) + "/xyz2bincoords_template.vmd"
def get_pdb2xsc_path():
    return str(os.path.dirname(os.path.abspath(__file__))) + "/pdb2xsc_template.vmd"
