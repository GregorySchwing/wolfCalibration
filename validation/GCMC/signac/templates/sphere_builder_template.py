#  Creates a NAMD control file, from a template file, which cools and equilibrates the system using a find/replace.
import os

def get_sphere_builder_path():
    return str(os.path.dirname(os.path.abspath(__file__))) + "/makeSphere_template.tcl"
