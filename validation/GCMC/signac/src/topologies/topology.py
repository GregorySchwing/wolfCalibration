#  Creates a NAMD control file, from a template file, which cools and equilibrates the system using a find/replace.
import os

def get_topology_path(topologyFile):
    return str(os.path.dirname(os.path.abspath(__file__))) + "/" + topologyFile
