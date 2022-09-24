#  Creates a NAMD control file, from a template file, which cools and equilibrates the system using a find/replace.
import os

def generate_namd_equilb_control_file(template_path_filename=None,
                                      namd_path_conf_filename=None,
                                      namd_path_file_output_names=None,
                                      namd_uses_water=None,
                                      namd_water_model=None,
                                      namd_electrostatics_bool=None,
                                      namd_vdw_geometric_sigma_bool=None,
                                      namd_psf_path_filename=None,
                                      namd_pdb_path_filename=None,
                                      namd_ff_path_filename=None,
                                      namd_production_temp_K=None,
                                      namd_production_pressure_bar=None,
                                      electrostatic_1_4=None,
                                      non_bonded_cutoff=None,
                                      non_bonded_switch_distance=None,
                                      pairlist_distance=None,
                                      box_lengths=None,
                                      ):
    """
    Creates a NAMD control file, from a template file, which cools and equilibrates the system using a find/replace.

    template_path_filename : str
        The NAMD control file template file path and name
    namd_path_conf_filename : str
        The path and file name for the NAMD control file name, with
        the .conf extension, or no extension.  If no extension is provided, the
        code will add the .conf extension to the provided file name.
    namd_path_file_output_names; str
        The name or partial name of the NAMD created output files
    namd_uses_water: bool
        Does it use a water model, True or False
    namd_water_model: str or None,  ( tip3, tip4, swm4) or None
        select the water model if used or None
    namd_electrostatics_bool: bool
        Whether to use electrostatics in NAMD
    namd_vdw_geometric_sigma_bool: bool
         Whether to use geometric mixing.  True give geometric mixing, and False give LB mixing
    namd_psf_path_filename : str
        The namd psf file name including the path
    namd_pdb_path_filename : str
        The namd pdb file name including the path
    namd_ff_path_filename : str
        The namd force field file name including the path
    namd_production_temp_K : int
        The production simulation temperature in Kelvin
    namd_production_pressure_bar : float or int
        The production simulation pressure in bar
    electrostatic_1_4: float
        The electrostatic 1-4 scalar
    non_bonded_cutoff: float
        The non-bonded cutoff distance in angstrom
    non_bonded_switch_distance: float
        The non-bonded switch distance in angstrom
    pairlist_distance: float
        The pairlist distance in angstrom
    box_lengths : list or 3 floats,  [x_dist, y_dist, z_dist]
        The box lengths of the system, which assumes that it is an orthogonal box.
    """
    variable_dict = {template_path_filename: template_path_filename,
                     'namd_path_conf_filename': namd_path_conf_filename,
                     'namd_path_file_output_names': namd_path_file_output_names,
                     'namd_electrostatics_bool': namd_electrostatics_bool,
                     'namd_vdw_geometric_sigma_bool': namd_vdw_geometric_sigma_bool,
                     'namd_psf_path_filename': namd_psf_path_filename,
                     'namd_pdb_path_filename': namd_pdb_path_filename,
                     'namd_ff_path_filename': namd_ff_path_filename,
                     'namd_production_temp_K': namd_production_temp_K,
                     'namd_production_pressure_bar': namd_production_pressure_bar,
                     'electrostatic_1_4': electrostatic_1_4,
                     'non_bonded_cutoff': non_bonded_cutoff,
                     'non_bonded_switch_distance': non_bonded_switch_distance,
                     'pairlist_distance': pairlist_distance,
                     'box_lengths': box_lengths,
                     }

    for variable_i in list(variable_dict.keys()):
        if variable_dict[variable_i] is None:
            print()
            raise ValueError(f"ERROR: The {variable_i} variable needs to be provided.")

    if os.path.splitext(namd_path_conf_filename)[1] == '.conf':
        print(f'INFO: the correct extension for the control file was provided in the file name, .conf '
              'with control file name = {namd_path_conf_filename}'
              )
    elif os.path.splitext(namd_path_conf_filename)[1] == '':
        namd_path_conf_filename = namd_path_conf_filename + '.conf'
        print(f'INFO: No extension name was provided for the control file. Therefore, the proper '
              'extension, .conf, was added.  The new total control file name = {namd_path_conf_filename}'
              )

    for check_string_i in [template_path_filename,
                           namd_path_conf_filename,
                           namd_path_file_output_names,
                           namd_psf_path_filename,
                           namd_pdb_path_filename,
                           namd_ff_path_filename
                           ]:
        if not isinstance(check_string_i, str):
            raise TypeError("ERROR: The template_path_filename, "
                            "namd_path_conf_filename, "
                            "namd_path_file_output_names, "
                            "namd_psf_path_filename, "
                            "namd_pdb_path_filename, "
                            "or namd_ff_path_filename "
                            "is not a string.")

    if not isinstance(namd_uses_water, bool):
        raise TypeError(f"ERROR: Select True of False for the namd_uses_water variable.")
    elif namd_uses_water is True and namd_water_model not in ['tip3', 'tip4', 'swm4']:
        raise ValueError(f"ERROR: If the namd_uses_water in not None, then the "
                          f"namd_water_model must be provide one of these ['tip3', 'tip4', 'swm4'].")

    if not isinstance(namd_vdw_geometric_sigma_bool, bool):
        raise TypeError(f"ERROR: Select the a bool for the namd_vdw_geometric_sigma_bool variable.")
    elif namd_vdw_geometric_sigma_bool is True:
        namd_vdw_geometric_sigma_yes_no = 'yes'
    elif namd_vdw_geometric_sigma_bool is False:
        namd_vdw_geometric_sigma_yes_no = 'no'

    if not isinstance(namd_electrostatics_bool, bool):
        raise TypeError(f"ERROR: Select the a bool for the namd_electrostatics_bool variable.")

    for check_positive_i in [namd_production_temp_K]:
        if isinstance(check_positive_i, int):
            if check_positive_i < 0:
                raise ValueError(f"ERROR: The  namd_production_temp_K, "
                                  f"is a negative value.")
        else:
            raise ValueError(f"ERROR: The namd_production_temp_K, "
                             f"is not a integer.  "
                             f"NOTE: The namd_production_temp_K are "
                             f"need to be integers so they can be ramped down from hot to cool "
                             f"to equilibrate the system.")

    for check_positive_i in [namd_production_pressure_bar]:
        if isinstance(check_positive_i, float) or isinstance(check_positive_i, int):
            if check_positive_i < 0:
                raise ValueError(f"ERROR: The namd_production_pressure_bar "
                                  f"is a negative value.")
        else:
            raise ValueError(f"ERROR: The namd_production_pressure_bar "
                              f"is not a float or integer.")

    if isinstance(electrostatic_1_4, float) or isinstance(electrostatic_1_4, int):
        if not (electrostatic_1_4 >= 0 and electrostatic_1_4 <= 1):
            raise ValueError(f"ERROR: The electrostatic_1_4 not a float or int from 0 to 1.")
    else:
        raise TypeError(f"ERROR: The electrostatic_1_4 value is not a float or integer.")

    for check_positive_j in [non_bonded_cutoff, non_bonded_switch_distance, pairlist_distance]:
        if isinstance(check_positive_j, float) or isinstance(check_positive_j, int):
            if check_positive_j <= 0:
                raise ValueError(f"ERROR: The non_bonded_cutoff, non_bonded_switch_distance, or pairlist_distance "
                                  f"is <= 0.")
        else:
            raise TypeError(f"ERROR: The non_bonded_cutoff, non_bonded_switch_distance, or pairlist_distance "
                              f"is not a float or integer.")

    if non_bonded_switch_distance >= non_bonded_cutoff or non_bonded_cutoff >= pairlist_distance \
            or non_bonded_switch_distance >= pairlist_distance:
        raise ValueError(f"ERROR: This must be true, "
                          f"non_bonded_switch_distance < non_bonded_cutoff < pairlist_distance")

    if isinstance(box_lengths, list):
        if len(box_lengths) !=3:
            raise ValueError(f"ERROR: The box_lengths is not a list of length 3.")
        for check_positive_k in box_lengths:
            if isinstance(check_positive_k, float) or isinstance(check_positive_k, int):
                if check_positive_k <= 0:
                    raise ValueError(f"ERROR: The box_lengths has a value that is <= 0.")
            else:
                raise TypeError(f"ERROR: The box_lengths has a value that is not a float or integer.")
    else:
        raise TypeError(f"ERROR: The box_lengths is not a list.")


    # build input file from template
    file_in = open(template_path_filename, 'r')
    filedata = file_in.read()
    file_in.close()

    written_to_file = open(namd_path_conf_filename, 'w')
    if namd_uses_water is None or namd_uses_water is False:
        newdata = filedata.replace("NAMD_WATER_MODEL", "")
    else:
        newdata = filedata.replace("NAMD_WATER_MODEL", str(f"waterModel {namd_water_model}"))

    newdata = newdata.replace("NAMD_PSF", str(namd_psf_path_filename))
    newdata = newdata.replace("NAMD_PDB", str(namd_pdb_path_filename))

    newdata = newdata.replace("NAMD_FF", str(namd_ff_path_filename))

    newdata = newdata.replace("NAMD_MIXING_GEOMETRIC", str(namd_vdw_geometric_sigma_yes_no))
    newdata = newdata.replace("NAMD_USE_ELECTROSTATICS", str(namd_electrostatics_bool))

    newdata = newdata.replace("NAMD_OUTPUT_FILE_NAMES", str(namd_path_file_output_names))

    newdata = newdata.replace("PRODUCTION_TEMP", str(namd_production_temp_K))

    newdata = newdata.replace("PRODUCTION_PRESSURE", str(namd_production_pressure_bar))

    newdata = newdata.replace("1_4_SCALING", str(electrostatic_1_4))
    newdata = newdata.replace("CUTOFF_DIST_ANG", str(non_bonded_cutoff))
    newdata = newdata.replace("SWITCH_DIST_ANG", str(non_bonded_switch_distance))
    newdata = newdata.replace("PAIR_LIST_DIST_ANG", str(pairlist_distance))

    newdata = newdata.replace("BOX_DIM_X_ANG", str(box_lengths[0]))
    newdata = newdata.replace("BOX_DIM_Y_ANG", str(box_lengths[1]))
    newdata = newdata.replace("BOX_DIM_Z_ANG", str(box_lengths[2]))

    newdata = newdata.replace("BOX_CENTER_X_ANG", str(box_lengths[0]/2))
    newdata = newdata.replace("BOX_CENTER_Y_ANG", str(box_lengths[1]/2))
    newdata = newdata.replace("BOX_CENTER_Z_ANG", str(box_lengths[2]/2))

    written_to_file.write(newdata)
    written_to_file.close()


"""

generate_namd_equilb_control_file(template_path_filename='NAMD_conf_template.conf',
                                  namd_path_conf_filename='../sNAMD.conf',
                                      namd_uses_water=True,
                                      namd_water_model='tip4',
                                      namd_psf_path_filename='NPT_water.psf',
                                      namd_pdb_path_filename='NPT_water.pdb',
                                      namd_ff_path_filename='NPT_water.inp',
                                      namd_production_temp_K=298,
                                      namd_production_pressur_bare=1,
                                      electrostatic_1_4=0.5,
                                      non_bonded_cutoff=16,
                                      non_bonded_switch_distance=14,
                                      pairlist_distance=20,
                                      box_lengths=[40, 40, 40]
                                  )


"""

