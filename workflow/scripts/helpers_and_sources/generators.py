"""
    Script generators for external calls
"""

import os

#  Constant template file paths and the AIRLOCALIZE external script path
AIRLOCALIZE_template_ini_path = "workflow/scripts/helpers_and_sources/airlocalize_template.ini"
AIRLOCALIZE_template_m_path = "workflow/scripts/helpers_and_sources/airlocalize_template.m"
AIRLOCALIZE_root = "resources/AIRLOCALIZE-main"


def AIRLOCALIZE_gen(config, input_tifs, output_mpath, output_dir_path):
    """
        Generates MATLAB code for AIRLOCALIZE run of a list of input tifs.

    :param output_dir_path: where should AIRLOCALIZE output to
    :param config: smk config
    :param input_tifs: list of tif
    :param output_mpath: string to a .m file, configs of AIRLOCALIZE will be saved to the same folder
    :return: none
    """
    output_base = os.path.split(output_mpath)[0]
    ini_paths = []
    for i in range(len(input_tifs)):
        # Generates ini file for each input tif
        tif_name = os.path.split(input_tifs[i])[1]
        ini_paths.append(os.path.join(output_base, tif_name + ".ini"))
        f = open(AIRLOCALIZE_template_ini_path, "r")  # read template
        t = f.read()
        f.close()
        # Replace the predefined values
        t = t.replace("{INPUT_SINGLE_FILE}", input_tifs[i])
        t = t.replace("{OUTPUT_DIRECTORY}", output_dir_path)
        t = t.replace("{PSF_XY}", config["fishdot"]["psf_xy"])
        t = t.replace("{PSF_Z}", config["fishdot"]["psf_z"])
        t = t.replace("{THRESHOLD}", config["fishdot"]["threshold"])
        f = open(ini_paths[i], "w")  # write template
        f.write(t)
        f.close()
    # Now, generate the master .m for MATLAB execution
    f = open(AIRLOCALIZE_template_m_path, "r")  # read template
    t = f.read()
    f.close()
    #   Note - when saving the ini paths to the .m script, get absolute paths
    ini_paths_abspath = [os.path.abspath(p) for p in ini_paths]
    #   When doing str(list), strings will be quoted single.
    #       Making sure the original file paths do NOT have any single quotes
    quote_check = [p.find("'") != -1 for p in ini_paths_abspath]
    if any(quote_check):
        print(ini_paths_abspath)
        raise ValueError("Please remove all single quotes in your absolute paths.")
    #       Convert single quotes to double
    t = t.replace("{input_inis}", str(ini_paths_abspath).replace("'", '"'))
    #   Note - also get absolute path for the airlocalize root directory
    airlocalize_root = AIRLOCALIZE_root
    airlocalize_root_abspath = os.path.join(os.getcwd(), airlocalize_root)
    t = t.replace("{airlocalize_root}", airlocalize_root_abspath)
    f = open(output_mpath, "w")  # write master .m file
    f.write(t)
    f.close()
