""" heavy-lifting workflow functions"""

# General I/O
import os
import subprocess
import time
import image_base

# fishdot
import generators

# segmentation
from tqdm import tqdm
import tifffile as tf
import numpy as np
import image_operators
import image_augmentations
import image_segmentation_model
import torch


def print_current_time(message=None):
    """
        Prints current time.
    :param message: optional message to show
    :return: none
    """
    current_time = time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(time.time()))
    if message is None:
        print(current_time)
    else:
        print(current_time + " " + message)


def print_check_run_info(run):
    """
        Print and checks run result.
    :param run: subprocess.run object
    """
    print("STDOUT------")
    print(run.stdout)
    print("STDERR------")
    print(run.stderr)
    run.check_returncode()


def fishdot(inputs, output_dir, wildcards, config):
    """
        Given input single channel (fish signal) image full paths, run dot analysis and write to output_dir
        snakemake config params:
            data_roots/generated_scripts, tools/{airlocalize, airlocalize_template_ini, airlocalize_template_m},
            fishdot, threads_per_sample/fishdot, mem_mb_per_sample/fishdot
    :param inputs: FULL paths to single channel images (fish signal)
    :param output_dir: Output directory
    :param wildcards: wildcards involved; this is used to generate script folder path
    :param config: snakemake config
    :return: none
    """
    #  Note - AIRLOCALIZE by default output file names to be the same as inputs, but extension swapped to par4 and loc4.
    os.mkdir(output_dir)
    #  First, generate scripts
    params = "+".join([wildcards.probe, wildcards.sample, wildcards.params])
    params = wildcards.light_train + "++" + params
    script_dir = os.path.join(config["data_roots"]["generated_scripts"], "AIRLOCALIZE", params)
    output_mpath = os.path.join(script_dir, "run.m")
    os.makedirs(script_dir, exist_ok=True)  # For generated script, snakemake will NOT be involved, need to make dirs.
    inputs_abspath = [os.path.abspath(inp) for inp in inputs]  # get absolute path
    output_dir_abspath = os.path.abspath(output_dir)
    print_current_time(f"Generating MATLAB scripts for {params}")
    generators.AIRLOCALIZE_gen(config, inputs_abspath, output_mpath, output_dir_abspath)
    #  Next, call MATLAB for execution
    print_current_time(f"Running MATLAB for {params}")
    matlab_run_string = f"run('{output_mpath}'); exit;"
    run = subprocess.run(["matlab", "-nodisplay -nosplash -nodesktop", "-r", matlab_run_string])
    print_check_run_info(run)
    print_current_time("Done fishdot")


def segmentation(inputs, output_dir, config):
    """
        Given input single channel (marker signal) image full paths, run model and write masks to output_dir.
        This function is a wrapper of the U-Net segmentation workflow.
        snakemake config params:
            segmentation, threads_per_sample/segmentation, mem_mb_per_sample/segmentation
    :param inputs: input image full paths (single channel of the marker)
    :param output_dir: output directory
    :param config: snakemake config
    :return: none
    """
    #   Get segmentation parameters
    try:
        model_name = config["segmentation"]["model"]
        settings = config["segmentation_models"][model_name]

        model_path = os.path.join(config["segmentation"]["model_root"], model_name + ".pt")
        device = config["segmentation"]["device"]
        encoder = settings["encoder"]
        weights = settings["weights"]
        target_pixels = settings["target_pixels"]
        MIP_method = settings["MIP"]
        projection_physical = settings["z_projection_physical"]
        augmentation_function = settings["augmentation_function"]
        augmentation_params = settings["augmentation_function_params"]
        fill = settings["fill"]
    except:
        raise ValueError("Unsupported segmentation model. Check internal configs.")
    #   Initialize segmentation model
    print_current_time("Loading pytorch model")
    model = image_segmentation_model.SegmentationModel(encoder, weights)
    model.to(device)
    model.load_state_dict(torch.load(model_path, map_location=torch.device(device)))
    print_current_time("Model loaded, building output")
    original_names = [os.path.split(inp)[1] for inp in inputs]
    output_paths = [os.path.join(output_dir, oname) for oname in original_names]
    os.mkdir(output_dir)
    print_current_time("Performing prediction")
    for i in tqdm(range(len(output_paths))):
        print_current_time("...input = " + inputs[i])
        # Load tifffile
        image = tf.TiffFile(inputs[i])
        # Get physicalSizeZ information
        try:
            physicalSizeZ = float(image.shaped_metadata[0]["PhysicalSizeZ"])
        except:
            raise ValueError("Input file does not hold valid physical size information.")
        # Perform shape transpose and get original sizes
        image = image.asarray()  # shape = zyx
        image = np.transpose(image, axes=(2, 1, 0))  # new shape = xyz
        image = np.expand_dims(image, axis=2)  # expand channel shape = xycz where c=1
        print_current_time("......input shape after channel expansion=" + str(image.shape))
        original_x_pixels = image.shape[0]  # STOPPOINT----------
        original_y_pixels = image.shape[1]  # original x & y will be useful for resizing back.
        # Perform augmentation
        try:
            augmentations_f = getattr(image_augmentations, augmentation_function)
        except:
            raise ValueError("Augmentation function undefined in image_augmentations.")
        if len(augmentation_params) == 0:
            image = augmentations_f(image, target_pixels=target_pixels)
        else:
            try:
                image = augmentations_f(image, target_pixels=target_pixels, **augmentation_params)
            except:
                raise ValueError("use dict for augmentation function parameters.")
        image = np.squeeze(image, axis=2)
        print_current_time("augmentation completed & channel collapsed, shape=" + str(image.shape))
        # Perform MIP of the whole stack
        image = image_operators.Z_MIP_Physical_Stack(img=image, z_spacing_physical=projection_physical,
                                                     method=MIP_method, physicalSizeZ=physicalSizeZ)
        print_current_time("MIP completed")
        # Run model
        pred_mask_stack = image_segmentation_model.getMask(model=model, image=image, device=device,
                                                           original_xy_pixels={"x": original_x_pixels,
                                                                               "y": original_y_pixels})
        print_current_time("Prediction completed")
        # Optionally fill the image
        if fill:
            for j in range(pred_mask_stack.shape[2]):
                pred_mask_stack[:, :, j] = image_operators.imfill(pred_mask_stack[:, :, j])
            print_current_time("Fill completed")
        tf.imwrite(output_paths[i], np.transpose(pred_mask_stack, axes=(2, 1, 0)), photometric='minisblack')
        print_current_time("Output mask written")
