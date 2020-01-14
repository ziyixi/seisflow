"""
plot_events.py: plot events beachballs, the color represents the value to show.
"""

from jinja2 import Environment, PackageLoader
from os.path import join, basename
import tempfile
from ...utils.load_files import load_pickle

env = Environment(
    loader=PackageLoader('seisflow', 'plot/gmt/template')
)


def prepare_environment(cpt_type, value_pickle_path, cpt_color, cpt_range, flag_J, flag_R,
                        colorbar_name, colorbar_unit, output_path):
    """
    convert the template to the bash file, prepare all the essencial files to run the gmt.
    """
    value_used = load_pickle(value_pickle_path)
    # * firstly we prepare all the required parameters in the template
    cpt_1, cpt_2, cpt_3 = False, False, False
    if (cpt_type == 1):
        cpt_1 = True
    elif(cpt_type == 2):
        cpt_2 = True
    elif(cpt_type == 3):
        cpt_3 = True
    else:
        raise Exception("no such cpt_type")
    tmpdir = tempfile.mkdtemp()
    cpt_min, cpt_max = cpt_range.split(",")
    # determine cpt_step value

    # * required parameters for the template
    # tmpdir, cpt_1, cpt_2, cpt_3, cpt_color, cpt_min, cpt_max, cpt_step,
    # output_path, flag_J, flag_R, value_min, value_max, colorbar_name, colorbar_unit
    template = env.get_template('plot_events.j2')
    template.render({
        "tmpdir": tmpdir,
        "cpt_1": cpt_1,
        "cpt_2": cpt_2,
        "cpt_3": cpt_3,
        "cpt_color": cpt_color,
        "cpt_min": cpt_min,
        "cpt_max": cpt_max,
        "cpt_step": cpt_step,
        "output_path": output_path,
        "flag_J": flag_J,
        "flag_R": flag_R,
        "value_min": value_min,
        "value_max": value_max,
        "colorbar_name": colorbar_name,
        "colorbar_unit": colorbar_unit
    })
