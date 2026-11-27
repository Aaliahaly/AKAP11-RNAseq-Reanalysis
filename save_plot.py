"""
===================================================================
save_plot.py
===================================================================

WHAT THIS FILE DOES
    It contains ONE small function, save_plot(), that saves a figure
    to a .png file and then closes it.

WHERE IT FITS
    Every plotting script ends with a line like:

        save_plot(fig, output / "03_pca.png")

    Putting the saving code in one place means every plot in the
    project is saved the same way (same sharpness, same margins).
    To change how ALL plots are saved, you only change this file.

PYTHON IDEAS USED IN THIS FILE
    import ... as ...   load a library and give it a short nickname
    def                 start a new function (a named, reusable
                        block of code)
    arguments           the values you hand to a function when you
                        call it (here: the figure and the file path)
===================================================================
"""

# matplotlib is Python's main plotting library.
# "pyplot" is the part of it we use to make and close figures.
# "as plt" lets us type the short name plt instead of
# matplotlib.pyplot every time.
import matplotlib.pyplot as plt

# The image sharpness (dots per inch) comes from settings.py, so every
# figure in the project is saved at the same resolution.
from settings import image_resolution    # 200


def save_plot(fig, path):
    """Save a figure to a file, then close it."""

    # This function receives two things (its "arguments"):
    #
    #   fig  -> the figure to save. Each plotting script makes one with
    #           a line like:  fig, ax = plt.subplots(...)
    #
    #   path -> where to save it: a folder plus a file name, for example
    #           akap11_go_bp_results/03_pca.png
    #           The file ending (.png) tells matplotlib which image
    #           format to write.

    # dpi means "dots per inch": how many pixels for each inch of the
    # figure. image_resolution is set in settings.py (200 gives a sharp
    # image that still has a modest file size).

    # Write the figure to the file.
    #
    # bbox_inches="tight" trims the extra white border around the
    # figure. It also makes sure labels or legends that stick out past
    # the plot's edge are kept, instead of being cut off.
    fig.savefig(path, dpi=image_resolution, bbox_inches="tight")

    # Close the figure now that it is saved.
    #
    # Every open figure uses computer memory. Our analysis makes more
    # than 15 figures, so without closing them they would all stay in
    # memory until the program finished (matplotlib also prints a
    # warning when more than 20 figures are open).
    plt.close(fig)
