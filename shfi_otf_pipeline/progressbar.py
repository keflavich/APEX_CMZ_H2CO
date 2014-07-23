"""
No progress bars in notebooks
"""
try:
    cfg = get_ipython().config
    if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
        ProgressBar = lambda x: x
    else:
        from astropy.utils.console import ProgressBar
except NameError: # Not in a notebook
    from astropy.utils.console import ProgressBar

