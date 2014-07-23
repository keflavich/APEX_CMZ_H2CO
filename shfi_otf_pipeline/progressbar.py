"""
No progress bars in notebooks
"""
try:
    cfg = get_ipython().config
    if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
        # monkeypatch update so it won't do anything
        print "Running from a notebook.  Progress bar deactivated."
        from astropy.utils import console
        def update(self, *args, **kwargs):
            pass
        console.ProgressBar.update = update
        ProgressBar = lambda x: x
    else:
        from astropy.utils.console import ProgressBar
except NameError: # Not in a notebook
    from astropy.utils.console import ProgressBar
