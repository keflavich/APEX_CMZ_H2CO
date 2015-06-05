from astropy import units as u
from astropy.io import ascii

def exp_to_tex(st):
    if st == 'nan':
        return '-'
    elif 'e' in st:
        pt1,pt2 = st.split('e')
        return "{0}\\ee{{{1:d}}}".format(pt1,int(pt2))
    return st

def format_float(st):
    return exp_to_tex("{0:0.2g}".format(st))


latexdict = ascii.latex.latexdicts['AA']
latexdict['tabletype'] = 'table*'
latexdict['tablealign'] = 'htp'

