from masked_cubes import cube303m,cube321m,cube303,cube321,mask
from noise import noise, noise_cube

#noise = fits.getdata(hpath('APEX_H2CO_303_202_noise.fits'))
#noise = cube303[:50].std(axis=0).value
#noise = fits.getdata(mpath('APEX_H2CO_merge_high_sub_noise.fits'))
#nhits = nhits = fits.getdata(paths.mpath('APEX_H2CO_merge_high_nhits.fits'))
#noise[nhits<20] = np.nan

noise_flat = noise_cube[mask]
var_flat = noise_flat**2

ratio303321 = cube321m.flattened().value / cube303m.flattened().value
eratio303321 = (ratio303321**2 * (var_flat/cube303m.flattened().value**2 + var_flat/cube321m.flattened().value**2))**0.5
