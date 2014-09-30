files = ['dendro_mask',
         'dendro_temperature',
         'make_temperature_cube']

for pre in files:
    execfile(pre+".py")
