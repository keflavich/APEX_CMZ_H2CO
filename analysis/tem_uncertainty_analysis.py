from dendrograms import (catalog, catalog_sm, dend, dendsm)
from common_constants import logabundance,elogabundance

from constrain_parameters import paraH2COmodel
if 'mf' not in locals():
    mf = paraH2COmodel()

def sign(x):
    return -1 if x < 0 else 1

for row in catalog:
    if not row['is_leaf']: continue

    mf.set_constraints(ratio303321=row['r303321'], eratio303321=row['er303321'],
                       logh2column=row['logh2column'], elogh2column=row['elogh2column'],
                       logabundance=logabundance, elogabundance=elogabundance,
                       taline303=row['Smean303'], etaline303=row['e303'],
                       taline321=row['Smean321'], etaline321=row['e321'],
                       linmindens=row['dustmindens'],
                       linewidth=5)
    base_row_data = mf.get_parconstraints()

    # lower linewidth = higher abundance
    mf.set_constraints(ratio303321=row['r303321'], eratio303321=row['er303321'],
                       logh2column=row['logh2column'], elogh2column=row['elogh2column'],
                       logabundance=logabundance+np.log10(4.), elogabundance=elogabundance,
                       taline303=row['Smean303'], etaline303=row['e303'],
                       taline321=row['Smean321'], etaline321=row['e321'],
                       linmindens=row['dustmindens'],
                       linewidth=5)
    low_row_data = mf.get_parconstraints()

    # higher linewidth = lower abundance
    mf.set_constraints(ratio303321=row['r303321'], eratio303321=row['er303321'],
                       logh2column=row['logh2column'], elogh2column=row['elogh2column'],
                       logabundance=logabundance-np.log10(4.), elogabundance=elogabundance,
                       taline303=row['Smean303'], etaline303=row['e303'],
                       taline321=row['Smean321'], etaline321=row['e321'],
                       linmindens=row['dustmindens'],
                       linewidth=5)
    high_row_data = mf.get_parconstraints()

    print "Leaf: {5:7} S/N={0:8.3f}: LOW: {1:8.3f}%  {2:8.3f}K HIGH: {3:8.3f}%  {4:8.3f}K  BASE={6:8.3f}K".format(row['Smean303']/row['e303'],
                                      (low_row_data['temperature_chi2']-base_row_data['temperature_chi2'])/base_row_data['temperature_chi2']*100,
                                      low_row_data['temperature_chi2'],
                                      (high_row_data['temperature_chi2']-base_row_data['temperature_chi2'])/base_row_data['temperature_chi2']*100,
                                      high_row_data['temperature_chi2'],
                                                                                                 str(row['is_leaf']),
                                                                                                                base_row_data['temperature_chi2'])
    if sign((low_row_data['temperature_chi2']-base_row_data['temperature_chi2'])) == sign((high_row_data['temperature_chi2']-base_row_data['temperature_chi2'])):
        if ((low_row_data['temperature_chi2']-base_row_data['temperature_chi2']) == 0
            or 
            (high_row_data['temperature_chi2']-base_row_data['temperature_chi2']) == 0):
            continue
        pass
        #print "base: ",base_row_data
        #print "low: ",low_row_data
        #print "high: ",high_row_data
        #break

