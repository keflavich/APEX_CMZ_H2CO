import pylab as pl
# name, tdust, tgas, etgas
galaxy_data = [('NGC 253', 34, 78, 22),
               ('NGC 660', 37, 160, 96),
               ('NGC 891', 28, 30, 30),
               ('Maffei 2', 40, 43, 6),
               ('Maffei 2', 40, 79, 6),
               ('NGC 1365', 32, 50, 11),
               ('IC 342', 30, 150, 0),
               ('M82 SW', 45, 58, 19),
               #('NGC 3079', 32, 150, 150),
               ('IC 860', 40, 206, 79),
               ('M 83', 31, 56, 15),
               ('IR 15107+0724', 40, 189, 57),
               ('Arp 220', 44, 234, 52),
               ('NGC 6946', 30, 47, 8),
               ('CMZ', 23, 64, 16),
              ]

pl.clf()
for galaxy in galaxy_data:
    name, td, tg, etg = galaxy
    pl.errorbar(td, tg, yerr=etg, marker='s')
    pl.text(td, tg, name, va='center', ha='center', size=14)

pl.plot([0,50],[0,50],'k--', linewidth=3, zorder=-5, alpha=0.5)
# just the CMZ
pl.text(td, tg, " " + name, va='center', ha='left', size=14)
pl.xlabel("$T_{dust}$")
pl.ylabel("$T_{gas}$")
pl.xlim(20,48)

pl.show()
