import sys
import pandas as pd
import numpy as np
argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}
nLoci=int(argv['nLoci'])
mu=float(argv['mu'])
Nref=int(argv['Nref'])
L=int(argv['locus_length'])
nameA=argv['nameA']
nameB=argv['nameB']
size_popA=int(argv['size_popA'])
size_popB=int(argv['size_popB'])
rho_over_theta=float(argv['rho_over_theta'])
timeStamp=argv['timeStamp']

locus_data = pd.DataFrame({'locus_length':np.full(nLoci,L),'size_popA':np.full(nLoci,size_popA),
        'size_popB':np.full(nLoci,size_popB),'totpopsize':np.full(nLoci,size_popA + size_popB),
        'theta':np.full(nLoci, 4 * Nref * L * mu),'rho':np.full(nLoci,4 * Nref * L * mu * rho_over_theta)})

locus_data.to_csv('{0}/locus_datafile'.format(timeStamp),sep='\t',index=False)
with open('{0}/Nref.txt'.format(timeStamp),"w") as o:
    o.write("{0}\n".format(Nref))

with open('{0}/nLoci.txt'.format(timeStamp),"w") as o:
    o.write("{0}\n".format(nLoci))
