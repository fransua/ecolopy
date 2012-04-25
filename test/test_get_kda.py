import sys
from ecolopy_dev import Abundance

abd = Abundance ('../dataset_trial/bci_full.txt')

bip = sys.argv[1]=='bip'

if bip:
    print 'bip'
    abd._get_kda_bip()
else:
    abd._get_kda()

abd.etienne_optimal_params()
abd.set_current_model('etienne')
print abd

#print abd._kda[:10]
#print abd._kda[20872:]
#print '\n'
#print len (abd._kda)
