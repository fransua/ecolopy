import sys
from ecolopy import Abundance

abd = Abundance ('../dataset_trial/bci_full.txt')

bip = sys.argv[1]=='bip'

if bip:
    print 'bip'
    abd._get_kda_bip()
else:
    abd._get_kda()

