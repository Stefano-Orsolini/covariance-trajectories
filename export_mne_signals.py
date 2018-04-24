import numpy as np
import scipy.io as spio
from mne.datasets import eegbci
from mne.io.edf import read_raw_edf
from mne.io import concatenate_raws
from mne.event import find_events
#reference: http://www.martinos.org/mne/stable/generated/mne.datasets.eegbci.load_data.html

# subject range: 1-109
subject=1

# Motor imagery: left vs right hand
runs=[4, 8, 12]
## Motor imagery: hands vs feet
#runs=[6, 10, 14]

# sampling frequency
fs = 160

# eeg signals
raw_fnames = eegbci.load_data(subject, runs)
raws = [read_raw_edf(f, preload=True) for f in raw_fnames]
raw = concatenate_raws(raws)
# experimental paradigm sequence
S_events = find_events(raw, shortest_event=0, stim_channel='STI 014')

# electrodes signals
S_1=raw.get_data(picks=[raw.ch_names.index('C3..')])*1e6
S_2=raw.get_data(picks=[raw.ch_names.index('Cz..')])*1e6
S_3=raw.get_data(picks=[raw.ch_names.index('C4..')])*1e6
S_4=raw.get_data(picks=[raw.ch_names.index('Fc3.')])*1e6
S_5=raw.get_data(picks=[raw.ch_names.index('Fcz.')])*1e6
S_6=raw.get_data(picks=[raw.ch_names.index('Fc4.')])*1e6

# export paradigm
np.savetxt('S_events.dat', S_events, fmt='%d\t%d\t%d', delimiter='\n')
# export signals matrix
S_signals = np.transpose(np.concatenate((S_1,S_2,S_3,S_4,S_5,S_6)))
spio.savemat('S_signals.mat', mdict={'S_signals': S_signals})
