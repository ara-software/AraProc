import argparse
import numpy as np
import h5py
from tqdm import tqdm as tq

from araproc.framework import dataset
from araproc.analysis import snr, rpr
from araproc.analysis import standard_reco as sr
from araproc.framework import waveform_utilities as wfu

# Argument parser
parser = argparse.ArgumentParser()

parser.add_argument("--input_file", type=str, required=True, help="full path to the input file")
parser.add_argument("--ped_file", type=str, default=None, required=False, help="path to pedestal file")
parser.add_argument("--is_simulation", type=int, default=0, required=True, help="is simulation; 0 = not simulation (default), 1 = is simulation")
parser.add_argument("--station", type=int, required=True, help="station id")

args = parser.parse_args()
args.is_simulation = bool(args.is_simulation)

# Raise errors for incompatible configurations
if args.is_simulation and (args.ped_file is not None):
    raise KeyError("You cannot mix a simulation with a pedestal file")
if (not args.is_simulation) and (args.ped_file is None):
    raise KeyError("If you are analyzing data, you must provide a pedestal file")

# Setup input
d = dataset.AnalysisDataset(
    station_id=args.station,
    path_to_data_file=args.input_file,
    path_to_pedestal_file=args.ped_file,
    is_simulation=args.is_simulation
)
reco = sr.StandardReco(d.station_id)
# Read run number and total events
run_number = d.run_number
num_evts = d.num_events

# Setup output arrays
save_snr = np.full((16, num_evts), np.nan, dtype=float)
save_rpr = np.full((16, num_evts), np.nan, dtype=float)  # New array to store rpr with same shape as snr

# Process first 5 events
for e in tq(range(num_evts)):
    print(f'Processing event {e}')
    useful_event = d.get_useful_event(e)
    this_wave_bundle = d.get_waveforms(useful_event)
    #print(useful_event)
    chans = None
    chans = list(this_wave_bundle.keys())
    for ch in chans:
        # Store the snr value for this channel and event
        save_snr[ch, e] = snr.get_snr(this_wave_bundle[ch])

        time, trace = wfu.tgraph_to_arrays(this_wave_bundle[ch])
        # Call the RPRCalculator for each channel
        rpr_calc = rpr.RPRCalculator()
        rpr_calc.run_rpr_calculation(trace, time, pad_num = len(trace))
        # Store the rpr value for this channel and event
        save_rpr[ch, e] = np.nanmax(rpr_calc.rpr_arr)  # Example: save the max RPR value, modify as needed
        print(save_rpr[ch, e])
        del time,trace
# Save results to HDF5 file
print('plotting ')
file_name = f'/data/ana/ARA/ARA0{args.station}/scratch/snr2_A{args.station}_run{run_number}.h5'
print('defined file',file_name)
with h5py.File(file_name, 'w') as hf:
    hf.create_dataset('snr', data=save_snr, compression="gzip", compression_opts=9)
    hf.create_dataset('rpr', data=save_rpr, compression="gzip", compression_opts=9)

