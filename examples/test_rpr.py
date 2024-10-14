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
save_rpr = np.full((16, num_evts), np.nan, dtype=float)
save_av_rpr = np.full(( num_evts), np.nan, dtype=float)
# Process first 5 events
for e in tq(range(5)):
    print(f'Processing event {e}')
    useful_event = d.get_useful_event(e)
    this_wave_bundle = d.get_waveforms(useful_event)
    rpr_arr ,average_rpr = rpr.get_avg_rpr(this_wave_bundle,individual_antenna = True)
    save_rpr[:,e] = rpr_arr
    save_av_rpr[e] = average_rpr

# Save results to HDF5 file
file_name = f'/data/ana/ARA/ARA0{args.station}/scratch/rpr2_A{args.station}_run{run_number}.h5'
with h5py.File(file_name, 'w') as hf:
    hf.create_dataset('per_channel_rpr', data=save_rpr, compression="gzip", compression_opts=9)
    hf.create_dataset('av_rpr', data=save_av_rpr, compression="gzip", compression_opts=9)
