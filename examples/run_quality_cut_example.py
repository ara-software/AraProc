import argparse
from array import array
import numpy as np
import h5py
import araproc 
from tqdm import tqdm as tq
from araproc.framework import dataset
from araproc.analysis import standard_reco as sr
import araproc.analysis.daq_quality_cut as daq_qual_cut
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

# Read run number and total events
run_number = d.run_number
num_evts = d.num_events

# Setup output arrays
save_final_cuts = np.full((num_evts), 0, dtype=float)

#### investigate events with quality issue ######

# Process first 5 events
for e in range(5):
    useful_event = d.get_useful_event(e)
    evt_num = useful_event.eventNumber
    combined_errors = daq_qual_cut.check_daq_quality(useful_event, args.station, run_number)
    save_final_cuts[e] = combined_errors  # If this sum value is >0 for an event then it's a bad event 
    del useful_event,evt_num
    print(combined_errors)
# Save results to HDF5 file
file_name = f'daq_err_A{args.station}_run{run_number}.h5'
with h5py.File(file_name, 'w') as hf:
    hf.create_dataset('final_errors', data=save_final_cuts, compression="gzip", compression_opts=9)
