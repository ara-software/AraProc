import ctypes
import argparse
import ROOT
from array import array

import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.analysis import standard_reco as sr
from araproc.analysis import snr

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", type=str, required=True,
    help="full path to the input file")
parser.add_argument("--is_simulation", type=int, default=0, required=True,
	help="is simulation; 0 = not simulation (default), 1 = is simulation")
parser.add_argument("--station", type=int, required=True,
    help="station id")
parser.add_argument("--output_file", type=str, required=True,
    help="full path to the output file")

args = parser.parse_args()
args.is_simulation = bool(args.is_simulation)

# set up input 
d = dataset.AnalysisDataset(
    station_id = args.station,
    path_to_data_file = args.input_file,
    is_simulation = args.is_simulation
)
reco = sr.StandardReco(d.station_id)

# set up outputs
f = ROOT.TFile( args.output_file, 
               "RECREATE")
tree = ROOT.TTree("results_tree", "results_tree")

# this shows you how to store an integer
output_event_number = ctypes.c_int()
tree.Branch("event_number", output_event_number, "event_number/I")

# this shows you how to store a double
output_avg_snr = ctypes.c_double()
tree.Branch("avg_snr", output_avg_snr, "avg_snr/D")

# this shows you how to store an array
output_reco_result_pulser_v = array( "d", [0]*3 )
tree.Branch("reco_result_pulser_v", output_reco_result_pulser_v, "reco_result_pulser_v[3]/D")


for e in range(0, 5, 1):

    print(e)

    useful_event = d.get_useful_event(e)
    wavepacket = d.get_wavepacket(useful_event)
    this_wave_bundle = wavepacket["waveforms"]
    this_reco_results = reco.do_standard_reco(wavepacket)
    this_avg_snr = snr.get_avg_snr(this_wave_bundle)
    # stash the output results
    output_event_number.value = useful_event.eventNumber
    output_reco_result_pulser_v[0] = this_reco_results["pulser_v"]["corr"]
    output_reco_result_pulser_v[1] = this_reco_results["pulser_v"]["theta"]
    output_reco_result_pulser_v[2] = this_reco_results["pulser_v"]["phi"]
    output_avg_snr.value = this_avg_snr

    tree.Fill()

tree.Write()
f.Close()
