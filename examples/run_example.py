import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.analysis import standard_reco as sr
from araproc.framework import data_visualization as dv
from araproc.analysis import snr

import logging
logging.getLogger().setLevel(logging.INFO)

################################################################
################################################################
# uncomment this dataset loader to see some nice A5 cal pulsers
################################################################
################################################################

d = dataset.AnalysisDataset(
    station_id=5,
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    do_not_calibrate = False, # turn this flag to True if you only need raw atri events
)
iter_start = 0
iter_stop = 200

################################################################
################################################################
# uncomment this one to see some simulated events
################################################################
################################################################

# d = dataset.AnalysisDataset(
#     station_id = 4,
#     path_to_data_file="/data/ana/ARA/ARA04/sim_signal_and_noise_updated_121523/AraOut.setup_A4_C2.txt.run20.root",
#     is_simulation=True
# )
# iter_start = 0
# iter_stop = d.num_events

################################################################
################################################################
# uncomment this dataset loader to see a nice CW contaminated A2 event
# (eventNumber 213179, which is TTree index 20695 in the burn sample file)
################################################################
################################################################

# d = dataset.AnalysisDataset(
#     station_id = 2,
#     path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
# )
# iter_start = 20690
# iter_stop = 20700

################################################################
################################################################
# uncomment this one to see an example of CW removal using identified CW 
################################################################
################################################################

# d = dataset.AnalysisDataset(
#     station_id = 100,
#     path_to_data_file="/data/ana/ARA/processing/data/10pct/L2/station_100/year_2020/station_100_run_00019706.root",
#     path_to_cw_ids = "/home/brianclark/ARA/FiveStation/scripts/phase_var/CWID_station_100_run_19706.root",
# )
# iter_start = 0
# iter_stop = d.num_events

print(f"Config is {d.config}")
print(f"Excluded channels are {d.excluded_channels}")

# instantiate the standard reconstruction methods
reco = sr.StandardReco(d.station_id, 
                       excluded_channels=d.excluded_channels)

for e in range(iter_start, iter_stop, 1):

    if not d.is_simulation:
        # if this is not simulation, we can get raw events
        raw_event = d.get_raw_event(e)
    
    if not d.do_not_calibrate:
        # only if we have calibrated data
        useful_event = d.get_useful_event(e)
        
        if d.is_simulation:
            # if and only the dataset is simulation, we can ask for access
            # to some of the MC truth information 
            sim_info = d.get_event_sim_info(e)
            print(f"This event has weight {sim_info['weight']:.1e}, Energy {sim_info['enu']:.1e}")

        process_event = False
        if d.is_simulation:
            # for simulated events, I'm willing to plot everything
            process_event = True
        else:
            # but in the other case, I'd like to showcase either the CW contaminated event for A2
            # or just cal pulsers for A5
            station2_case = (d.station_id==2) and (useful_event.eventNumber == 213179)
            station5_case = (d.station_id==5) and (useful_event.isCalpulserEvent())
            process_event = station2_case or station5_case

        if process_event:
        
            # by default, you get all the bells and whistles
            # (interpolated, dedispersed, cw filtered, and bandpassed)
            wavepacket = d.get_wavepacket(useful_event)
            wave_bundle = wavepacket["waveforms"]

            # print the average snr across channels 
            avg_snr = snr.get_avg_snr(wave_bundle, excluded_channels=d.excluded_channels)
            print(f"The Average SNR is {avg_snr:.1f}")
    
            # run our standard suite of reconstructions
            reco_results = reco.do_standard_reco(wavepacket)

            pair_idx = reco.get_pair_index(1, 2, reco.pairs_v)

            # here's how we can lookup arrival times given a reconstructed direction
            arrival_time = reco.lookup_arrival_time(channel = 0, 
                                                    theta = reco_results["pulser_v"]["theta"], 
                                                    phi = reco_results["pulser_v"]["phi"],
                                                    which_distance="nearby",
                                                    solution=0,
                                                    )
            print(f"Arrival time at ch 0 is {arrival_time:.1f} ns")

            # plot the waveforms
            dv.plot_waveform_bundle(wave_bundle, 
                        time_or_freq="time",
                        output_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_waves.png",
                        )
            
            # and also plot one of the skymaps (in this case, the vpol map the distance of the pulser)
            dv.plot_skymap(reco_results["pulser_v"]["map"],
                        output_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_map.png",
                        )



