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
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)
iter_start = 0
iter_stop = 20700

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
################################################################
################################################################

# (eventNumber 213179, which is TTree index 20695 in the burn sample file)
# d = dataset.AnalysisDataset(
#     station_id = 2,
#     path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
#     path_to_pedestal_file="/cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/ara_build/share/araCalib/ATRI/araAtriStation2Pedestals.txt"
# )
# iter_start = 20690
# iter_stop = 20700


# instantiate the standard reconstruction methods
reco = sr.StandardReco(d.station_id)

for e in range(iter_start, iter_stop, 1):

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
        wave_bundle = d.get_waveforms(useful_event)

        # print the average snr across channels 
        print(snr.get_avg_snr(wave_bundle))      
 
        # run our standard suite of reconstructions
        reco_results = reco.do_standard_reco(wave_bundle)

        # plot the waveforms
        dv.plot_waveform_bundle(wave_bundle, 
                    time_or_freq="time",
                    ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_waves.png",
                    )
        
        # and also plot one of the skymaps (in this case, the vpol map the distance of the pulser)
        dv.plot_skymap(reco_results["pulser_v"]["map"],
                       ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{e}_map.png",
                       )
