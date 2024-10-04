import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.analysis import standard_reco as sr
from araproc.framework import data_visualization as dv

import logging
logging.getLogger().setLevel(logging.INFO)

# uncomment this dataset loader to see a nice CW contaminated A2 event
# (eventNumber 213179, which is TTree index 20695 in the burn sample file)
# d = dataset.AraDataset(
#     path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
#     path_to_pedestal_file="/cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/ara_build/share/araCalib/ATRI/araAtriStation2Pedestals.txt"
# )

# uncomment this dataset loader to see some nice A5 cal pulsers
d = dataset.AraDataset(
    path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
    path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
)

# instantiate the standard reconstruction methods
reco = sr.StandardReco(d.station_id)

for e in range(20690, 20700, 1):

    useful_event = d.get_useful_event(e)

    station2_case = (d.station_id==2) and (useful_event.eventNumber == 213179)
    station5_case = (d.station_id==5) and (useful_event.isCalpulserEvent())

    if station2_case or station5_case:
      
        # by default, you get all the bells and whistles
        # (interpolated, dedispersed, cw filtered, and bandpassed)
        wave_bundle = d.get_waveforms(useful_event, which_traces="filtered")
        
        # run our standard suite of reconstructions
        reco_results = reco.do_standard_reco(wave_bundle)

        # plot the waveforms
        dv.plot_waveform_bundle(wave_bundle, 
                    time_or_freq="time",
                    ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_waves.png",
                    )
        
        # and also plot one of the skymaps (in the case the vpol map at the distance of the pulser)
        dv.plot_skymap(reco_results["pulser_v"]["map"],
                       ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_map.png",
                       )
