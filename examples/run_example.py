import araproc # noqa: F401
from araproc.framework import dataset 
from araproc.framework import data_visualization as dv

import logging
logging.getLogger().setLevel(logging.INFO)

# uncomment this dataset loade to see a nice CW contaminated A2 event
# a nice A2 dataset with CW (go after eventNumber 213179, which is TTree index 20695 in the burn sample file)
d = dataset.AraDataset(
    path_to_data_file="/data/exp/ARA/2013/filtered/burnSample1in10/ARA02/root/run1548/event1548.root",
    path_to_pedestal_file="/data/user/brianclark/ARA/framework_dev/software_deps/cvmfs/ara.opensciencegrid.org/trunk/ara_build/share/araCalib/ATRI/araAtriStation2Pedestals.txt"
)

# uncomment this dataset loader to see some nice A5 cal pulsers
# d = dataset.AraDataset(
#     path_to_data_file="/data/exp/ARA/2019/blinded/L1/ARA05/0701/run005626/event005626.root",
#     path_to_pedestal_file="/data/ana/ARA/ARA05PA/ARA05_pedestals/ped_full_values_A5_run005626.dat"
# )

for e in range(20690, 20700, 1):

    useful_event = d.get_useful_event(e)

    station2_case = (d.station_id==2) and (useful_event.eventNumber == 213179)
    station5_case = (d.station_id==5) and (useful_event.isCalpulserEvent())

    if station2_case or station5_case:
        print(e)
        
        # by default, you get all the bells and whistles
        # (interpolated, dedispersed, cw filtered, and bandpassed)
        wave_bundle = d.get_waveforms(useful_event)

        dv.plot_waveform_bundle(wave_bundle, 
                    time_or_freq="freq",
                    ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}.png",
                    )

        # we can also ask for just the interpolated traces, if we were really interested
        # (because we want to see the frequency domain, it's important we *at least*
        # as for interpolated traces)
        wave_bundle_interp_only = d.get_waveforms(useful_event, "interpolated")

        dv.plot_waveform_bundle(wave_bundle_interp_only, 
                            time_or_freq="freq",
                            ouput_file_path=f"./station_{d.station_id}_run_{d.run_number}_event_{useful_event.eventNumber}_interp_only.png",
                            )
