import numpy as np

def get_list_from_txt(logs_file_name):

        bad_run = []
        logs_file = open(logs_file_name, "r")
        for lines in logs_file:
            try:
                run_num = int(lines.split('\t')[0])
            except ValueError:
                run_num = int(lines.split()[0])
            bad_run.append(run_num)
        bad_run = np.asarray(bad_run)
        logs_file.close()

        return bad_run

def get_bad_surface_run(station_id):

        # masked run(2014~2016) from brian's analysis
        # https://github.com/clark2668/a23_analysis_tools/blob/a7093ab2cbd6b743e603c23b9f296bf2bcce032f/tools_Cuts.h#L782
        # array for bad run
        bad_run = np.array([], dtype=int)

        if station_id == 2:

            # Runs shared with Ming-Yuan
            # http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889

            bad_run = np.append(bad_run, 2090)
            bad_run = np.append(bad_run, 2678)
            bad_run = np.append(bad_run, 4777)
            bad_run = np.append(bad_run, 5516)
            bad_run = np.append(bad_run, 5619)
            bad_run = np.append(bad_run, 5649)
            bad_run = np.append(bad_run, 5664)
            bad_run = np.append(bad_run, 5666)
            bad_run = np.append(bad_run, 5670)
            bad_run = np.append(bad_run, 5680)
            bad_run = np.append(bad_run, 6445)
            bad_run = np.append(bad_run, 6536)
            bad_run = np.append(bad_run, 6542)
            bad_run = np.append(bad_run, 6635)
            bad_run = np.append(bad_run, 6655)
            bad_run = np.append(bad_run, 6669)
            bad_run = np.append(bad_run, 6733)

            # Runs identified independently

            bad_run = np.append(bad_run, 2091)
            bad_run = np.append(bad_run, 2155)
            bad_run = np.append(bad_run, 2636)
            bad_run = np.append(bad_run, 2662)
            bad_run = np.append(bad_run, 2784)
            bad_run = np.append(bad_run, 4837)
            bad_run = np.append(bad_run, 4842)
            bad_run = np.append(bad_run, 5675)
            bad_run = np.append(bad_run, 5702)
            bad_run = np.append(bad_run, 6554)
            bad_run = np.append(bad_run, 6818)
            bad_run = np.append(bad_run, 6705)
            bad_run = np.append(bad_run, 8074)

        elif station_id == 3:

            # Runs shared with Ming-Yuan
            # http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=2041

            bad_run = np.append(bad_run, 977)
            bad_run = np.append(bad_run, 1240)
            bad_run = np.append(bad_run, 3158)
            bad_run = np.append(bad_run, 3431)
            bad_run = np.append(bad_run, 3432)
            bad_run = np.append(bad_run, 3435)
            bad_run = np.append(bad_run, 3437)
            bad_run = np.append(bad_run, 3438)
            bad_run = np.append(bad_run, 3439)
            bad_run = np.append(bad_run, 3440)
            bad_run = np.append(bad_run, 3651)
            bad_run = np.append(bad_run, 3841)
            bad_run = np.append(bad_run, 4472)
            bad_run = np.append(bad_run, 4963)
            bad_run = np.append(bad_run, 4988)
            bad_run = np.append(bad_run, 4989)

            # Runs identified independently

            bad_run = np.append(bad_run, 1745)
            bad_run = np.append(bad_run, 3157)
            bad_run = np.append(bad_run, 3652)
            bad_run = np.append(bad_run, 3800)
            bad_run = np.append(bad_run, 6193)
            bad_run = np.append(bad_run, 6319)
            bad_run = np.append(bad_run, 6426)

            # Runs I am sure we will exclude...

            bad_run = np.append(bad_run, 2000)
            bad_run = np.append(bad_run, 2001)

        else:
            pass

        return bad_run

def get_bad_run(station_id):

        # masked run(2014~2016) from brian's analysis
        # https://github.com/clark2668/a23_analysis_tools/blob/a7093ab2cbd6b743e603c23b9f296bf2bcce032f/tools_Cuts.h#L881

        # array for bad run
        bad_run = np.array([], dtype=int)

        if station_id == 2:

            ## 2013 ##
            # L1 data is already excluded the calibration run written in wiki
            bad_run = np.append(bad_run, 1918) # burn samples.....

            ## 2014 ##
            # 2014 rooftop pulsing, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, 3120)

            # 2014 surface pulsing
            # originally flagged by 2884, 2895, 2903, 2912, 2916
            # going to throw all runs jan 14-20
            bad_run = np.append(bad_run, np.arange(2884, 2918+1))
            """
            bad_run = np.append(bad_run, 2884) # jan 14 2014 surface pulser runs. actual problem causer
            bad_run = np.append(bad_run, [2885, 2889, 2890, 2891, 2893]) # exclusion by proximity

            bad_run = np.append(bad_run, 2895) # jan 16 2014 surface pulser runs. actual problem causer
            bad_run = np.append(bad_run, 2898) # exclusion by proximity
            bad_run = np.append(bad_run, [2900, 2901, 2902]) # jan 17 2014. exclusion by proximity

            bad_run = np.append(bad_run, 2903) # # jan 18 2014 surface pulser runs. actual problem causer
            bad_run = np.append(bad_run, [2905, 2906, 2907]) # exclusion by proximity

            bad_run = np.append(bad_run, 2912) # # jan 19 2014 surface pulser runs. actual problem causer
            bad_run = np.append(bad_run, 2915) # exclusion by proximity

            bad_run = np.append(bad_run, 2916) # jan 20 2014 surface pulser runs. actual problem causer
            bad_run = np.append(bad_run, 2918) # exclusion by proximity
            """

            # surface pulsing from m richman (identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14)
            bad_run = np.append(bad_run, [2938, 2939])

            # 2014 Cal pulser sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(3139, 3162+1))
            bad_run = np.append(bad_run, np.arange(3164, 3187+1))
            
            # 2014 rooftop pulsing
            bad_run = np.append(bad_run, 3242)

            # 2014 Cal pulser sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(3289, 3312+1))

            # ARA02 stopped sending data to radproc. Alert emails sent by radproc.
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            # http://ara.icecube.wisc.edu/wiki/index.php/Drop_29_3_2014_ara02
            bad_run = np.append(bad_run, 3336)

            # 2014 L2 Scaler Masking Issue.
            # Cal pulsers sysemtatically do not reconstruct correctly, rate is only 1 Hz
            # Excluded because configuration was not "science good"
            bad_run = np.append(bad_run, np.arange(3464, 3504+1))
            bad_run = np.append(bad_run, 3505)

            # 2014 Trigger Length Window Sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(3577, 3598+1))

            
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            # 2014, 4th June, Checking the functionality of the L1Scaler mask.
            bad_run = np.append(bad_run, 3695) # Masiking Ch0,1, 14
            bad_run = np.append(bad_run, 3700) # Masiking Ch2, 14
            bad_run = np.append(bad_run, 3701) # Masiking Ch4,5, 14
            bad_run = np.append(bad_run, 3702) # Masiking Ch6,7, 14
            bad_run = np.append(bad_run, 3703) # Masiking Ch8,9, 14
            bad_run = np.append(bad_run, 3704) # Masiking Ch10,11, 14
            bad_run = np.append(bad_run, 3705) # Masiking Ch12,13, 14
            bad_run = np.append(bad_run, 3706) # Masiking Ch14, 15

            # 2014, 16th June, Software update on ARA02 to fix the L1triggers.
            bad_run = np.append(bad_run, 3768)

            # 2014, 31st July, Testing new software to change trigger and readout window, pre-trigger samples.
            bad_run = np.append(bad_run, np.arange(3988, 3994+1))

            # 2014, 5th Aug, More tests on the pre-trigger samples.
            bad_run = np.append(bad_run, np.arange(4019, 4022+1))

            # 2014, 6th Aug, Switched to new readout window: 25 blocks, pre-trigger: 14 blocks.
            bad_run = np.append(bad_run, 4029)

            # 2014, 14th Aug, Finally changed trigger window size to 170ns.
            # http://ara.icecube.wisc.edu/wiki/index.php/File:Gmail_-_-Ara-c-_ARA_Operations_Meeting_Tomorrow_at_0900_CDT.pdf
            bad_run = np.append(bad_run, 4069)
            
            
            ## 2015 ##
            # ??
            bad_run = np.append(bad_run, 4004)

            # 2015 icecube deep pulsing
            # 4787 is the "planned" run
            # 4795,4797-4800 were accidental
            bad_run = np.append(bad_run, 4785) # accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001719/003/181001_ARA02AnalysisUpdate.pdf, slide 38)
            bad_run = np.append(bad_run, 4787) # deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
            bad_run = np.append(bad_run, np.arange(4795, 4800+1))

            # 2015 noise source tests, Jan, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
            bad_run = np.append(bad_run, np.arange(4820, 4825+1))
            bad_run = np.append(bad_run, np.arange(4850, 4854+1))
            bad_run = np.append(bad_run, np.arange(4879, 4936+1))
            bad_run = np.append(bad_run, np.arange(5210, 5277+1))

            # 2015 surface pulsing, Jan, http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
            bad_run = np.append(bad_run, [4872, 4873])
            bad_run = np.append(bad_run, 4876) # Identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14

            # 2015 Pulser Lift, Dec, http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 2)
            # Run number from private communication with John Kelley
            bad_run = np.append(bad_run, 6513)

            # 2015 ICL pulsing, Dec, http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
            bad_run = np.append(bad_run, 6527)


            ## 2016 ##
            # 2016 03/11/2016 Switched to D5 Vpol pulser for timing check.
            bad_run = np.append(bad_run, 7005)           
 
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
            # 2016, 21st July, Reduced trigger delay by 100ns.
            bad_run = np.append(bad_run, 7623)

            # 2016 cal pulser sweep, Jan 2015?, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
            bad_run = np.append(bad_run, np.arange(7625, 7686+1))

            ## other ##
            # D1 Glitches, Identified by MYL as having glitches after long periods of downtime
            bad_run = np.append(bad_run, 3)
            bad_run = np.append(bad_run, 11)
            bad_run = np.append(bad_run, 59)
            bad_run = np.append(bad_run, 60)
            bad_run = np.append(bad_run, 71)

            # Badly misreconstructing runs
            # run 8100. Loaded new firmware which contains the individual trigger delays which were lost since PCIE update in 12/2015.
            bad_run = np.append(bad_run, np.arange(8100, 8246+1))


            ## 2017 ##
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2017
            # 01/16/2017, Rooftop pulser run, Hpol ran for 30 min at 1 Hz starting 22:13:06. Vpol ran for 30 min at 1 Hz starting 22:44:50.
            bad_run = np.append(bad_run, 8530)

            # 01/24/2017, Deep pulser run, IC string 1 shallow pulser ~23:48-00:00. IC string 22 shallow pulser (Jan 25) ~00:01-00:19.
            bad_run = np.append(bad_run, 8573)

            # 01/25/2017, A2D6 pulser lift, Ran in continuous noise mode with V&Hpol Tx.
            bad_run = np.append(bad_run, [8574, 8575])

            # 01/25/2017, Same configuration as run8575, Ran in continuous noise mode with Hpol Tx. Forgot to switch back to normal configuration. No pulser lift in this period.
            bad_run = np.append(bad_run, [8576, 8578])

            # Cal pulser attenuation sweep
            """
            bad_run = np.append(bad_run, 8953) # 04/10/2017, Data runs incorrectly tagged as calibration (this is real data!)
            bad_run = np.append(bad_run, np.arange(8955, 8956+1)) # 04/10/2017, Data runs incorrectly tagged as calibration (this is real data!)
            bad_run = np.append(bad_run, np.arange(8958, 8962+1)) # 04/10/2017, Data runs incorrectly tagged as calibration (this is real data!)
            """
            bad_run = np.append(bad_run, np.arange(8963, 9053+1)) # 04/10/2017, System crashed on D5 (D6 completed successfully); D6 VPol 0 dB is 8963...D6 VPol 31 dB is 8974...D6 HPol 0 dB is 8975...D5 VPol 0 dB is 9007...crashed before D5 HPol
            bad_run = np.append(bad_run, np.arange(9129, 9160+1)) # 04/25/2017, D6 VPol: 9129 is 0 dB, 9130 is 1 dB, ... , 9160 is 31 dB
            bad_run = np.append(bad_run, np.arange(9185, 9216+1)) # 05/01/2017, D6 HPol: 9185 is 0 dB, 9186 is 1 dB, ... , 9216 is 31 dB
            bad_run = np.append(bad_run, np.arange(9231, 9262+1)) # 05/04/2017, D5 VPol: 9231 is 0 dB, ... , 9262 is 31 dB
            bad_run = np.append(bad_run, np.arange(9267, 9298+1)) # 05/05/2017, D5 HPol: 9267 is 0 dB, ... , 9298 is 31 dB


            ## 2018 ##
            # now using ARARunLogDataBase


            ## 2019 ##
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2019
            # D5 Calpulser sweep, 01/25/2019
            bad_run = np.append(bad_run, np.arange(12842, 12873+1)) # D5 Vpol attenuation sweep 0 to 31 dB with a step of 1 dB.
            bad_run = np.append(bad_run, np.arange(12874, 12905+1)) # D5 Hpol attenuation sweep 0 to 31 dB with a step of 1 dB. Wanted to verify if D5 Hpol actually fires or not. Conclusion was that D5 Hpol does not fire and ARA02 defaults to firing D5 Vpol instead.

            # D6 Vpol fired at 0 dB attenuation. Trigger delays of ARA2 ch. adjusted.
            # 03/22/2019 ~ 04/11/2019
            bad_run = np.append(bad_run, np.arange(13449, 13454+1))
            bad_run = np.append(bad_run, np.arange(13455, 13460+1))
            bad_run = np.append(bad_run, np.arange(13516, 13521+1))
            bad_run = np.append(bad_run, np.arange(13522, 13527+1))
            bad_run = np.append(bad_run, np.arange(13528, 13533+1))
            bad_run = np.append(bad_run, 13542)
            bad_run = np.append(bad_run, np.arange(13543, 13547+1))
            bad_run = np.append(bad_run, 13549)
            bad_run = np.append(bad_run, np.arange(13550, 13554+1))
            bad_run = np.append(bad_run, np.arange(13591, 13600+1))
            bad_run = np.append(bad_run, np.arange(13614, 13628+1))
            bad_run = np.append(bad_run, np.arange(13630, 13644+1))
            bad_run = np.append(bad_run, np.arange(13654, 13663+1))
            bad_run = np.append(bad_run, np.arange(13708, 13723+1))
            bad_run = np.append(bad_run, np.arange(13732, 13746+1))
            bad_run = np.append(bad_run, np.arange(13757, 13771+1))
            bad_run = np.append(bad_run, np.arange(13772, 13775+1))

            # Trigger delays of ARA2 ch.
            # 04/18/2019 ~ 05/2/2019
            bad_run = np.append(bad_run, np.arange(13850, 13875+1))
            bad_run = np.append(bad_run, np.arange(13897, 13898+1))
            bad_run = np.append(bad_run, np.arange(13900, 13927+1))
            bad_run = np.append(bad_run, np.arange(13967, 13968+1))
            bad_run = np.append(bad_run, np.arange(13970, 13980+1))
            bad_run = np.append(bad_run, np.arange(13990, 14004+1))
            bad_run = np.append(bad_run, np.arange(14013, 14038+1))
            bad_run = np.append(bad_run, np.arange(14049, 14053+1))
            bad_run = np.append(bad_run, np.arange(14055, 14060+1))
            bad_run = np.append(bad_run, np.arange(14079, 14087+1))
            bad_run = np.append(bad_run, np.arange(14097, 14105+1))
            bad_run = np.append(bad_run, np.arange(14115, 14123+1))
            bad_run = np.append(bad_run, np.arange(14133, 14141+1))
            bad_run = np.append(bad_run, np.arange(14160, 14185+1))
            bad_run = np.append(bad_run, np.arange(14194, 14219+1))
            bad_run = np.append(bad_run, np.arange(14229, 14237+1))

        elif station_id == 3:

            ## 2013 ##
            # Misc tests: http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2013
            # bad_run = np.append(bad_run, np.arange(22, 62+1))

            # ICL rooftop: http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
            # bad_run = np.append(bad_run, np.arange(63, 70+1))
            # bad_run = np.append(bad_run, np.arange(333, 341+1))

            # Cal sweep: http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
            # bad_run = np.append(bad_run, np.arange(72, 297+1))
            # bad_run = np.append(bad_run, np.arange(346, 473+1))

            # Eliminate all early data taking (all runs before 508)
            bad_run = np.append(bad_run, np.arange(508+1))

            bad_run = np.append(bad_run, 965) # burn samples.....
            bad_run = np.append(bad_run, 1652) # burn samples.....

            # Cal sweep: http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
            # ??

            ## 2014 ##
            # 2014 Rooftop Pulser, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, 2235)

            # 2014 Cal Pulser Sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(2251, 2274+1))

            # 2014 Rooftop Pulser
            bad_run = np.append(bad_run, 2328)

            # 2014 Cal Pulser Sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(2376, 2399+1))
            
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            # 2014, 6th Aug, Switched to new readout window: 25 blocks, pre-trigger: 14 blocks.
            bad_run = np.append(bad_run, 3063)

            # 2014, 14th Aug, Finally changed trigger window size to 170ns.
            bad_run = np.append(bad_run, 3103)
            
            ## 2015 ##
            # 2015 surface or deep pulsing
            # got through cuts
            # happened jan 5-6, some jan 8
            # waveforms clearly show double pulses or things consistent with surface pulsing
            bad_run = np.append(bad_run, 3811)
            bad_run = np.append(bad_run, [3810, 3820, 3821, 3822]) # elminated by proximity to deep pulser run
            bad_run = np.append(bad_run, 3823) # deep pulser, observation of 10% iterator event numbers 496, 518, 674, 985, 1729, 2411

            # 2015 noise source tests, Jan, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
            bad_run = np.append(bad_run, np.arange(3844, 3860+1))
            bad_run = np.append(bad_run, np.arange(3881, 3891+1))
            bad_run = np.append(bad_run, np.arange(3916, 3918+1))
            bad_run = np.append(bad_run, np.arange(3920, 3975+1))
            bad_run = np.append(bad_run, np.arange(4009, 4073+1))

            # 2015 surface pulsing, Jan, http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
            bad_run = np.append(bad_run, [3977, 3978])

            # 2015 ICL pulsing, Dec, http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
            bad_run = np.append(bad_run, 6041)

            # 2015 station anomaly
            # see moni report: http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1213
            # identified by MYL: http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
            bad_run = np.append(bad_run, np.arange(4914, 4960+1))

            ## 2016 ##
            # 2016 03/11/2016 Switched to D5 Vpol pulser for timing check.
            bad_run = np.append(bad_run, 6508)

            #http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
            # 2016, 21st July, Reduced trigger delay by 100ns.
            bad_run = np.append(bad_run, 7124)

            # More events with no RF/deep triggers, seems to precede coming test
            bad_run = np.append(bad_run, 7125)

            # 2016 Cal Pulser Sweep, http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
            bad_run = np.append(bad_run, np.arange(7126, 7253+1))
            
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
            # 2016 Loaded new firmware which contains the individual trigger delays which were lost since PCIE update in 12/2015.
            bad_run = np.append(bad_run, 7658)


            ## 2018 ##
            # now using ARARunLogDataBase

        elif station_id == 5:

            ## 2018 ##
            # http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2018
            # Calibration pulser lowered, http://ara.physics.wisc.edu/docs/0015/001589/002/ARA5CalPulser-drop-Jan-2018.xlsx

            # SPICEcore Run, http://ara.physics.wisc.edu/docs/0015/001589/002/SPICEcore-drop-log-8-Jan-2018.xlsx
            bad_run = np.append(bad_run)

        else:
            pass

        return bad_run

def get_L0_to_L1_Processing_run(station_id):

        # from http://ara.icecube.wisc.edu/wiki/index.php/Data_processing_and_storage_plan
        # array for bad run
        bad_run = np.array([], dtype=int)

        Year = np.arange(2014, 2020, dtype = int) #ASG: What years?
        Type = 'small'
        small_run_path = '/home/alansalgo/ARA/AraProc/araproc/analysis/data/Data_processing_and_storage_plan/'
       
        if station_id == 2 or station_id == 3: 
            for yrs in Year:
                small_run_file_name = f'a{station_id}_{int(yrs)}_{Type}_run.txt'
                small_run = np.loadtxt(small_run_path + small_run_file_name, dtype = int)    

                bad_run = np.append(bad_run, small_run)

        return bad_run

def get_ARARunLogDataBase(station_id):

        # https://github.com/yvpan/ARARunLogDataBase/tree/master/logs
        # array for bad run
        bad_run = np.array([], dtype=int)

        logs_path = '/home/alansalgo/ARA/AraProc/araproc/analysis/data/ARARunLogDataBase/'
        logs_file_name = f'{logs_path}a{station_id}_log.txt'

        if station_id == 2 or station_id == 3:            
            bad_run = np.append(bad_run, get_list_from_txt(logs_file_name))
            if station_id == 3:
                bad_run = bad_run[3:]
        else:
            pass

        return bad_run

def get_software_dominant_run(station_id):

        # https://github.com/clark2668/a23_analysis_tools/tree/master/data

        # array for bad run
        bad_run = np.array([], dtype=int)

        if station_id == 2:
            pass
        elif station_id == 3:
            bad_run = np.append(bad_run, [1795, 1796, 1797, 1799, 1800, 1801, 1802, 1804, 1805, 1806, 1807, 1809, 1810,
                                        1811, 1812, 1814, 1126, 1129, 1130, 1131, 1132, 1133, 1139, 1140, 1143, 1228, 1231,
                                        1322, 1428, 2000, 2037, 2038, 2042, 2043, 2466, 2467, 2468, 2469, 2471, 2472, 2473,
                                        3421, 3422, 3423, 3424, 3426, 3427, 3428, 3429, 3788, 3861, 3892, 3919, 4978, 5014,
                                        5024, 7756, 7757, 7758, 7760, 7761, 7762, 7763, 7765, 7766, 7767, 7768, 7770, 7771,
                                        7772, 7125, 7312, 7561, 7570])
        else:
            pass

        return bad_run

def get_obviously_bad_run(station_id):

         # array for bad run
        bad_run = np.array([], dtype=int)

        if station_id == 2:
            bad_run = np.append(bad_run, 2428) # 2428 2013/09/27 unknown signal
            bad_run = np.append(bad_run, 2631) # 2631 2013/11/06 unknown signal
            bad_run = np.append(bad_run, 2814) # 2814 2014/02/08 unknown signal
            bad_run = np.append(bad_run, 2849) # 2849 2014/02/10 unknown signal
            bad_run = np.append(bad_run, np.arange(2861, 2869 + 1, dtype = int)) # 2861 ~ 2869 2014/02/11 unknown signal
            bad_run = np.append(bad_run, 3442) # 3442 2014/04/19 unknown signal
            bad_run = np.append(bad_run, np.arange(3917, 3919 + 1, dtype = int)) # 2014/07/17 unknown signal
            bad_run = np.append(bad_run, 4107) # 4316 2014/08/22 unknown signal
            bad_run = np.append(bad_run, 4316) # 4316 2014/10/28 unknown signal
            bad_run = np.append(bad_run, 4645) # 4645 2014/12/08 unknown daq error
            bad_run = np.append(bad_run, 5911) # 5911 2015/07/21 unknown signal
            bad_run = np.append(bad_run, 6694) # 6694 2016/01/10 unknown signal
            bad_run = np.append(bad_run, 7100) # 7100 daq error
            bad_run = np.append(bad_run, 7405) # 7405 2016/06/06 temperature dependance
            bad_run = np.append(bad_run, 8518) # 8518 2017/01/12 unknown signal
            bad_run = np.append(bad_run, 8521) # 8521 2017/01/13 unknown signal
            bad_run = np.append(bad_run, np.arange(8523, 8527 + 1, dtype = int)) # 2017/01/15 unknown signal, temp issue
            bad_run = np.append(bad_run, 8545) # 8545 2017/01/19 unknown signal
            bad_run = np.append(bad_run, np.arange(8549, 8551 + 1, dtype = int)) # 8549 ~ 8551 2017/01/20 unknown signal
            bad_run = np.append(bad_run, np.arange(8558, 8574 + 1, dtype = int)) # 2017/01/22 ~ 3 only software, unknown signal
            bad_run = np.append(bad_run, np.arange(8591, 8596 + 1, dtype = int)) # 8591 ~ 8595 2017/01/28 unknown signal
            bad_run = np.append(bad_run, 8622) # 2017/02/03 unknown signal
            bad_run = np.append(bad_run, np.arange(8640, 8652 + 1, dtype = int)) # 8640 ~ 8652 2017/02/07 ~ 09 unknown signal
            bad_run = np.append(bad_run, np.arange(8680, 8687 + 1, dtype = int)) # 8680 ~ 8687 2017/02/14 ~ 16 unknown signal
            bad_run = np.append(bad_run, np.arange(8720, 8721 + 1, dtype = int)) # 2017/02/23 unknown signal
            bad_run = np.append(bad_run, 8763) # 8763 2017/03/01 unknown signal
            bad_run = np.append(bad_run, np.arange(8787, 8792 + 1, dtype = int)) # 2017/03/06 ~ 7 unknown signal
            bad_run = np.append(bad_run, 8821) # 2017/03/13 unknown signal
            bad_run = np.append(bad_run, np.arange(8933, 8934 + 1, dtype = int)) # 8933 ~ 8934 2017/04/05 unknown signal
            bad_run = np.append(bad_run, np.arange(9402, 9848 + 1, dtype = int)) # short... and pole season...
            bad_run = np.append(bad_run, 9916) # 9916 2018/02/07 unknown signal
            bad_run = np.append(bad_run, 9941) # 9941 bad noise model
            bad_run = np.append(bad_run, 10252) # 10252 2018/03/30 unknown signal
            bad_run = np.append(bad_run, np.arange(11071, 11076 + 1, dtype = int)) # 11071 ~ 11076 2018/06/25 unknown signal
            bad_run = np.append(bad_run, np.arange(11654, 11655 + 1, dtype = int)) # 11654 ~ 11655 2018/09/10 unknown signal
            bad_run = np.append(bad_run, np.arange(12128, 12131 + 1, dtype = int)) # 12131 2018/11/01 unknown signal
            bad_run = np.append(bad_run, 12248) # 12259 2018/11/15 unknown signal
            bad_run = np.append(bad_run, 12259) # 12259 2018/11/16 unknown signal
            bad_run = np.append(bad_run, 12345) # 12345 2018/11/23 unknown signal
            bad_run = np.append(bad_run, np.arange(12364, 12375 + 1, dtype = int)) # 12364 ~ 12375 2018/11/25 ~ 26 unknown signal
            bad_run = np.append(bad_run, 12383) # long!!!
            bad_run = np.append(bad_run, np.arange(12421, 12422 + 1, dtype = int)) # 12421 2018/12/05 unknow daq error
            bad_run = np.append(bad_run, np.arange(12429, 12430 + 1, dtype = int)) # 12430 2018/12/07 unknown daq error
            bad_run = np.append(bad_run, 12431) # short...
            bad_run = np.append(bad_run, np.arange(12435, 12436 + 1, dtype = int)) # 2018/12/07 unknown signal
            bad_run = np.append(bad_run, np.arange(12441, 12449 + 1, dtype = int)) # 2018/12/08 ~ 10 unknown signal
            bad_run = np.append(bad_run, np.arange(12463, 12465 + 1, dtype = int)) # 12464 2018/12/12 unknown signal
            bad_run = np.append(bad_run, np.arange(12470, 12474 + 1, dtype = int)) # 12470 ~ 12473 2018/12/14 unknown signal
            bad_run = np.append(bad_run, np.arange(12502, 12504 + 1, dtype = int)) # 12502 2018/12/17 unknown signal
            bad_run = np.append(bad_run, np.arange(12540, 12541 + 1, dtype = int)) # 12540 2018/12/21 unknown signal
            bad_run = np.append(bad_run, np.arange(12547, 12550 + 1, dtype = int)) # 12547 ~ 12548 2018/12/22 unknown signal
            bad_run = np.append(bad_run, np.arange(12560, 12567 + 1, dtype = int)) # 12560 ~ 12567 2018/12/24 unknown signal
            bad_run = np.append(bad_run, np.arange(12574, 12576 + 1, dtype = int)) # 12574 ~ 12576 2018/12/25 unknown signal
            bad_run = np.append(bad_run, np.arange(12583, 12598 + 1, dtype = int)) # 2018/12/26 ~ 27 unknown signal
            bad_run = np.append(bad_run, np.arange(12600, 12610 + 1, dtype = int)) # 12608 2018/12/28 unknown signal
            bad_run = np.append(bad_run, np.arange(12616, 12618 + 1, dtype = int)) # 12616 ~ 12618 2018/12/29 unknown signal
            bad_run = np.append(bad_run, np.arange(12633, 12647 + 1, dtype = int)) # 2018/12/31 ~ 01 unknown signal
            bad_run = np.append(bad_run, np.arange(12680, 12682 + 1, dtype = int)) # 2019/01/05 unknown signal
            bad_run = np.append(bad_run, np.arange(12691, 12692 + 1, dtype = int)) # 2019/01/06 unknown signal
            bad_run = np.append(bad_run, np.arange(12717, 12719 + 1, dtype = int)) # 2019/01/09 unknown signal
            bad_run = np.append(bad_run, np.arange(12725, 12727 + 1, dtype = int)) # 2019/01/10 unknown signal
            bad_run = np.append(bad_run, 12733) # 2019/01/11 unknown signal
            bad_run = np.append(bad_run, np.arange(12751, 12753 + 1, dtype = int)) # 12752 2019/01/13 unknown signal
            bad_run = np.append(bad_run, np.arange(12770, 12771 + 1, dtype = int)) # 12770 2019/01/17 unknown signal
            bad_run = np.append(bad_run, np.arange(12778, 12779 + 1, dtype = int)) # 12778 ~ 12779 2019/01/18 unknown signal
            bad_run = np.append(bad_run, np.arange(12835, 12905 + 1, dtype = int)) # short...
            bad_run = np.append(bad_run, np.arange(12922, 12923 + 1, dtype = int)) # 12922 ~ 12923 2019/01/27 unknown signal
            bad_run = np.append(bad_run, np.arange(12930, 12932 + 1, dtype = int)) # 12930 2019/01/28 unknown signal
            bad_run = np.append(bad_run, np.arange(12948, 12952 + 1, dtype = int)) # 12948 ~ 12952 2019/01/30 ~ 31 unknown signal
            bad_run = np.append(bad_run, 12961) # 12967 2019/02/01 unknown signal
            bad_run = np.append(bad_run, 12967) # 12967 2019/02/01 unknown signal
            bad_run = np.append(bad_run, 12988) # 12988 2019/02/04 unknown signal
            bad_run = np.append(bad_run, np.arange(12994, 12995 + 1, dtype = int)) # 12994 ~ 12995 2019/02/04 unknown signal
            bad_run = np.append(bad_run, np.arange(13001, 13002 + 1, dtype = int)) # 13001 2019/02/05 ev8850 unknown signal 
            bad_run = np.append(bad_run, 13253) # 13253 2019/03/02 unknown signal
            bad_run = np.append(bad_run, 13583) # 13583 2019/04/01 unknown signal
            bad_run = np.append(bad_run, 15246) # 15246 2019/07/16 unknown signal
            bad_run = np.append(bad_run, 15333) # 15333 2019/07/26 unknown signal
            bad_run = np.append(bad_run, np.arange(15340, 15343 + 1, dtype = int)) # 15340 ~ 15343 unknown signal
            bad_run = np.append(bad_run, 15587) # 15587 2020/01/13 unknown signal
            bad_run = np.append(bad_run, np.arange(15604, 15606 + 1, dtype = int)) # 15604 15605 15606 wrong l1 goal
            bad_run = np.append(bad_run, 16885) # 16885 wrong l1 goal

            bad_run = np.append(bad_run, np.array([4979, 5169, 5173, 5184, 5406, 5517, 5660, 5661, 5662, 5665, 5677, 5762, 5766, 5772, 5975], dtype = int)) # high event rate period

        elif station_id == 3:

            bad_run = np.append(bad_run, np.arange(515 + 1, dtype = int))
            bad_run = np.append(bad_run, np.arange(1124, 1144 + 1, dtype = int)) # short...
            #bad_run = np.append(bad_run, np.arange(1669, 1674 + 1, dtype = int)) # wrong l1 goal. It is from 2018 data....
            bad_run = np.append(bad_run, 1751) # 2013/11/08 unknown signal
            bad_run = np.append(bad_run, 1770) # unknown signal
            bad_run = np.append(bad_run, np.arange(1796, 1814 + 1, dtype = int)) # dda issue
            bad_run = np.append(bad_run, 1840) # unknown signal
            bad_run = np.append(bad_run, 1894) # unknown signal
            bad_run = np.append(bad_run, 2079) # unknown signal
            bad_run = np.append(bad_run, 2148) # wrong l1 goal
            bad_run = np.append(bad_run, 2628) # surface signal
            bad_run = np.append(bad_run, 3043) # dda issue
            bad_run = np.append(bad_run, np.arange(3843, 3861 + 1, dtype = int)) # noise source test
            bad_run = np.append(bad_run, np.arange(3881, 3892 + 1, dtype = int)) # noise source test
            bad_run = np.append(bad_run, np.arange(3916, 3975 + 1, dtype = int)) # noise source test
            bad_run = np.append(bad_run, 3982) # st1 spike
            bad_run = np.append(bad_run, np.arange(4008, 4073 + 1, dtype = int)) # noise source test
            bad_run = np.append(bad_run, 6270) # bad noise noise model
            bad_run = np.append(bad_run, np.arange(7122, 7153 + 1, dtype = int)) # 2016 Cal Pulser Sweep
            bad_run = np.append(bad_run, np.arange(10000, 10102 + 1, dtype = int)) # trim short runs
            bad_run = np.append(bad_run, np.arange(10158, 10160 + 1, dtype = int)) # 10158 2018/02/07 unknown signal
            bad_run = np.append(bad_run, 10167) # dda issue
            bad_run = np.append(bad_run, np.arange(10435, 10437 + 1, dtype = int)) # 10436 ~ 10437 2018/03/22 unknown signal
            bad_run = np.append(bad_run, np.arange(10659, 10666 + 1, dtype = int)) # dda issue
            bad_run = np.append(bad_run, 10684) # short...
            bad_run = np.append(bad_run, 10974) # 10974 2018/05/18 unknown signal
            bad_run = np.append(bad_run, np.arange(11108, 11110 + 1, dtype = int)) # 11108 ~ 11110 2018/06/01 ~ 2 unknown signal
            bad_run = np.append(bad_run, 11122) # short...
            bad_run = np.append(bad_run, 11295) # short...
            bad_run = np.append(bad_run, np.arange(11325, 11326 + 1, dtype = int)) # dda issue
            bad_run = np.append(bad_run, np.arange(11330, 11334 + 1, dtype = int)) # short...
            bad_run = np.append(bad_run, 11335) # possibly pulser
            bad_run = np.append(bad_run, np.arange(11428, 11432 + 1, dtype = int)) # 2018/07/04 ~ 5 unknown signal
            bad_run = np.append(bad_run, np.arange(12660, 12661 + 1, dtype = int)) # possibly pulser
            bad_run = np.append(bad_run, 12705) # 12705 2018/11/13 unknown daq error
            bad_run = np.append(bad_run, 12734) # 12734 2018/11/16 unknown signal
            bad_run = np.append(bad_run, np.arange(12884, 12885 + 1, dtype = int)) # possibly pulser
            bad_run = np.append(bad_run, np.arange(12903, 12905 + 1, dtype = int)) # 12903 ~ 12905 2018/12/12 unknown signal
            bad_run = np.append(bad_run, np.arange(12927, 12929 + 1, dtype = int)) # possibly pulser
            bad_run = np.append(bad_run, 13016) # 13016 2018/12/24 unknown signal 
            bad_run = np.append(bad_run, np.arange(13019, 13024 + 1, dtype = int)) # possibly spice core
            bad_run = np.append(bad_run, np.arange(13029, 13030 + 1, dtype = int)) # possibly spice core
            bad_run = np.append(bad_run, np.arange(13039, 13065 + 1, dtype = int)) # possibly spice core
            bad_run = np.append(bad_run, np.arange(13085, 13090 + 1, dtype = int)) # unknown signal
            bad_run = np.append(bad_run, np.arange(13113, 13117 + 1, dtype = int)) # 13113 ~ 13114 2019/01/03 unknown signal
            bad_run = np.append(bad_run, np.arange(13184, 13185 + 1, dtype = int)) # 13185 2019/01/11 unknown signal
            bad_run = np.append(bad_run, np.arange(13200, 13205 + 1, dtype = int)) # 13200 2019/01/13 unknown signal
            bad_run = np.append(bad_run, 13211) # 13211 2019/01/16 unknown signal
            bad_run = np.append(bad_run, 13265) # 2019/01/22 unknown signal
            bad_run = np.append(bad_run, np.arange(13272, 13274 + 1, dtype = int)) # 13272 ~ 13274 2019/01/23 unknown signal
            bad_run = np.append(bad_run, np.arange(13281, 13284 + 1, dtype = int)) # 13281 ~ 13284 2019/01/24 unknown signal
            bad_run = np.append(bad_run, np.arange(13290, 13291 + 1, dtype = int)) # 13291 2019/01/25 unknown signal
            bad_run = np.append(bad_run, 13333) # 13333 2019/01/30 unknown signal
            bad_run = np.append(bad_run, np.arange(13338, 13340 + 1, dtype = int)) # 13338 2019/01/30 unknown signal
            bad_run = np.append(bad_run, 13356) # 13356 2019/02/01 unknown signal
            bad_run = np.append(bad_run, np.arange(13443, 13444 + 1, dtype = int)) # 2019/02/11 unknown signal
            bad_run = np.append(bad_run, 13482) # 13482 2019/02/15 unknown signal
            bad_run = np.append(bad_run, 13482) # 2019/2/15 unknown signal
            bad_run = np.append(bad_run, 13995) # 2019/4/13 unknown signal
            bad_run = np.append(bad_run, np.arange(13964, 13969 + 1, dtype = int)) # low power...
            bad_run = np.append(bad_run, 16307) # unkown signal
            bad_run = np.append(bad_run, np.arange(16347, 16349 + 1, dtype = int)) # 16347 ~ 16348 2019/11/27 possible noise mode
            bad_run = np.append(bad_run, np.arange(16365, 16368 + 1, dtype = int)) # 16365 ~ 16366 2019/11/29 possible pulser signal
            bad_run = np.append(bad_run, np.arange(16374, 16376 + 1, dtype = int)) # obvious short run
            bad_run = np.append(bad_run, np.arange(16415, 16419 + 1, dtype = int)) # obvious short run
            bad_run = np.append(bad_run, np.arange(16450, 16451 + 1, dtype = int)) # 16450 ~ 16451 2019/12/08 possible pulser signal
            bad_run = np.append(bad_run, 16460) # 16460 2019/12/09 possible pulser signal
            bad_run = np.append(bad_run, np.arange(16518, 16520 + 1, dtype = int)) # unkown signal
            bad_run = np.append(bad_run, np.arange(16531, 16533 + 1, dtype = int)) # unkown signal
            bad_run = np.append(bad_run, np.arange(16538, 16540 + 1, dtype = int)) # unkown signal

            bad_run = np.append(bad_run, np.arange(924, 933 + 1, dtype = int)) # small Hpol runs
            bad_run = np.append(bad_run, np.array([6270, 13346, 13485, 13592, 13716, 13835, 14096, 14126, 14254, 15332, 16075, 16299, 
                                                    16358, 16790, 16919, 17310, 17311, 17312, 17313, 17315, 17317, 17318, 17320, 17410, 17800, 17948], dtype = int))
            bad_run = np.append(bad_run, np.array([16487, 16488, 16489, 16490, 16491, 16492, 16493, 16494, 16495, 16497, 16499, 16500,
                                                    16501, 16502, 16503, 16504, 16505, 16506, 16507, 16510, 16511, 16512, 16513, 16514,
                                                    16515, 16516, 16517, 16521, 16522, 16523, 16526, 16527, 16528, 16529, 16530, 16534,
                                                    16535, 16536, 16537], dtype = int)) # duplication ...urgh....

        else:
            pass

        return bad_run

def get_unchecked_unixtime(unix_time, station_id):

        bad_unix_time = False

        if station_id == 2:

            #if ((unix_time>=1364291997 and unix_time<=1364313598) or # 1605 2013/04/26 unknown signal
            if ((unix_time>=1365005690 and unix_time<=1365007171) or # 1653 2013/04/04 unknown signal
            (unix_time==1365158576) or # 1662 2013/04/05 ev56559 unknown signal
            (unix_time>=1365180309 and unix_time<=1365181221) or # 1663 2013/04/06 unknown signal
            (unix_time>=1365200930 and unix_time<=1365203643) or # 1664 2013/04/06 unknown signal
            (unix_time==1365240316) or # 1667 2013/04/06 ev7081 unknown signal
            (unix_time==1365267729) or # 1668 2013/04/07 ev68563 unknown signal
            (unix_time==1365278344) or # 1668 2013/04/07 ev180301 unknown signal
            (unix_time>=1365398088 and unix_time<=1365401803) or # 1675 2013/04/08 unknown signal
            (unix_time>=1365430776 and unix_time<=1365432390) or # 1677 2013/04/08 unknown signal
            (unix_time==1365610940) or # 1688 2013/04/11 ev38255 unknown signal 
            (unix_time==1365765373) or # 1697 2013/04/12 ev69677 unknown signal 
            (unix_time==1365772205) or # 1697 2013/04/12 ev140832 unknown signal 
            (unix_time>=1366128010 and unix_time<=1366128542) or # 1709 2013/04/17 unknown signal           
            (unix_time==1366329871) or # 1717 2013/04/19 ev56543 unknown signal
            (unix_time==1366387708) or # 1719 2013/04/19 ev209963 unknown signal
            (unix_time==1366389625) or # 1721 2013/04/20 ev2145 unknown signal
            (unix_time>=1367820783 and unix_time<=1367823667) or # 1760 2013/05/06 unknown signal
            (unix_time==1368153002) or # 1781 2013/05/10 ev67608 unknown signal
            (unix_time>=1368269973 and unix_time<=1368271966) or # 1788 2013/05/11 unknown signal
            (unix_time==1368976302) or # 1829 2013/05/19 ev76701 unknown signal    
            #(unix_time>=1368985684 and unix_time<=1369007282) or # 1830 2013/05/20 unknown signal
            #(unix_time>=1369014843 and unix_time<=1369017671) or # 1831 2013/05/20 unknown signal
            (unix_time==1370993889) or # 1902 2013/06/12 ev89323 unknown signal
            (unix_time==1371062135) or # 1905 2013/06/12 ev110714 unknown signal
            (unix_time==1378386757) or # 2320 2013/09/05 ev48411 unknown signal
            (unix_time==1379042233) or # 2359 2013/09/13 ev27058 unknown signal
            (unix_time==1379348598) or # 2377 2013/09/16 ev50469 unknown signal
            #(unix_time>=1380227664 and unix_time<=1380249257) or # 2428 2013/09/27 unknown signal
            (unix_time==1381889044) or # 2527 2013/10/16 ev75114 unknown signal
            (unix_time==1377413436) or # 2669 2013/08/25 ev120615 unknown signal
            (unix_time==1382919018) or # 2586 2013/10/28 ev94938 unknown signal
            #(unix_time>=1383683339 and unix_time<=1383704926) or # 2631 2013/11/06 unknown signal
            (unix_time==1385061252) or # 2711 2013/11/22 ev59463 unknown signal
            (unix_time==1385095045) or # 2712 2013/11/22 ev143445 unknown signal
            (unix_time==1388718619) or # 2829 2014/02/09 ev57693 unknown signal
            #(unix_time>=1389075034 and unix_time<=1389078167) or # 2849 2014/02/10 unknown signal 
            #(unix_time==1389277730) or # 2861 2014/02/11 ev34856 unknown signal
            #(unix_time>=1389316086 and unix_time<=1389337687) or # 2864 2014/02/11 unknown signal
            #(unix_time>=1389380987 and unix_time<=1389402587) or # 2868 2014/02/11 unknown signal
            #(unix_time==1389406778) or # 2869 2014/02/11 ev29309 unknown signal
            (unix_time>=1394225932 and unix_time<=1394247529) or # 3206 2014/03/08 unknown signal
            #(unix_time>=1397847956 and unix_time<=1397869555) or # 3442 2014/04/19 unknown signal
            (unix_time==1397959819) or # 3448 2014/04/20 ev29089 unknown signal
            (unix_time==1398013597) or # 3450 2014/04/20 ev110977 unknown signal
            (unix_time>=1402626427 and unix_time<=1402635831) or # 3750 2014/06/13 unknown signal
            (unix_time==1402732600) or # 3756 2014/06/14 ev94030 unknown signal       
            (unix_time>=1404275398 and unix_time<=1404276477) or # 3845 2014/07/02 unknown signal 
            (unix_time>=1405005392 and unix_time<=1405015777) or # 3886 2014/07/10 unknown signal 
            #(unix_time>=1405540485 and unix_time<=1405562082) or # 3917 2014/07/17 unknown signal 
            #(unix_time>=1405563124 and unix_time<=1405584720) or # 3919 2014/07/17 unknown signal 
            #(unix_time>=1412273315 and unix_time<=1412294914) or # 4316 2014/10/28 unknown signal 
            #(unix_time>=1418085579 and unix_time<=1418086324) or # 4645 2014/12/08 unknown daq error 
            (unix_time==1420546039) or # 4789 2015/01/06 ev29317 unknown signal           
            (unix_time>=1427755456 and unix_time<=1427757054) or # 5344 2015/03/30 unknown signal
            (unix_time==1427786641) or # 5345 2015/03/30 ev62457 unknown signal
            (unix_time>=1428163002 and unix_time<=1428184597) or # 5369 2015/04/04 unknown signal
            (unix_time==1428192620) or # 5370 2015/04/04 ev45146 unknown signal
            (unix_time==1435689453) or # 5805 2015/06/30 ev33371 unknown signal
            #(unix_time>=1437535952 and unix_time<=1437557546) or # 5911 2015/07/21 unknown signal
            (unix_time>=1452458202 and unix_time<=1452458308) or # 6694 2016/01/10 unknown signal
            (unix_time>=1453806013 and unix_time<=1453827605) or # 6770 2016/01/26 unknown signal
            (unix_time>=1456487612 and unix_time<=1456500604) or # 6933 2016/02/26 unknown signal
            #(unix_time>=1460489979 and unix_time<=1460511573) or # 7100 2016/04/12 unknown daq error
            (unix_time==1463141756) or # 7251 2016/05/13 ev89961 unknown signal  
            (unix_time>=1467446591 and unix_time<=1467468181) or # 7526 2016/07/02 unknown signal
            (unix_time>=1468892061 and unix_time<=1468897890) or # 7610 2016/07/18 unknown signal
            #(unix_time>=1469121769 and unix_time<=1469121836) or # 7624 2016/07/21 unknown signal
            (unix_time==1472429614) or # 7836 2016/08/28 ev66045 unknown signal
            (unix_time>=1472439514 and unix_time<=1472439944) or # 7836 2016/08/28 unknown signal
            (unix_time==1472507046) or # 7841 2016/08/29 ev7308 unknown signal
            (unix_time==1474313716) or # 7953 2016/09/19 ev104430 unknown signal
            #(unix_time>=1484277608 and unix_time<=1484277772) or # 8518 2017/01/12 unknown signal
            #(unix_time>=1484339916 and unix_time<=1484340121) or # 8521 2017/01/13 unknown signal
            #(unix_time>=1484515348 and unix_time<=1484517471) or # 8525 2017/01/15 unknown signal
            #(unix_time>=1484533404 and unix_time<=1484533604) or # 8526 2017/01/15 unknown signal
            #(unix_time==1484539534) or # 8527 2017/01/15 ev2905 unknown signal
            (unix_time==1484586356) or # 8529 2017/01/16 ev109945 unknown signal
            (unix_time==1484780668) or # 8540 2017/01/18 ev94709 unknown signal
            (unix_time==1484786214) or # 8540 2017/01/18 ev133130 unknown signal
            (unix_time==1484794191) or # 8541 2017/01/18 ev39070 unknown signal
            (unix_time==1484849482) or # 8544 2017/01/19 ev113601 unknown signal
            #(unix_time>=1484858727 and unix_time<=1484858793) or # 8545 2017/01/19 unknown signal
            #(unix_time==1484873486) or # 8545 2017/01/19 ev129213 unknown signal           
            #(unix_time>=1484920230 and unix_time<=1484941825) or # 8549 2017/01/20 unknown signal
            #(unix_time>=1484941834 and unix_time<=1484985027) or # 8550 ~ 8551 2017/01/20 unknown signal            
            #(unix_time>=1485127509 and unix_time<=1485166161) or # 8562 ~ 8563 2017/01/22 only software
            #(unix_time>=1485166183 and unix_time<=1485217628) or # 8564 ~ 8566 2017/01/23 unknown signal            
            #(unix_time>=1485217646 and unix_time<=1485239239) or # 8567 2017/01/23 only software
            #(unix_time>=1485239686 and unix_time<=1485261281) or # 8569 2017/01/24 unknown signal
            #(unix_time>=1485275601 and unix_time<=1485277009) or # 8570 2017/01/24 unknwon signal
            #(unix_time>=1485595456 and unix_time<=1485690018) or # 8591 ~ 8595 2017/01/28 unknown signal
            #(unix_time>=1485690030 and unix_time<=1485695000) or # 8596 2017/01/29 unknown signal
            (unix_time==1485752464) or # 8598 2017/01/29 ev128769 unknown signal
            (unix_time==1486068673) or # 8617 2017/02/02 ev51160 unknown signal
            #(unix_time>=1486168830 and unix_time<=1486190425) or # 8622 2017/02/03 unknown signal
            #(unix_time>=1486475432 and unix_time<=1486699107) or # 8640 ~ 8652 2017/02/07 ~ 09 unknown signal
            #(unix_time>=1487136780 and unix_time<=1487202291) or # 8680 ~ 8683 2017/02/14 ~ 15 unknown signal
            #(unix_time>=1487209220 and unix_time<=1487230815) or # 8684 2017/02/15 unknown signal
            #(unix_time==1487287677) or # 8687 2017/02/16 ev93962 unknown signal
            (unix_time==1487737935) or # 8713 2017/02/21 ev91749 unknown signal 
            (unix_time==1487828177) or # 8717 2017/02/22 ev125435 unknown signal
            #(unix_time>=1487854672 and unix_time<=1487876266) or # 8720 2017/02/23 unknown signal
            #(unix_time==1487891823) or # 8721 2017/02/23 ev105986 unknown signal 
            #(unix_time>=1488408765 and unix_time<=1488430359) or # 8763 2017/03/01 unknown signal
            (unix_time>=1488749304 and unix_time<=1488749444) or # 8782 2017/03/05 unknown signal
            (unix_time==1488789009) or # 8785 2017/03/06 ev20944 unknown signal
            (unix_time==1488803885) or # 8785 2017/03/06 ev123098 unknown signal
            (unix_time==1488819624) or # 8786 2017/03/06 ev82291 unknown signal
            #(unix_time>=1488829221 and unix_time<=1488850797) or # 8787 2017/03/06 unknown signal
            #(unix_time>=1488851843 and unix_time<=1488873438) or # 8789 2017/03/06 unknown signal            
            #(unix_time>=1488873446 and unix_time<=1488895038) or # 8790 2017/03/07 unknown signal
            #(unix_time>=1488916652 and unix_time<=1488938248) or # 8792 2017/03/07 unknown signal
            (unix_time>=1489226852 and unix_time<=1489240906) or # 8810 2017/03/11 unknown signal
            #(unix_time>=1489427707 and unix_time<=1489429596) or # 8821 2017/03/13 unknown signal
            (unix_time==1490012910) or # 8856 2017/03/20 ev4283 unknown signal           
            (unix_time>=1490481308 and unix_time<=1490481616) or # 8881 2017/03/25 unknown signal
            #(unix_time==1491385142) or # 8933 2017/04/05 ev45771 unknown signal 
            #(unix_time>=1491398237 and unix_time<=1491404695) or # 8933 ~ 8934 2017/04/05 unknown signal
            (unix_time>=1491636026 and unix_time<=1491642052) or # 8947 2017/04/07 unknown signal
            (unix_time>=1494298072 and unix_time<=1494310350) or # 9315 2017/05/08 unknown signal
            (unix_time==1495264665) or # 9369 2017/05/19 ev78691 unknown signal
            (unix_time>=1495282918 and unix_time<=1495283656) or # 9371 2017/05/20 unknown signal
            #(unix_time>=1514952103 and unix_time<=1514973698) or # 9511 2018/01/02 unknown signal
            #(unix_time>=1515057615 and unix_time<=1515058703) or # 9516 2018/01/04 only software
            #(unix_time>=1515279938 and unix_time<=1515301534) or # 9530 2018/01/06 unknown signal
            #(unix_time>=1515352635 and unix_time<=1515356375) or # 9557 ~ 9559 2018/01/07 unknown signal
            #(unix_time>=1516135871 and unix_time<=1516140648) or # 9755 ~ 9759 2018/01/16 unknown signal 
            #(unix_time==1516399622) or # 9813 2018/01/19 ev534 unknown signal
            #(unix_time>=1518035003 and unix_time<=1518056598) or # 9916 2018/02/07 unknown signal
            (unix_time>=1519366238 and unix_time<=1519368402) or # 10001 2018/02/22 unknown signal
            (unix_time==1521786112) or # 10185 2018/03/22 ev36449 unknown signal
            (unix_time>=1528963594 and unix_time<=1528967442) or # 10948 2018/06/14 unknown signal
            #(unix_time>=1529919692 and unix_time<=1529967349) or # 11071 ~ 11076 2018/06/25 unknown signal
            (unix_time>=1530292340 and unix_time<=1530294903) or # 11109 2018/06/29 unknown signal
            (unix_time>=1530432024 and unix_time<=1530442821) or # 11124 2018/07/01 unknown signal
            (unix_time>=1532659100 and unix_time<=1532662001) or # 11214 2018/07/26 unknown signal
            (unix_time>=1536270217 and unix_time<=1536270417) or # 11606 2018/09/06 unknown signal
            (unix_time==1536594391) or # 11651 2018/09/10 ev46296 unknown signal
            #(unix_time>=1536615278 and unix_time<=1536620571) or # 11654 ~ 11655 2018/09/10 unknown signal
            (unix_time>=1537528259 and unix_time<=1537529685) or # 11750 2018/09/21 unknown signal
            (unix_time>=1537653570 and unix_time<=1537654170) or # 11765 2018/09/22 unknown signal
            (unix_time==1537817429) or # 11781 2018/09/24 ev54262 unknown signal
            (unix_time>=1540683879 and unix_time<=1540686065) or # 12084 2018/10/27 unknown signal
            (unix_time>=1540879759 and unix_time<=1540882966) or # 12104 2018/10/29 unknown signal
            #(unix_time>=1541133032 and unix_time<=1541143829) or # 12131 2018/11/01 unknown signal
            #(unix_time>=1542394856 and unix_time<=1542405653) or # 12259 2018/11/16 unknown signal
            #(unix_time>=1543000089 and unix_time<=1543000501) or # 12345 2018/11/23 unknown signal
            #(unix_time==1543004788) or # 12345 2018/11/23 ev39161 unknown signal 
            #(unix_time>=1543174423 and unix_time<=1543193164) or # 12364 ~ 12365 2018/11/25 unknown signal
            #(unix_time>=1543272737 and unix_time<=1543281657) or # 12374 ~ 12375 2018/11/26 unknown signal
            #(unix_time>=1544060237 and unix_time<=1544061739) or # 12421 2018/12/05 unknow daq error
            #(unix_time==1544173852) or # 12430 2018/12/07 ev2921 unknown signal
            #(unix_time>=1544178891 and unix_time<=1544179107) or # 12430 2018/12/07 unknown daq error
            #(unix_time>=1544216543 and unix_time<=1544217194) or # 12435 2018/12/07 unknown signal
            #(unix_time>=1544287259 and unix_time<=1544315436) or # 12441 ~ 12443 2018/12/08 unknown signal
            #(unix_time>=1544473656 and unix_time<=1544512049) or # 12446 ~ 12449 2018/12/10 unknown signal
            #(unix_time>=1544654007 and unix_time<=1544654489) or # 12464 2018/12/12 unknown signal
            #(unix_time>=1544789738 and unix_time<=1544790291) or # 12470 2018/12/14 unknown signal
            #(unix_time>=1544814958 and unix_time<=1544825752) or # 12473 2018/12/14 unknown signal
            #(unix_time>=1545096700 and unix_time<=1545107496) or # 12502 2018/12/17 unknown signal
            #(unix_time>=1545433692 and unix_time<=1545433902) or # 12540 2018/12/21 unknown signal
            #(unix_time>=1545441082 and unix_time<=1545442425) or # 12540 2018/12/21 unknown signal
            #(unix_time>=1545508000 and unix_time<=1545519104) or # 12547 ~ 12548 2018/12/22 unknown signal
            #(unix_time>=1545763692 and unix_time<=1545786218) or # 12574 ~ 12576 2018/12/25 unknown signal
            #(unix_time>=1545862246 and unix_time<=1545870745) or # 12585 2018/12/26 unknown signal
            #(unix_time>=1545881645 and unix_time<=1545892246) or # 12592 2018/12/26 unknown signal
            #(unix_time>=1545941091 and unix_time<=1545950845) or # 12598 2018/12/27 unknown signal
            #(unix_time>=1546033422 and unix_time<=1546040180) or # 12608 2018/12/28 unknown signal
            #(unix_time>=1546112069 and unix_time<=1546144477) or # 12616 ~ 12618 2018/12/29 unknown signal
            #(unix_time>=1546144477 and unix_time<=1546295392) or # 12633 ~ 12634 2018/12/31 unknown signal
            #(unix_time>=1546317014 and unix_time<=1546327811) or # 12637 2018/12/31 unknown signal
            #(unix_time==1544062179) or # 12422 2018/12/05 ev3809 unknown signal
            #(unix_time>=1546372078 and unix_time<=1546393678) or # 12643 ~ 12644 2019/01/01 unknown signal
            (unix_time==1546527774) or # 12660 2019/01/03 ev2268 unknown signal
            #(unix_time>=1546723827 and unix_time<=1546734624) or # 12680 2019/01/01 unknown signal
            #(unix_time==1546743335) or # 12681 2019/01/05 ev48667 unknown signal
            #(unix_time>=1546745437 and unix_time<=1546756234) or # 12682 2019/01/05 unknown signal
            #(unix_time>=1546832915 and unix_time<=1546843712) or # 12691 2019/01/06 unknown signal
            #(unix_time>=1547078029 and unix_time<=1547088827) or # 12717 2019/01/09 unknown signal
            #(unix_time>=1547154704 and unix_time<=1547165498) or # 12725 2019/01/10 unknown signal
            #(unix_time>=1547176314 and unix_time<=1547187110) or # 12727 2019/01/10 unknown signal
            (unix_time==1547238967) or # 12733 2019/01/11 ev42737 unknown signal
            #(unix_time>=1547409505 and unix_time<=1547417218) or # 12752 2019/01/13 unknown signal
            #(unix_time>=1547778399 and unix_time<=1547779000) or # 12770 2019/01/17 unknown signal
            #(unix_time>=1547859028 and unix_time<=1547868920) or # 12778 ~ 12779 2019/01/18 unknown signal
            #(unix_time>=1548619137 and unix_time<=1548640740) or # 12922 ~ 12923 2019/01/27 unknown signal
            #(unix_time>=1548700875 and unix_time<=1548706609) or # 12930 2019/01/28 unknown signal
            #(unix_time>=1548870785 and unix_time<=1548903191) or # 12948 ~ 12950 2019/01/30 unknown signal
            #(unix_time>=1549064973 and unix_time<=1549067340) or # 12967 2019/02/01 unknown signal
            #(unix_time>=1549324753 and unix_time<=1549325527) or # 12994 2019/02/04 unknown signal
            #(unix_time==1549392031) or # 13001 2019/02/05 ev8850 unknown signal
            (unix_time==1550713678) or # 13167 2019/02/20 ev24484 unknown signal
            #(unix_time==1551536370) or # 13253 2019/03/02 ev6148 unknown signal
            #(unix_time>=1554112083 and unix_time<=1554112319) or # 13583 2019/04/01 unknown signal
            (unix_time==1554558561) or # 13685 2019/04/06 ev2205 unknown signal
            (unix_time==1554616905) or # 13691 2019/04/06 ev20905 unknown signal
            (unix_time==1554726152) or # 13702 2019/04/08 ev30542 unknown signal
            (unix_time==1555253840) or # 13806 2019/04/14 ev33746 unknown signal
            (unix_time==1555307059) or # 13811 2019/04/14 ev29231 unknown signal
            (unix_time==1555330618) or # 13814 2019/04/15 ev34245 unknown signal
            #(unix_time>=1560373936 and unix_time<=1560374098) or # 14805 2019/06/12 unknown signal
            #(unix_time>=1563301523 and unix_time<=1563312320) or # 15246 2019/07/16 unknown signal
            #(unix_time==1564205064) or # 15340 2019/07/26 ev14551 unknown signal
            #(unix_time==1564217347) or # 15341 2019/07/27 ev21787 unknown signal
            #(unix_time==1564234689) or # 15343 2019/07/27 ev51097 unknown signal
            (unix_time>=1564787556 and unix_time<=1564792632) or # 15401 2019/08/02 unknown signal
            (unix_time>=1564814402 and unix_time<=1564819157)): # 15404 2019/08/02 unknown signal

                bad_unix_time = True
        
        elif station_id == 3:

            if ((unix_time==1368588355) or # 805 2013/05/15 ev95159 unknown signal
            (unix_time>=1369906668 and unix_time<=1369906672) or # 892 2013/05/30 unknown signal
            (unix_time==1373818145) or # 1115 2013/07/14 ev76709 unknown signal
            (unix_time==1374308750) or # 1164 2013/07/20 ev68655 unknown signal
            (unix_time==1374399536) or # 1169 2013/07/21 ev97563 unknown signal
            (unix_time==1375584751) or # 1239 2013/08/04 ev84081 unknown signal
            (unix_time==1375631381) or # 1241 2013/08/04 ev107265 unknown signal
            (unix_time==1375717169) or # 1246 2013/08/05 ev102055 unknown signal
            #(unix_time>=1383884594 and unix_time<=1383886064) or # 1751 2013/11/08 unknown signal
            (unix_time>=1396880478 and unix_time<=1396888693) or # 2473 2014/04/07 remaining bias voltage
            (unix_time==1417740758) or # 3650 2014/12/05 ev140802 unknown signal
            (unix_time==1462176835) or # 6696 2016/05/01 ev89505 unknown signal
            (unix_time==1462203893) or # 6697 2016/05/02 ev131864 unknown signal
            #(unix_time>=1516349134 and unix_time<=1516349760) or # 10034 2018/01/19 unknown hpol signal
            #(unix_time>=1518047608 and unix_time<=1518049306) or # 10158 2018/02/07 unknown signal 
            #(unix_time>=1521727312 and unix_time<=1521748910) or # 10436 ~ 10437 2018/03/22 unknown signal
            (unix_time==1521867138) or # 10453 2018/03/23 ev57689 unknown signal
            #(unix_time>=1526664593 and unix_time<=1526665370) or # 10974 2018/05/18 unknown signal
            #(unix_time==1527925450) or # 11108 2018/06/01 ev57654 unknown signal
            #(unix_time>=1527926151 and unix_time<=1527943352) or # 11109 ~ 11110 2018/06/02 unknown signal
            #(unix_time>=1530765991 and unix_time<=1530776786) or # 11428 2018/07/04 unknown signal
            #(unix_time>=1530777833 and unix_time<=1530788625) or # 11430 2018/07/05 unknown signal
            #(unix_time>=1530788637 and unix_time<=1530789085) or # 11431 2018/07/05 unknown signal
            #(unix_time>=1530799442 and unix_time<=1530810237) or # 11432 2018/07/05 unknown signal
            (unix_time==1534243026) or # 11802 2018/08/14 evt7 possible untagged calpulser
            (unix_time>=1538405062 and unix_time<=1538405098) or # 12252 2018/10/01 unknown signal
            #(unix_time>=1542394027 and unix_time<=1542404823) or # 12734 2018/11/16 unknown signal
            (unix_time==1544579661) or # 12896 2018/12/11 evt10965 possible untagged calpulser
            #(unix_time>=1544651797 and unix_time<=1544674421) or # 12903 ~ 12905 2018/12/12 unknown signal
            #(unix_time>=1545944323 and unix_time<=1545947134) or # 13049 ~ 13050 2018/12/16 unknown signal
            #(unix_time>=1546022641 and unix_time<=1546025431) or # 13059 2018/12/28 unknown signal
            #(unix_time>=1546276704 and unix_time<=1546287501) or # 13086 2018/12/31 possible noise mode
            #(unix_time>=1546544447 and unix_time<=1546566048) or # 13113 ~ 13114 2019/01/03 unknown signal
            #(unix_time>=1547270034 and unix_time<=1547271615) or # 13185 2019/01/11 unknown signal
            #(unix_time>=1547416408 and unix_time<=1547416572) or # 13200 2019/01/13 unknown signal
            #(unix_time>=1547674833 and unix_time<=1547685628) or # 13211 2019/01/16 unknown signal
            (unix_time>=1547780392 and unix_time<=1547780462) or # 13221 2019/01/17 unknown signal
            #(unix_time>=1548202563 and unix_time<=1548213359) or # 13265 2019/01/22 unknown signal
            #(unix_time>=1548285602 and unix_time<=1548300827) or # 13273 ~ 13274 2019/01/23 unknown signal
            #(unix_time>=1548358209 and unix_time<=1548377141) or # 13282 ~ 13283 2019/01/24 unknown signal
            #(unix_time>=1548382029 and unix_time<=1548382641) or # 13284 2019/01/24 unknown signal
            #(unix_time>=1548444038 and unix_time<=1548454827) or # 13291 2019/01/25 unknown signal
            #(unix_time>=1548903559 and unix_time<=1548903984) or # 13338 2019/01/30 unknown signal
            #(unix_time>=1549067558 and unix_time<=1549067752) or # 13356 2019/02/01 unknown signal
            #(unix_time>=1549916886 and unix_time<=1549918464) or # 13443 ~ 13444 2019/02/11 unknown signal
            #(unix_time>=1550288014 and unix_time<=1550288043) or # 13482 2019/02/15 unknown signal
            (unix_time>=1571810125 and unix_time<=1571813629)): # 16020 ~ 16021 2019/10/22 ~ 23 possible surface activity
            #(unix_time>=1574486360 and unix_time<=1574497155) or # 16307 2019/11/22 possible surface activity
            #(unix_time>=1574883653 and unix_time<=1574898070) or # 16347 ~ 16348 2019/11/27 possible noise mode
            #(unix_time>=1575062329 and unix_time<=1575073417) or # 16365 ~ 16366 2019/11/29 possible pulser signal
            #(unix_time>=1575824828 and unix_time<=1575839953) or # 16450 ~ 16451 2019/12/08 possible pulser signal
            #(unix_time>=1575916845 and unix_time<=1575927640) or # 16460 2019/12/09 possible pulser signal
            #(unix_time>=1576503639 and unix_time<=1576533683)): # 16518 ~ 16520 2019/12/16 possible pulser signal

                bad_unix_time = True

        elif station_id == 5:
            pass

        return bad_unix_time

def get_bad_unixtime(unix_time, station_id):

        # masked unixtime(2014~2016) from brian's analysis
        # https://github.com/clark2668/a23_analysis_tools/blob/a7093ab2cbd6b743e603c23b9f296bf2bcce032f/tools_Cuts.h#L503

        bad_unix_time = False

        if station_id == 2:

            # Livetime flagged as bad by Biran
            if((unix_time>=1389381600 and unix_time<=1389384000) or # from run 2868
            (unix_time>=1420317600 and unix_time<=1420318200) or # from run 4775
            # (unix_time>=1449189600 and unix_time<=1449190200) or # from run 6507
            (unix_time>=1449187200 and unix_time<=1449196200) or # from run 6507

            #Livetime flagged as bad by Biran's undergrads
            #config 1
            # (unix_time>=1380234000 and unix_time<=1380236400) or # from run 2428 22 hour balloon launch
            # (unix_time>=1382046000 and unix_time<=1382047500) or # from run 2536 22 hour balloon launch
            (unix_time>=1382712900 and unix_time<=1382713500) or # from run 2575
            (unix_time>=1382972700 and unix_time<=1382973300) or # from run 2589
            # (unix_time>=1383689100 and unix_time<=1383690900) or # from run 2631 22 hour balloon launch
            (unix_time>=1383884400 and unix_time<=1383886200) or # from run 2642
            (unix_time>=1384060200 and unix_time<=1384061100) or # from run 2652
            (unix_time>=1384487400 and unix_time<=1384489800) or # from run 2677
            (unix_time>=1384489980 and unix_time<=1384491060) or # from run 2678 at start may be glitch or continued from 2677
            (unix_time>=1384856520 and unix_time<=1384856640) or # from run 2698 super zoomed in two minute window
            # (unix_time>=1385674200 and unix_time<=1385675100) or # from run 2744 22 hour balloon launch
            (unix_time>=1389381600 and unix_time<=1389383700) or # from run 2868 first of two from run 2868
            (unix_time>=1389398700 and unix_time<=1389400200) or # from run 2868 second of two from run 2868
            (unix_time>=1389665100 and unix_time<=1389666300) or # from run 2884
            (unix_time>=1393288800 and unix_time<=1393289400) or # from run 3099
            # (unix_time>=1397856600 and unix_time<=1397858400) or # from run 3442 22 hour balloon launch

            #config 2
            (unix_time>=1376731800 and unix_time<=1376733000) or # from run 2235

            #conifg 3
            (unix_time>=1400276700 and unix_time<=1400277300) or # from run 3605 mainly looks like glitch at end

            #config 4
            (unix_time>=1409986500 and unix_time<=1409988000) or # from run 4184
            # (unix_time>=1412026200 and unix_time<=1412027100) or # from run 4301 22 hr balloon
            # (unix_time>=1412285400 and unix_time<=1412288100) or # from run 4316 weird 22hr balloon
            # (unix_time>=1412544600 and unix_time<=1412545500) or # from run 4331 22hr balloon
            # (unix_time>=1412803800 and unix_time<=1412804700) or # from run 4346 22hr balloon
            (unix_time>=1413898200 and unix_time<=1413899100) or # from run 4408
            (unix_time>=1414083900 and unix_time<=1414086000) or # from run 4418
            (unix_time>=1414350300 and unix_time<=1414351200) or # from run 4434 pt 1
            # (unix_time>=1414358700 and unix_time<=1414359900) or # from run 4434 pt 2 22hr balloon
            (unix_time>=1414674300 and unix_time<=1414674780) or # from run 4452
            (unix_time>=1414986600 and unix_time<=1414987200) or # from run 4471
            (unix_time>=1415223000 and unix_time<=1415223900) or # from run 4483
            (unix_time>=1415380500 and unix_time<=1415381400) or # from run 4493
            (unix_time>=1415558100 and unix_time<=1415559000) or # from run 4503
            (unix_time>=1415742300 and unix_time<=1415743800) or # from run 4513
            (unix_time>=1416207000 and unix_time<=1416212100) or # from run 4541
            (unix_time>=1420978200 and unix_time<=1420978800) or # from run 4814
            (unix_time>=1416905100 and unix_time<=1416910500) or # from run 4579 two spikes about an hour apart
            # (unix_time>=1416950700 and unix_time<=1416951600) or # from run 4582 22 hour balloon launch
            (unix_time>=1417677000 and unix_time<=1417678200) or # from run 4621  weird and cool
            (unix_time>=1417836000 and unix_time<=1417837500) or # from run 4631
            (unix_time>=1420097100 and unix_time<=1420098300) or # from run 4763
            (unix_time>=1420293300 and unix_time<=1420294200) or # from run 4774
            (unix_time>=1420317600 and unix_time<=1420318200) or # from run 4775
            (unix_time>=1420978200 and unix_time<=1420978800) or # from run 4814
            (unix_time>=1421024400 and unix_time<=1421025300) or # from run 4817
            (unix_time>=1421713200 and unix_time<=1421718600) or # from run 4872 looks full of errors and not spiky but could have a spiky
            (unix_time>=1421718000 and unix_time<=1421725800) or # from run 4873 definitely an error but also has spiky boy, part 1 of 2
            (unix_time>=1421733300 and unix_time<=1421733900) or # from run 4873 spiky boy alone but in a run with errors, part 2 of 2
            (unix_time>=1421783400 and unix_time<=1421794200) or # from run 4876 definitely an error but not a spikey boy
            # (unix_time>=1428529800 and unix_time<=1428530700) or # from run 5389 22 hour balloon launch
            (unix_time>=1435623000 and unix_time<=1435623600) or # from run 5801
            # (unix_time>=1436394000 and unix_time<=1436395200) or # from run 5845 22 hour balloon launch
            (unix_time>=1437601200 and unix_time<=1437602700) or # from run 5915 looks like error at the start
            # (unix_time>=1439933700 and unix_time<=1439934960) or # from run 6048 22 hour balloon launch
            (unix_time>=1440581700 and unix_time<=1440582480) or # from run 6086
            # (unix_time>=1441489200 and unix_time<=1441490280) or # from run 6137 22 hour balloon launch
            # (unix_time>=1444685400 and unix_time<=1444687080) or # from run 6322 22 hour balloon launch
            # (unix_time>=1445722020 and unix_time<=1445723220) or # from run 6383 22 hour balloon launch
            (unix_time>=1445934900 and unix_time<=1445935500) or # from run 6396
            (unix_time>=1445960400 and unix_time<=1445961000) or # from run 6397
            # (unix_time>=1445982120 and unix_time<=1445982900) or # from run 6398 22 hour balloon launch
            (unix_time>=1446165600 and unix_time<=1446166200) or # from run 6408
            # (unix_time>=1446327300 and unix_time<=1446328200) or # from run 6418 22 hour balloon launch
            (unix_time>=1446607800 and unix_time<=1446608640) or # from run 6433 looks like an error at end
            (unix_time>=1446784200 and unix_time<=1446784800) or # from run 6445
            # (unix_time>=1476739800 and unix_time<=1476741000) or # from run 8100 22 hour balloon launch
            # (unix_time>=1476999000 and unix_time<=1476999900) or # from run 8114 22 hour balloon launch but barely noticeable
            # (unix_time>=1477258200 and unix_time<=1477259100) or # from run 8129 22 hour balloon launch
            (unix_time>=1477511700 and unix_time<=1477512600) or # from run 8143 weird possible balloon launch
            (unix_time>=1477950300 and unix_time<=1477951500) or # from run 8168 22 hour balloon launch
            # (unix_time>=1478033400 and unix_time<=1478034000) or # from run 8173 22 hour balloon launch
            # (unix_time>=1478295300 and unix_time<=1478296200) or # from run 8188 22 hour balloon launch
            # (unix_time>=1478728500 and unix_time<=1478729400) or # from run 8213 22 hour balloon launch
            (unix_time>=1479231900 and unix_time<=1479232500) or # from run 8241

            # config 5
            (unix_time>=1449280500 and unix_time<=1449281100) or # from run 6513
            (unix_time>=1449610200 and unix_time<=1449612000) or # from run 6531
            (unix_time>=1450536000 and unix_time<=1450537200) or # from run 6584
            # (unix_time>=1450906200 and unix_time<=1450907100) or # from run 6606    22hr
            # (unix_time>=1451423700 and unix_time<=1451424600) or # from run 6635   22hr
            (unix_time>=1452008100 and unix_time<=1452009000) or # from run 6669
            # (unix_time>=1452115800 and unix_time<=1452116700) or # from run 6675    22hr
            (unix_time>=1452197700 and unix_time<=1452198600) or # from run 6679
            (unix_time>=1452213600 and unix_time<=1452214200) or # from run 6680
            (unix_time>=1452282000 and unix_time<=1452282600) or # from run 6684
            (unix_time>=1452298200 and unix_time<=1452298800) or # from run 6685    possible error
            (unix_time>=1452385500 and unix_time<=1452386400) or # from run 6690
            # (unix_time>=1452457800 and unix_time<=1452458700) or # from run 6694   22 hr
            (unix_time>=1452494100 and unix_time<=1452495000) or # from run 6696   possible error
            # (unix_time>=1452545100 and unix_time<=1452545880) or # from run 6700    could be error or 22hr
            # (unix_time>=1452636900 and unix_time<=1452637500) or # from run 6705   could be error or 22hr
            (unix_time>=1452715200 and unix_time<=1452716100) or # from run 6709   possible error
            (unix_time>=1452972300 and unix_time<=1452973440) or # from run 6724   possible error
            # (unix_time>=1453325400 and unix_time<=1453326600) or # from run 6743   22 hr
            (unix_time>=1453408500 and unix_time<=1453409400) or # from run 6747
            (unix_time>=1453930200 and unix_time<=1453931400) or # from run 6776
            # (unix_time>=1454535000 and unix_time<=1454536500) or # from run 6818   22 hr
            # (unix_time>=1455746400 and unix_time<=1455747900) or # from run 6889   22 hr
            (unix_time>=1456200900 and unix_time<=1456201800) or # from run 6916
            (unix_time>=1456392600 and unix_time<=1456393800) or # from run 6927
            (unix_time>=1456997400 and unix_time<=1456999200) or # from run 6962
            # (unix_time>=1457559000 and unix_time<=1457560800) or # from run 6994   22 hr
            (unix_time>=1460842800 and unix_time<=1460844600) or # from run 7119   22 hr // has CW contam cal pulsers
            # (unix_time>=1461620100 and unix_time<=1461621900) or # from run 7161   22 hr
            (unix_time>=1463002200 and unix_time<=1463004000) or # from run 7243  22 hr // has CW contam cal pulsers
            (unix_time>=1466501400 and unix_time<=1466503200) or # from run 7474
            (unix_time>=1466721900 and unix_time<=1466724600) or # from run 7486 22 hr // has CW contam cal pulsers
            (unix_time>=1466805600 and unix_time<=1466808300) or # from run 7489 22 hr // has CW contam cal pulsers
            (unix_time>=1466890200 and unix_time<=1466892000) or # from run 7494   22 hr // has CW contam cal pulsers
            (unix_time>=1467927600 and unix_time<=1467929700) or # from run 7552   22 hr
            # (unix_time>=1472333400 and unix_time<=1472335200) or # from run 7831   22 hr
            (unix_time>=1473111300 and unix_time<=1473112800) or # from run 7879    22 hr // has CW contam cal
            # (unix_time>=1473370500 and unix_time<=1473372900) or # from run 7899   22 hr
            # (unix_time>=1475011500 and unix_time<=1475013600) or # from run 7993   22 hr
            (unix_time>=1475185200 and unix_time<=1475187900) or # from run 8003 balloon 22hr // has CW contam cal pulsers
            # (unix_time>=1475358000 and unix_time<=1475359800) or # from run 8013 balloon 22h
            (unix_time>=1475529900 and unix_time<=1475531400) or # from run 8023 balloon 22hr // has CW contam cal pulsers
            # (unix_time>=1475702700 and unix_time<=1475704200) or # from run 8033 balloon 22hr
            (unix_time>=1476221400 and unix_time<=1476222300)): # from run 8069 balloon 22hr // has CW contam cal pulsers
            # (unix_time>=1476479700 and unix_time<=1476481800) # from run 8084 balloon 22hr

                bad_unix_time = True

        elif station_id == 3:

            # config 1 from undergrads
            if((unix_time>=1380234300 and unix_time<=1380235500) or # from run 1538, 22 hour balloon launch
            (unix_time>=1381008600 and unix_time<=1381010400) or # from run 1584, 22 hour balloon launch
            (unix_time>=1382476200 and unix_time<=1382477400) or # from run 1670, 22 hour balloon launch-ish
            (unix_time>=1382687400 and unix_time<=1382688600) or # from run 1682
            (unix_time>=1382712600 and unix_time<=1382713800) or # from run 1684, 15 hour spike
            (unix_time>=1382972700 and unix_time<=1382973300) or # from run 1698, 15 hour spike
            (unix_time>=1383688800 and unix_time<=1383691500) or # from run 1739, 22 hour balloon launch
            (unix_time>=1384060200 and unix_time<=1384060800) or # from run 1761
            (unix_time>=1384208700 and unix_time<=1384209900) or # from run 1770, 22 hour balloon launch
            (unix_time>=1384486200 and unix_time<=1384492800) or # from run 1786, repeated bursts over ~2 hrs
            (unix_time>=1389399600 and unix_time<=1389400800) or # from run 1980
            (unix_time>=1389744000 and unix_time<=1389747600) or # from run 2001, lots of activity, sweeps in phi
            (unix_time>=1390176600 and unix_time<=1390182000) or # from run 2025
            (unix_time>=1391027700 and unix_time<=1391028900) or # from run 2079, 22 hour balloon launch, but early?
            (unix_time>=1393652400 and unix_time<=1393660800) or # from run 2235, repeated bursts over ~2 hrs
            (unix_time>=1394846400 and unix_time<=1394856000) or # from run 2328, repeated bursts over ~2.5 hours
            (unix_time>=1395437400 and unix_time<=1395438600) or # from run 2363, 22 hour balloon launch
            (unix_time>=1397856300 and unix_time<=1397857800) or # from run 2526, 22 hour balloon launch

            # config 2
            (unix_time>=1390176600 and unix_time<=1390182000) or # from run 3533

            # config 3
            (unix_time>=1409954100 and unix_time<=1409956200) or # from run 3216, 22 hour balloon launch
            (unix_time>=1409986800 and unix_time<=1409988600) or # from run 3217
            (unix_time>=1412026200 and unix_time<=1412028000) or # from run 3332
            (unix_time>=1412284920 and unix_time<=1412287020) or # from run 3347, 22 hour balloon launch
            (unix_time>=1412544120 and unix_time<=1412546400) or # from run 3362, 22 hour balloon launch
            (unix_time>=1412803620 and unix_time<=1412805780) or # from run 3377, 22 hour balloon launch
            (unix_time>=1413897900 and unix_time<=1413899100) or # from run 3439
            (unix_time>=1413914400 and unix_time<=1413922200) or # from run 3440 big wide weird above ground
            (unix_time>=1414083600 and unix_time<=1414086300) or # from run 3449 , 2 spikes
            (unix_time>=1413550800 and unix_time<=1413552600) or # from run 3419, end of the run, before a software dominated run starts
            (unix_time>=1414674000 and unix_time<=1414675500) or # from run 3478
            (unix_time>=1415380500 and unix_time<=1415381400) or # from run 3520
            (unix_time>=1415460600 and unix_time<=1415461500) or # from run 3524
            (unix_time>=1415742000 and unix_time<=1415744100) or # from run 3540 22hr balloon
            (unix_time>=1416207300 and unix_time<=1416209700) or # from run 3568 2 small spikes
            (unix_time>=1416457800 and unix_time<=1416459000) or # from run 3579
            (unix_time>=1416909600 and unix_time<=1416910680) or # from run 3605
            (unix_time>=1416951000 and unix_time<=1416952500) or # from run 3608 22hr balloon
            (unix_time>=1417676400 and unix_time<=1417679400) or # from run 3647
            (unix_time>=1417742400 and unix_time<=1417743600) or # from run 3651
            (unix_time>=1417836600 and unix_time<=1417839300) or # from run 3656
            (unix_time>=1420317000 and unix_time<=1420318200) or # from run 3800
            (unix_time>=1420493700 and unix_time<=1420494600) or # from run 3810 22hr balloon
            (unix_time>=1420513200 and unix_time<=1420515000) or # from run 3811
            (unix_time>=1420598700 and unix_time<=1420600500) or # from run 3816
            (unix_time>=1420857900 and unix_time<=1420859700) or # from run 3830
            (unix_time>=1421019000 and unix_time<=1421020200) or # from run 3840 22hr balloon maybe?
            (unix_time>=1421101800 and unix_time<=1421103600) or # from run 3863 22hr balloon
            (unix_time>=1421723400 and unix_time<=1421723940) or # from run 3910
            (unix_time>=1421750700 and unix_time<=1421751720) or # from run 3912
            (unix_time>=1421868600 and unix_time<=1421881200) or # from run 3977 looks intentional
            (unix_time>=1421881200 and unix_time<=1421884680) or # from run 3978 continuation of thing above
            (unix_time>=1422048900 and unix_time<=1422049800) or # from run 3987 , 22 hour balloon launch
            (unix_time>=1422307200 and unix_time<=1422308100) or # from run 3995 22hr balloon
            (unix_time>=1423660800 and unix_time<=1423661700) or # from run 4132
            (unix_time>=1424819880 and unix_time<=1424820720) or # from run 4200
            (unix_time>=1428529500 and unix_time<=1428531000) or # from run 4412, 22 hour balloon launch
            (unix_time>=1429094400 and unix_time<=1429095600) or # from run 4445
            (unix_time>=1429615800 and unix_time<=1429617600) or # from run 4473
            (unix_time>=1429616700 and unix_time<=1429627500) or # from run 4474
            (unix_time>=1429733400 and unix_time<=1429734600) or # from run 4482
            (unix_time>=1431034500 and unix_time<=1431036900) or # from run 4557 , 22 hour balloon launch
            (unix_time>=1433365500 and unix_time<=1433367900) or # from run 4693
            (unix_time>=1435755600 and unix_time<=1435756500) or # from run 4829
            (unix_time>=1435791000 and unix_time<=1435791600) or # from run 4832
            (unix_time>=1436393700 and unix_time<=1436395500) or # from run 4867
            (unix_time>=1476740100 and unix_time<=1476741300) or # from run 7658
            (unix_time>=1477511400 and unix_time<=1477518300) or # from run 7704, big spike followed by nothing at all
            (unix_time>=1477604700 and unix_time<=1477605900) or # from run 7709,  22 hour balloon launch
            (unix_time>=1477950300 and unix_time<=1477951500) or # from run 7729
            (unix_time>=1479231600 and unix_time<=1479235800) or # from run 7802  , big spike followed by nothing at all

            # config 4
            (unix_time>=1448959200 and unix_time<=1448960100) or # from run 6009
            (unix_time>=1449610500 and unix_time<=1449611400) or # from run 6046 22 hour balloon launch
            (unix_time>=1450119900 and unix_time<=1450120500) or # from run 6077 possible 22 hour balloon launch
            (unix_time>=1450536360 and unix_time<=1450536720) or # from run 6098 spike is at end of time
            (unix_time>=1452116100 and unix_time<=1452116700) or # from run 6188 end of time and possible balloon launch
            (unix_time>=1452196800 and unix_time<=1452198600) or # from run 6193 could be balloon
            (unix_time>=1452213600 and unix_time<=1452214200) or # from run 6194
            (unix_time>=1452282300 and unix_time<=1452282900) or # from run 6198 could be balloon
            (unix_time>=1452298500 and unix_time<=1452299100) or # from run 6199 spike is at end of measured time
            (unix_time>=1452385800 and unix_time<=1452386400) or # from run 6203 spike is at end of measured time
            (unix_time>=1452457800 and unix_time<=1452458700) or # from run 6206 spike is at end of measured time, could be balloon
            (unix_time>=1452494100 and unix_time<=1452494700) or # from run 6208 spike is at end of measured time
            (unix_time>=1452544980 and unix_time<=1452545580) or # from run 6212 could be balloon
            (unix_time>=1452561120 and unix_time<=1452561480) or # from run 6213 spike is at end of measured time
            (unix_time>=1452637020 and unix_time<=1452637260) or # from run 6219 spike is at end of measured time, could be balloon
            (unix_time>=1452715320 and unix_time<=1452715680) or # from run 6223 spike is at end of measured time
            (unix_time>=1452972660 and unix_time<=1452973020) or # from run 6239 spike is at end of measured time
            (unix_time>=1453325400 and unix_time<=1453326300) or # from run 6259 could be balloon
            (unix_time>=1453930500 and unix_time<=1453931100) or # from run 6295 could be balloon
            (unix_time>=1454535000 and unix_time<=1454536200) or # from run 6328 could be balloon
            (unix_time>=1454911200 and unix_time<=1454911800) or # from run 6349 spike is at end of measured time could match below
            (unix_time>=1454911200 and unix_time<=1454912100) or # from run 6350 spike is at start of measured time could match above
            (unix_time>=1455746400 and unix_time<=1455747300) or # from run 6397 could be balloon
            (unix_time>=1456374300 and unix_time<=1456374900) or # from run 6433
            (unix_time>=1457559300 and unix_time<=1457560500) or # from run 6501 could be balloon
            (unix_time>=1460843100 and unix_time<=1460844600) or # from run 6618 spike is at start of measured time, could be balloon
            (unix_time>=1467927840 and unix_time<=1467929640) or # from run 7052 could be balloon
            (unix_time>=1473371280 and unix_time<=1473372180) or # from run 7458 could be balloon
            (unix_time>=1475186100 and unix_time<=1475187000) or # from run 7562 could be balloon
            (unix_time>=1475530500 and unix_time<=1475531700) or # from run 7584 could be balloon
            (unix_time>=1476221400 and unix_time<=1476222600)): # from run 7625 could be balloon

                bad_unix_time = True

        elif station_id == 5:
            pass

        return bad_unix_time

        