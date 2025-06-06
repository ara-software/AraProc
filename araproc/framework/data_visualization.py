import numpy as np
import matplotlib.pyplot as plt
from araproc.framework import waveform_utilities as wu
from araproc.framework import map_utilities as mu
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True) # tell root not to open gui windows

"""
It's important to use a non-GUI backend!
Otherwise using matplotlib functions in a loop can memory leak.
See this thread for an interesting discussion:
https://github.com/matplotlib/matplotlib/issues/20300.

But also, "agg" sometimes give segfaults, so instead we use "pdf".
"""
import matplotlib # noqa : E402
matplotlib.use("pdf")

def plot_waveform_bundle(
    waveform_dict,
    time_or_freq = "time",
    output_file_path = None
    ):
    """
    Code to plot a bundle of waveforms.
    Please be aware that because this is matplotlib,
    this function is very slow (about 2.5 s per event).
    As a result, you probably shouldn't use it to draw many events.
    It's meant for debugging and visualization.
    """
    
    ####################
    # sanitize inputs
    ####################

    # make sure they have requested something valid (time vs freq domain)
    time_or_freq_options = ["time", "freq"]
    if time_or_freq not in time_or_freq_options:
        raise KeyError(f"Requested option ({time_or_freq}) is not in approved list: {time_or_freq_options}")

    xlabel_options = {
        "time" : "Time (ns)",
        "freq" : "Frequency (GHz)"
    }
    ylabel_options = {
        "time" : "Voltage (mV)",
        "freq" : "Spectrum (dBm)"
    }    

    if not isinstance(output_file_path, str):
        raise TypeError("Path to output file must be a string")
        
    ####################
    # actually make plots
    ####################

    # get the waves to plot
    tgraphs_to_plot = waveform_dict

    # set up fig and axes
    fig, axd = plt.subplot_mosaic(
                            [["ch0", "ch1", "ch2", "ch3"],
                             ["ch4", "ch5", "ch6", "ch7"],
                             ["ch8", "ch9", "ch10", "ch11"],
                             ["ch12", "ch13", "ch14", "ch15"]],
                              figsize=(10, 7), 
                              sharex=True,
                              sharey=True,
                              layout="constrained")

    # draw the graphs, label each one appropriately
    ymin = 1e100
    ymax = -1e100
    spec2Tot = {}
    for wave_key in tgraphs_to_plot.keys():
        times, volts = wu.tgraph_to_arrays(tgraphs_to_plot[wave_key]) # 'times' in 'ns' and 'volts' in 'mV'
        xvals = times
        yvals = volts
     
        if time_or_freq == "freq":
            # if they frequested frequency domain, do the FFT
            freqs, spectrum = wu.time2freq(times, volts) # 'freqs' in 'GHz' and 'spectrum' (complex) in 'mV'
            xvals = freqs
            
            # Convert yvals to dBm
            # np.abs(spectrum)**2 makes spectrum in mV^2. For power, you do P = V^2/R , here R = Z_0 = 50 Ohm.
            # mV**2/50 =  mW * 1e-3. Power is always reperesented as dBm in dB scale which is 10*log10(P in mW)
            yvals = 10*np.log10(np.abs(spectrum)**2 / 50 / 1e3) # from mV to dBm
            
            # save total of square of absolute spectrum for reference
            spec2Tot[f"ch{wave_key}"] = (np.abs(spectrum)**2).sum()

        ymin = min(ymin, yvals.min())
        ymax = max(ymax, yvals.max())

        if wave_key < 8: # VPol channels
            axd[f"ch{wave_key}"].plot(xvals, yvals, color = 'purple', lw = 1) # make the plot
            axd[f"ch{wave_key}"].set_title(f"Channel {wave_key}")
        else: # HPol channels
            axd[f"ch{wave_key}"].plot(xvals, yvals, color = 'green', lw = 1) # make the plot
            axd[f"ch{wave_key}"].set_title(f"Channel {wave_key}")

    # label axes  
    for ax in [axd["ch12"], axd["ch13"], axd["ch14"], axd["ch15"]]:
        ax.set_xlabel(xlabel_options[time_or_freq])
    for ax in [axd["ch0"], axd["ch4"], axd["ch8"], axd["ch12"]]:
        ax.set_ylabel(ylabel_options[time_or_freq])
    
    # make limits look nice
    for ch in axd:
      if time_or_freq == "freq":
          axd[ch].set_xlim(0, 1)
          # limit y range downwards
          ymin = max(ymin, -25)
          axd[ch].set_ylim([ymin-5,ymax+5])

          # add secondary axis to help with CW filtering
          ax2 = axd[ch].twinx()
          y1min, y1max = axd[ch].get_ylim()
          ## convert back to just the spec^2
          spec2min = 50.*1e3*10.**(y1min/10.)
          spec2max = 50.*1e3*10.**(y1max/10.)
          ## normalize by total power in spectrum
          fmin = np.log10(spec2min/spec2Tot[ch])
          fmax = np.log10(spec2max/spec2Tot[ch])
          ax2.set_ylim(fmin, fmax)
          ax2.set_yticks(np.arange(np.ceil(fmin), np.floor(fmax)+0.1, 1)) # ticks MUST fall within the [fmin, fmax] range or y-axes will not correspond correctly
          ax2.axhline(y=0.0, c='lightgray', linestyle='--')
          if ch in ["ch3", "ch7", "ch11", "ch15"]:
            ax2.set_ylabel("lg(Fractional Power)")
      else:
          axd[ch].set_ylim(ymin, ymax)

    # save figure
    fig.savefig(output_file_path)

    # careful cleanup
    plt.close(fig)
    del fig, axd

def plot_skymap(the_map,
                plane_wave_elevation = None, 
                station_id = None, landmarks = None,
                calpulser_indices = None, spice_depth = None,
                output_file_path = None
                ):
    """
    This function returns skymap for given map type
    Parameters
    ----------
    the_map:
        reconstruction results for given map
    plane_wave_elevation : float
        elevation angle (in degrees) from plane wave reco 
    calpulser_indices : list ## example [0,1,2,3]
        which calpulsers you want to see in your skymap
    landmarks: list ## example ['ICL',IC22S','SPT','IC1S','SPIce','WT']
        which landmarks you want to see in your skymap
    spice_depth : int/float
        the depth of spice pulser ## example -1451.3
    output_file_path: str
        path to the output file
  
    Returns
    -------
    A plotted skymap in pdf
    """
    
    if the_map is None:
        raise Exception("the_map is None")

    if not isinstance(output_file_path, str):
        raise TypeError("Path to output file must be a string")
   
    corr_peak, peak_phi, peak_theta = mu.get_corr_map_peak(the_map)
    the_map.SetTitle(f"Peak #phi/#theta/Corr = {peak_phi:.1f}^{{o}}/ {peak_theta:.1f}^{{o}}/ {corr_peak:.2f}")
    the_map.GetXaxis().SetTitle("Azimuth, #phi (^{o})")
    the_map.GetYaxis().SetTitle("Elevation, #theta (^{o})")
    the_map.GetZaxis().SetTitle("Correlation")
    ROOT.gStyle.SetTitleFontSize(0.04)
 
    the_map.GetZaxis().SetRangeUser(0, corr_peak)
    the_map.GetXaxis().CenterTitle(1)
    the_map.GetYaxis().CenterTitle(1)
    the_map.GetZaxis().CenterTitle(1)

    c = ROOT.TCanvas("c", "c", 700, 500)
    c.cd()
    the_map.Draw("COLZ") # keeping this off for now: the_map.Draw("z aitoff")

    markers = []  # Keep references to markers to avoid garbage collection
    labels = []   # Keep references to labels
    ## Add plan wave elevation if requested
    if plane_wave_elevation is not None:
       horizontal_line0 = ROOT.TLine(-180, plane_wave_elevation, 180, plane_wave_elevation)  # Draw line from phi=-180 to phi=180
       horizontal_line0.SetLineColor(ROOT.kOrange)
       horizontal_line0.SetLineStyle(2)  # Dashed line
       horizontal_line0.SetLineWidth(2)
       horizontal_line0.Draw("SAME")
       label0 = ROOT.TLatex(110, plane_wave_elevation + 5, "plane wave")  # Offset for clarity
       label0.SetTextColor(ROOT.kOrange)
       label0.SetTextSize(0.03)
       label0.Draw("SAME")
       labels.append(label0)
      
    ## Add known locations to the skymap 
    landmark_dict = mu.AraGeom(station_id).get_known_landmarks(landmarks, calpulser_indices, spice_depth)

    for entry in landmark_dict.keys():
        if entry == 'critical_angle':
           critical_ang = landmark_dict[entry]
           horizontal_line = ROOT.TLine(-180, critical_ang, 180, critical_ang)  # Draw line from phi=-180 to phi=180 
           horizontal_line.SetLineColor(ROOT.kRed)
           horizontal_line.SetLineStyle(2)  # Dashed line
           horizontal_line.SetLineWidth(2)
           horizontal_line.Draw("SAME")
           label1 = ROOT.TLatex(150,  critical_ang + 5, "#theta_{c}")  # Offset for clarity
           label1.SetTextColor(ROOT.kRed)
           label1.SetTextSize(0.03)
           label1.Draw("SAME")
           labels.append(label1)
           continue

        phi = landmark_dict[entry][2]
        theta = landmark_dict[entry][1]
        # Draw the marker
        marker = ROOT.TMarker(phi, theta, 29)  # Style 29: Star
        color = ROOT.kBlack if "CP" in entry else ROOT.kRed
        marker.SetMarkerColor(color)

        marker.SetMarkerSize(2.0)
        marker.Draw("SAME")
        markers.append(marker)

        # Draw the label
        if entry in ['IC1S','SPT','CP1','CP3']:
           offset = 5 
        else:
           offset = -5
        label = ROOT.TLatex(phi + offset, theta - offset, entry)  # Offset for clarity
        label.SetTextColor(ROOT.kMagenta)
        label.SetTextSize(0.02)
        label.Draw("SAME")
        labels.append(label)

        if entry == "ICL":
           vertical_line = ROOT.TLine(phi, -90, phi, 90)  # Draw line from theta=-90 to theta=90
           vertical_line.SetLineColor(ROOT.kCyan)
           vertical_line.SetLineStyle(2)  # Dashed line
           vertical_line.SetLineWidth(2)
           vertical_line.Draw("SAME")
           label = ROOT.TLatex(phi + 2, theta + 30, "#phi_{SP}")  # Offset for clarity
           label.SetTextColor(ROOT.kCyan)
           label.SetTextSize(0.03)
           label.Draw("SAME")
           labels.append(label)

    ROOT.gStyle.SetPalette(112) # viridis
    ROOT.gPad.SetRightMargin(0.15) # make space for the z axis
    c.SaveAs(output_file_path)
    ROOT.gPad.SetRightMargin(0) # reset, so we don't affect settings globally
    c.Close()

    del c
