import numpy as np
import matplotlib.pyplot as plt
from araproc.framework import constants as const
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

def get_wf_color(channel_number, excluded_channels): 
    """
    Given a channel number and a list of excluded channels, return a purple 
      color for all VPols (channel_number<8) in purple, a green color for
      HPols (channel_number > 8), and a slightly greyer purple or green if the 
      channel is masked.
    """
    if channel_number in const.vpol_channel_ids: # Vpols
        if channel_number in excluded_channels: # Channel is masked
            return 'thistle'
        else:
            return 'purple'
    else: # HPols
        if channel_number in excluded_channels: # channel is masked
            return '#c8d9bf' # A lighter version of "xkcd:greenish grey"
        else: 
            return 'green'
        
def plot_waveform(
    waveform_dict,
    channel,
    time_or_freq = "time",
    output_file_path = None, 
    excluded_channels = (),
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
    fig, ax = plt.subplots()

    # draw the graphs, label each one appropriately
    ymin = 1e100
    ymax = -1e100
    spec2Tot = {}
    wave_key = channel
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

    # make the plot
    ax.plot(xvals, yvals, color=get_wf_color(wave_key, excluded_channels), lw = 1)
    ax.set_title(f"Channel {wave_key}")

    # Add warning to channels that are masked
    if wave_key in excluded_channels: 
        ax.annotate("Masked to Analysis", (0.03,0.87), xycoords='axes fraction', c='grey')

    # label axes  
    ax.set_xlabel(xlabel_options[time_or_freq])
    ax.set_ylabel(ylabel_options[time_or_freq])
    
    # make limits look nice
    if time_or_freq == "freq":
        ax.set_xlim(0, 1)
        # limit y range downwards
        ymin = max(ymin, -25)
        ax.set_ylim([ymin-5,ymax+5])

        # add secondary axis to help with CW filtering
        ax2 = ax.twinx()
        y1min, y1max = ax.get_ylim()
        ## convert back to just the spec^2
        spec2min = 50.*1e3*10.**(y1min/10.)
        spec2max = 50.*1e3*10.**(y1max/10.)
        ## normalize by total power in spectrum
        fmin = np.log10(spec2min/spec2Tot[f"ch{wave_key}"])
        fmax = np.log10(spec2max/spec2Tot[f"ch{wave_key}"])
        ax2.set_ylim(fmin, fmax)
        ax2.set_yticks(np.arange(np.ceil(fmin), np.floor(fmax)+0.1, 1)) # ticks MUST fall within the [fmin, fmax] range or y-axes will not correspond correctly
        ax2.axhline(y=0.0, c='lightgray', linestyle='--')
        ax2.set_ylabel("lg(Fractional Power)")
    else:
        ax.set_ylim(ymin, ymax)
    plt.tight_layout()

    # save figure
    fig.savefig(output_file_path)

    # careful cleanup
    plt.close(fig)
    del fig, ax


def plot_waveform_bundle(
    waveform_dict,
    time_or_freq = "time",
    output_file_path = None, 
    excluded_channels = (),
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

        # make the plot
        axd[f"ch{wave_key}"].plot(xvals, yvals, color=get_wf_color(wave_key, excluded_channels), lw = 1)
        axd[f"ch{wave_key}"].set_title(f"Channel {wave_key}")

        # Add warning to channels that are masked
        if wave_key in excluded_channels: 
            axd[f"ch{wave_key}"].annotate("Masked to Analysis", (0.03,0.87), xycoords='axes fraction', c='grey')

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
                plane_wave_elevation, 
                station_id, 
                map_type, 
                landmarks = None,
                landmark_labels = None,
                calpulser_indices = None,
                cp_labels = None,
                list_of_channels = None,
                spice_depth = None,
                aravertex_results = None,
                output_file_path = None,
                include_legend_helper=True,
                write_peak_on_map=True,
                include_critical_angle_rt=True,
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
    cp_labels : list of str, optional
        Custom display text for the calpulsers requested in calpulser_indices,
        matched 1-to-1 in order (assumes landmark_dict keys them as
        "CP" + index). Only works when calpulser_indices is listed
        individually, not "all".
    landmarks: list ## example ['ICL',IC22S','SPT','IC1S','SPIce','WT', 'SP']
        which landmarks you want to see in your skymap
        Selecting "SP" will get you only a vertical dashed line
    landmark_labels : list of str, optional
        Custom display text for the markers requested in landmarks, matched
        1-to-1 in order. Only works when landmarks is listed individually,
        not "all".
    list_of_channels : string ## example 'vpols', 'hpols' or 'all'
        which channels you want to see in your skymap
    spice_depth : int/float
        the depth of spice pulser ## example -1451.3
    output_file_path: str
        path to the output file
    aravertex_results : dict or None
        AraVertex reconstruction results; no marker is shown if hit-channel SNR < 5   
    include_legend_helper : bool
        Whether or not to write the little reminder of what star vs x means for reco
    write_peak_on_map : bool
        Whether you want the peak theta/phi/corr written on the map
        
    Returns
    -------
    A plotted skymap in pdf
    """
    
    if the_map is None:
        raise Exception("the_map is None")

    if not isinstance(output_file_path, str):
        raise TypeError("Path to output file must be a string")

    if landmark_labels is not None:
        custom_labels = _zip_custom_labels(
            landmarks, landmark_labels, "landmark_labels", "landmarks"
        )
    else:
        custom_labels = {}
    if cp_labels is not None:
        custom_labels.update(_zip_custom_labels(
            calpulser_indices, cp_labels, "cp_labels", "calpulser_indices",
            key_fn=lambda idx: f"CP{idx}",
        ))

    corr_peak, peak_phi, peak_theta = mu.get_corr_map_peak(the_map)
    if write_peak_on_map:
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

    # AraVertex marker (not shown if hit-channel SNR < 5)
    if aravertex_results is not None:
        for r in aravertex_results.values():
            if r["valid"]:
                phi = r["phi"]
                theta = r["theta"]
                break
        else:
            phi = theta = None

        # skip AraVertex default failure return (no valid hit pairs, e.g. SNR < 5)
        # when this happens, AraVetex saves placeholder reconstruction (theta = 90 degree and phi = 0 degree)
        if phi is not None and not (abs(theta - 90.0) < 1e-3 and phi == 0.0):
            av_marker = ROOT.TMarker(phi, theta, 34)
            av_marker.SetMarkerColor(ROOT.kBlue)
            av_marker.SetMarkerSize(1.0)
            av_marker.Draw("SAME")
            markers.append(av_marker)

            av_label = ROOT.TLatex(phi + 3, theta + 3, "AraVertex")
            av_label.SetTextColor(ROOT.kBlue)
            av_label.SetTextSize(0.025)
            av_label.Draw("SAME")
            labels.append(av_label)
 
    # Reconstruction sphere radius 
    if map_type in ["pulser_v", "pulser_h"]:
        radius_map = float(const.calpulser_r_library[station_id])
    else:
        radius_map = float(const.distant_events_r_library[station_id])
    
    solution = 1 if map_type in {"distant_v_ref", "distant_h_ref"} else 0 

    # Add known locations to the skymap 
    
    landmark_dict = mu.AraGeom(station_id).get_known_landmarks(landmarks, radius_map, calpulser_indices,list_of_channels, spice_depth, solution=solution)
    if include_critical_angle_rt == False:
        del landmark_dict["critical_angle_rt"]

    marker_status = landmark_dict.get("_marker_status", {})

    for entry in landmark_dict.keys():
    
        if entry == "_marker_status":
            continue
        if entry == 'critical_angle_rt' and include_critical_angle_rt:
            critical_ang_rt = landmark_dict[entry]
            horizontal_line_rt = ROOT.TLine(-180, critical_ang_rt, 180, critical_ang_rt) # Draw line from phi=-180 to phi=180 
            horizontal_line_rt.SetLineColor(ROOT.kRed)
            horizontal_line_rt.SetLineStyle(2) # Dashed line
            horizontal_line_rt.SetLineWidth(2) 
            horizontal_line_rt.Draw("SAME")
            label_rt = ROOT.TLatex(110, critical_ang_rt + 5, "#theta_{c}")  # Offset for clarity
            label_rt.SetTextColor(ROOT.kRed)
            label_rt.SetTextSize(0.03)
            label_rt.Draw("SAME")
            labels.append(label_rt)
            continue

        phi = landmark_dict[entry][2]
        theta = landmark_dict[entry][1]

        if entry=="SP":
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
           continue
             
        # Draw the marker   
        status = marker_status.get(entry, "raytraced")
        marker_style = 5 if status == "sl_fallback" else 29
        marker = ROOT.TMarker(phi, theta, marker_style)

        color = ROOT.kBlue if "CH" in entry else ROOT.kBlack if "CP" in entry else ROOT.kRed
        marker.SetMarkerColor(color)

        # Set Markersize for SL fallback marker and Ray-Traced Marker
        if status == "sl_fallback":
            marker.SetMarkerSize(1.6)
        else:
            marker.SetMarkerSize(2.0)

        marker.Draw("SAME")
        markers.append(marker)

        # Draw the label
        if entry in ['IC1S','SPT','CP1','CP3']:
           offset = 5 
        else:
           offset = -5
        label_text = custom_labels.get(entry, entry)
        label = ROOT.TLatex(phi + offset, theta - offset, label_text)  # Offset for clarity
        label.SetTextColor(ROOT.kMagenta)
        label.SetTextSize(0.02)
        label.Draw("SAME")
        labels.append(label)

    ROOT.gStyle.SetPalette(112) # viridis
    ROOT.gPad.SetRightMargin(0.15) # make space for the z axis

    if include_legend_helper:
        # Add note about markers
        caption_text = (
        "#lower[0.35]{#scale[2.0]{*}} : ray-traced direction;  "
        "#scale[1.0]{#times} : straight-line fallback when no ray-traced solution is found"
        )
        caption = ROOT.TLatex(0.5, 0.92, caption_text)
        caption.SetNDC(True)
        caption.SetTextAlign(22)
        caption.SetTextSize(0.030)
        caption.SetTextColor(ROOT.kRed)
        caption.Draw()
        labels.append(caption)

    c.SaveAs(output_file_path)
    ROOT.gPad.SetRightMargin(0) # reset, so we don't affect settings globally
    c.Close()

    del c


### starty

def _th2d_to_grid(hist):
    """
    Pull a ROOT TH2D's bin edges and contents out into numpy arrays.

    Returns
    -------
    phi_edges : (nx+1,) array, degrees
    theta_edges : (ny+1,) array, degrees
    z : (ny, nx) array of bin contents, matching pcolormesh's (row, col)
        convention
    """
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()

    phi_edges = np.array([xaxis.GetBinLowEdge(i) for i in range(1, nx + 2)])
    theta_edges = np.array([yaxis.GetBinLowEdge(j) for j in range(1, ny + 2)])

    z = np.empty((ny, nx))
    for iy in range(1, ny + 1):
        for ix in range(1, nx + 1):
            z[iy - 1, ix - 1] = hist.GetBinContent(ix, iy)

    return phi_edges, theta_edges, z


def _deg2rad(phi_deg, theta_deg):
    """Degrees -> radians, clipped to what a mollweide axes will accept."""
    phi_deg = np.asarray(phi_deg, dtype=float)
    theta_deg = np.asarray(theta_deg, dtype=float)
    phi = np.radians(np.clip(phi_deg, -180, 180))
    theta = np.radians(np.clip(theta_deg, -90, 90))
    return phi, theta


def _draw_constant_theta(ax, theta_val, n_samples=181, **line_kwargs):
    """
    Draw a curve of constant elevation. In Mollweide this comes out as a
    straight horizontal line (parallels don't curve), but we still sample
    many points across the full width for a clean, fully-clipped line.
    """
    phis = np.linspace(-180, 180, n_samples)
    thetas = np.full_like(phis, theta_val)
    x, y = _deg2rad(phis, thetas)
    ax.plot(x, y, **line_kwargs)


def _draw_constant_phi(ax, phi_val, n_samples=181, **line_kwargs):
    """
    Draw a curve of constant azimuth. In Mollweide this is a curved
    (elliptical) meridian arc, so it must be sampled point-by-point across
    theta and projected -- it is NOT a straight vertical line except at
    phi = 0.
    """
    thetas = np.linspace(-90, 90, n_samples)
    phis = np.full_like(thetas, phi_val)
    x, y = _deg2rad(phis, thetas)
    ax.plot(x, y, **line_kwargs)


def _zip_custom_labels(entries, custom, param_name, entries_param_name, key_fn=lambda x: x):
    """
    Build an {entry_in_landmark_dict: custom_label} lookup, matching
    `custom` 1-to-1 against `entries` in order. Shared by landmark_labels
    (zipped directly against landmarks) and cp_labels (zipped against
    calpulser_indices, with key_fn turning e.g. "1" into "CP1" to match
    how those entries are actually keyed in landmark_dict).
    """
    if custom is None:
        return {}
    if entries is None or len(custom) != len(entries):
        raise ValueError(
            f"{param_name} must be the same length as {entries_param_name}, "
            f"matched 1-to-1 in order"
        )
    return {key_fn(e): c for e, c in zip(entries, custom)}


def plot_skymap_mpl(
    the_map,
    plane_wave_elevation,
    station_id,
    map_type,
    landmarks=None,
    landmark_labels=None,
    calpulser_indices=None,
    cp_labels=None,
    list_of_channels=None,
    spice_depth=None,
    aravertex_results=None,
    output_file_path=None,
    include_legend_helper=True,
    write_peak_on_map=True,
    include_critical_angle_rt=True,
    flip_longitude=False,     # True for an RA-like convention (phi increasing to the left)
    cmap="viridis",
    dpi=200,                  # resolution of the rasterized correlation map (see below)
):
    """
    Drop-in matplotlib (w/ Mollweide projection) replacement for plot_skymap(). 
    This has the same signature as plot_skymap, but uses pcolormesh
    and matplotlib to do the plotting, rather than ROOT and TH2Ds.

    Parameters
    ----------
    flip_longitude : bool
        If True, mirrors phi (equivalent to viewing the sky from behind /
        an RA-style convention). Off by default, which matches the
        original TH2D's left-to-right phi orientation.
    dpi : int
        Resolution used for the rasterized correlation map (see note
        below). Higher = crisper, but you probably don't need more than 200.
    landmark_labels : list of str, optional
        Custom display text for the markers requested in `landmarks`,
        matched 1-to-1 in order. Doesn't work if `landmarks` is the "all"
        sentinel rather than an explicit list.
    cp_labels : list of str, optional
        Same idea, but for calpulser_indices. Assumes landmark_dict keys
        calpulser entries as "CP" + index (e.g. index "1" -> "CP1") --
        worth double-checking against what get_known_landmarks actually
        returns if labels don't land where you expect.
    """
    if the_map is None:
        raise Exception("the_map is None")
    if not isinstance(output_file_path, str):
        raise TypeError("Path to output file must be a string")
    if landmark_labels is not None:
        custom_labels = _zip_custom_labels(
            landmarks, landmark_labels, "landmark_labels", "landmarks"
        )
    else:
        custom_labels = {}
    if cp_labels is not None:
        custom_labels.update(_zip_custom_labels(
            calpulser_indices, cp_labels, "cp_labels", "calpulser_indices",
            key_fn=lambda idx: f"CP{idx}",
        ))

    # Create figure
    import matplotlib as mpl
    mpl.rcParams.update({'font.size':28})
    mpl.rcParams['font.family'] = 'STIXGeneral'
    mpl.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams.update({
        "axes.titlesize": 28,
        "axes.labelsize": 28,
        "xtick.labelsize": 28,
        "ytick.labelsize": 28,
        "legend.fontsize": 28
    })

    def flip(phi_deg):
        return -phi_deg if flip_longitude else phi_deg

    corr_peak, peak_phi, peak_theta = mu.get_corr_map_peak(the_map)
    phi_edges, theta_edges, z = _th2d_to_grid(the_map)

    if flip_longitude:
        # negating reverses the monotonic order of the edges, so the edges
        # and the corresponding columns of z both need to be flipped back
        # into increasing order together
        phi_edges = -phi_edges[::-1]
        z = z[:, ::-1]

    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_axes([0.05, 0.22, 0.90, 0.72], projection="mollweide")

    # --- the correlation map itself ---
    px, py = _deg2rad(phi_edges, theta_edges)
    mesh = ax.pcolormesh(
        px, py, z, cmap=cmap, vmin=0, vmax=corr_peak, shading="auto",
        rasterized=True,
    )
    fig.canvas.draw()
    pos = ax.get_position()

    ax.set_xlabel(r"Azimuth, $\phi$ (°)", labelpad=25)
    ax.set_ylabel(r"Elevation, $\theta$ (°)")
    fig.canvas.draw()   # re-draw so the xlabel's real rendered position is available

    xlabel_bbox = ax.xaxis.label.get_window_extent().transformed(fig.transFigure.inverted())
    cbar_gap = 0.02 # gap between the xlabel and the colorbar, in figure-fraction
    cbar_height = 0.03
    cax = fig.add_axes([pos.x0, xlabel_bbox.y0 - cbar_gap - cbar_height, pos.width, cbar_height])
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cbar.set_label("Correlation")

    if write_peak_on_map:
        ax.set_title(
            f"Peak \u03c6/\u03b8/Corr = {peak_phi:.1f}\u00b0/ {peak_theta:.1f}\u00b0/ {corr_peak:.2f}",
            fontsize=11,
        )
    ax.grid(True, alpha=0.5)
    ax.tick_params(axis='x', colors='white')   # or 'gray' / '0.7' / etc.

    # --- plane wave elevation ---
    if plane_wave_elevation is not None:
        _draw_constant_theta(
            ax, plane_wave_elevation,
            color="orange", linestyle="--", linewidth=2,
        )
        lx, ly = _deg2rad(flip(110), plane_wave_elevation - 10)
        ax.text(lx, ly, "plane wave", color="orange", fontsize=9)

    # --- AraVertex marker (skip the SNR<5 placeholder return) ---
    if aravertex_results is not None:
        phi = theta = None
        for r in aravertex_results.values():
            if r["valid"]:
                phi, theta = r["phi"], r["theta"]
                break
        if phi is not None and not (abs(theta - 90.0) < 1e-3 and phi == 0.0):
            mx, my = _deg2rad(flip(phi), theta)
            ax.plot(mx, my, marker="D", color="blue", markersize=7, linestyle="none")
            lx, ly = _deg2rad(flip(phi + 3), theta + 3)
            ax.text(lx, ly, "AraVertex", color="blue", fontsize=8)

    # --- reconstruction sphere radius / solution, unchanged from plot_skymap ---
    if map_type in ["pulser_v", "pulser_h"]:
        radius_map = float(const.calpulser_r_library[station_id])
    else:
        radius_map = float(const.distant_events_r_library[station_id])
    solution = 1 if map_type in {"distant_v_ref", "distant_h_ref"} else 0

    landmark_dict = mu.AraGeom(station_id).get_known_landmarks(
        landmarks, radius_map, calpulser_indices, list_of_channels, spice_depth,
        solution=solution,
    )
    if include_critical_angle_rt == False:
        del landmark_dict["critical_angle_rt"]
    marker_status = landmark_dict.get("_marker_status", {})

    for entry, val in landmark_dict.items():
        if entry == "_marker_status":
            continue

        if entry == 'critical_angle_rt':
            critical_ang_rt = val
            _draw_constant_theta(
                ax, critical_ang_rt,
                color="red", linestyle="--", linewidth=2,
            )
            lx, ly = _deg2rad(flip(110), critical_ang_rt + 5)
            ax.text(lx, ly, r"$\theta_c$", color="red", fontsize=9)
            continue
    
        phi = val[2]
        theta = val[1]

        if entry == "SP":
            _draw_constant_phi(
                ax, flip(phi),
                color="C9", linestyle="--", linewidth=2,
            )
            lx, ly = _deg2rad(flip(phi + 20), 50) # fix the theta value to 50 (total guess)
            ax.text(lx, ly, r"$SP$", color="C9", fontsize=28)
            continue

        status = marker_status.get(entry, "raytraced")
        marker = "x" if status == "sl_fallback" else "*"
        markersize = 9 if status == "sl_fallback" else 13
        color = "blue" if "CH" in entry else "C4" if "CP" in entry else "red"

        mx, my = _deg2rad(flip(phi), theta)
        import matplotlib.patheffects as pe
        ax.plot(mx, my, marker=marker, color="C6", markersize=markersize,
                linestyle="none", alpha=0.75)
        #markeredgewidth=1.5, alpha=0.75,
        #        markeredgecolor="white")


        offset = 5 if entry in ["IC1S", "SPT", "CP1", "CP3"] else -5
        lx, ly = _deg2rad(flip(phi + offset), theta - offset)
        label_text = custom_labels.get(entry, entry)
        txt = ax.text(lx, ly, label_text, color="C6", fontsize=28)
        # txt.set_path_effects([pe.withStroke(linewidth=1.5, foreground="white")])
    
    if include_legend_helper:
        caption = (
            "\u2605 : ray-traced direction     "
            "\u2715 : straight-line fallback when no ray-traced solution is found"
        )
        fig.text(0.5, 0.93, caption, ha="center", fontsize=8, color="red")
    
    ax.set_longitude_grid(60)
    ax.set_latitude_grid(30)

    fig.savefig(output_file_path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
