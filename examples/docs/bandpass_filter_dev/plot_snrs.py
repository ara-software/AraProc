import matplotlib.pyplot as plt
import numpy as np

station = 2

file = np.load(f"values_{station}.npz")

if station==2:
    bins_1d_snr = np.linspace(20,40,20)
    bins_2d_snr = np.linspace(20,40,40)
    bins_1d_corr = np.linspace(0.7,1.1,20)
    bins_2d_corr = np.linspace(0.7,1.1,40)
else:
    bins = np.linspace(10,20,20)

fig, (ax, ax2) = plt.subplots(1,2, figsize=(10,5))
ax.hist(file["dd_vals"], bins=bins_1d_snr, histtype="step", label=f"Dedis + CW: \nmed = {np.median(file['dd_vals']):.1f}", lw=3)
ax.hist(file["filt_vals"], bins=bins_1d_snr, histtype="step", label=f"Dedis + CW + Band:\n med = {np.median(file['filt_vals']):.1f}", lw=3)
ax.set_ylabel("Number of Events")
ax.set_xlabel("Vpol Average SNR")
ax.legend(loc="upper left")
ax2.hist2d(file["dd_vals"], file["filt_vals"], bins=[bins_2d_snr,bins_2d_snr], cmin=1)
ax2.set_xlabel("Vpol Average SNR (Dedis + CW )")
ax2.set_ylabel("Vpol Average SNR (Dedis + CW + Bandpass)")
ax2.set_aspect("equal")
ax2.plot(bins_2d_snr, bins_2d_snr, '--', color='C1')
fig.tight_layout()
fig.savefig(f"A{station}_snr.png", dpi=300)
del fig, ax, ax2


fig, (ax, ax2) = plt.subplots(1,2, figsize=(10,5))
ax.hist(file["dd_corr_vals"], bins=bins_1d_corr, histtype="step", label=f"Dedis + CW:\n med = {np.median(file['dd_corr_vals']):.2f}", lw=3)
ax.hist(file["filt_corr_vals"], bins=bins_1d_corr, histtype="step", label=f"Dedis + CW + Band:\n med = {np.median(file['filt_corr_vals']):.2f}", lw=3)
ax.set_ylabel("Number of Events")
ax.set_xlabel("Vpol Map Correlation")
ax.legend(loc="upper left")
ax2.hist2d(file["dd_corr_vals"], file["filt_corr_vals"], bins=[bins_2d_corr,bins_2d_corr], cmin=1)
ax2.set_xlabel("Vpol Map Corr (Dedis + CW )")
ax2.set_ylabel("Vpol Map Coor(Dedis + CW + Bandpass)")
ax2.set_aspect("equal")
ax2.plot(bins_2d_corr, bins_2d_corr, '--', color='C1')
fig.tight_layout()
fig.savefig(f"A{station}_corr.png", dpi=300)
del fig, ax, ax2