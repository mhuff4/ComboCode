"""
VISAR analysis

file: visar.py
author: JR Rygg
revisions:
    2016-1020: 
    2015-1213: substantial revisions
    2015-0217: created
"""
# TODO:
#  create VISAR analysis container (mashup of several dataimages)
#  gui should give choice of which view (image, phase, amplitude, etc)
#  optional save analysis settings in text file; load if present
#  separate form for analysis settings?
#  toggle setting annotations on/off (plot ROI, phase ROI, etc)
#  toggle showing subwindows (frequency, phase/amplitude lineouts)
#  overplot laser data?

import os
import os.path as osp
import sys
import imp
from glob import glob
from configparser import SafeConfigParser

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from PyQt4 import QtGui


try:
    imp.find_module("rypy")
except ImportError:
    sys.path.append(osp.abspath(osp.dirname(__file__) + "../../../.."))

import rypy.util.dataimage as dataimage
from rypy.util.h5nif import Hdf5
from rypy.rygui import rypl_qt
from rypy.ryplot.rycolor import get_cmap

TAU = 2 * np.pi

#============================= PARAMS =========================================
# default list of paths to search for H5 files
PATHLIST = [osp.join(os.environ['USERPROFILE'], 'Downloads'), # Win7 downloads
            osp.dirname(__file__)]

#============================= FUNCTIONS ======================================
def ftm_analyze(dat, flim=(1,-1)):
    """Analyze phase using a flavor of the Fourier transform method"""
    # TODO: f from inverse window to inverse pixels    
    
    #-- see P. M. Celliers et al, RSI, (2004), Section V --
    S = dat                   # (eq. 10), signal
    s = np.fft.fft(S, axis=0) # (eq. 11), 1-D fft of columns

    fmin, fmax = flim         # frequency filter
    if fmax is None:
        fmax = -s.shape[0]//2 # clip out all negative freqs
    d = np.zeros_like(s)
    d[fmin:fmax, :] = s[fmin:fmax, :] # filtered s
    D = np.fft.ifft(d, axis=0)        # (eq. 12) back transform

    W = np.angle(D)                   # (eq. 13) "Wrapped" phase angle...
    A = np.abs(D)                     # fringe modulation amplitude

    ratio = D[:,1:]/D[:,:-1]
    dphi = np.zeros_like(D)
    dphi[:,1:] = np.arctan(ratio.imag/np.abs(ratio)) # phase differences
    dphi2 = np.zeros_like(W)
    dphi2[:,1:] = np.angle(ratio) # angle differences

    U = -np.cumsum(dphi.real, 1) / TAU # unwrapped fringes
    U2 = -np.cumsum(dphi2, 1) / TAU    # unwrapped fringes alt

    return A, W, dphi2, U2

def mod_average(z, sigout=False):
    """Average columns of 2D array z, ignoring integer separation
    if sigmaout, return (modav, sigma) tuple
    """
    med = np.median(z, 0)
    delta = z - med
    delta = delta - np.round(delta) # np.round(delta, out=delta) # second option more mem efficient?
    modav = delta + med
    if not sigout: return modav
    else:          return (modav, np.std(delta))

def mod_ave2(z):
    """Column average, modulo 2pi. Under development"""
    x = np.mean(np.sin(TAU*z), 0) # col ave
    y = np.mean(np.cos(TAU*z), 0) # col ave
    phi = np.arctan(x/y) / TAU
    calc = (phi + np.where(y < 0, -0.5, 0) + 0.5) % 1 - 0.5
    return calc


#============================= BODY ===========================================
class VisarWindow(rypl_qt.MplWindowQT):
    def _on_open(self, fname=None):
        if fname is None:
            file_choices = ";;".join((
                "hdf5 files (*.h5; *.hdf5)",
                "Image Files (*.png; *.jpg; *.jpeg; *.bmp; *.tif; *.tiff)",
                "All files (*)",
                ))
            fname = str(QtGui.QFileDialog.getOpenFileName(self,
                            'Open file', None, file_choices))
        rypl_qt.MplWindowQT._on_open(self, fname, imageonly=True)


class VisarImage():
    """Representation of a visar image"""
    def __init__(self, fname, phase_roi=None, plot_roi=None,
                 carrier_period=0, carrier_band=(0.7, 1.4), vpf=1, leg=None):
        """
        ------------
        input params
        ------------
        fname: filename of .h5 image
        phase_roi: 4-tuple of phase roi in pixels: left, right, bottom, top
        plot_roi: 2-tuple of plot roi in pixels: bottom, top
        carrier_period: primary fringe carrier period
        carrier_band: fringe bandwidth to consider for fringe analysis
        vpf: vacuum velocity per fringe [km/s]
        leg: which streak leg: A, B, SOP
        """
        self.fname = fname
        self.phase_roi = phase_roi
        self.plot_roi = plot_roi
        self.carrier_period = carrier_period
        self.carrier_band = carrier_band
        self.vpf = vpf
        self.h5 = Hdf5(fname)
        self.image = self.h5.image
        self.leg = leg if leg is not None else 'BA'['-A-' in self.h5.identifier] # TODO: SOP

    def estimate_carrier_period(self):
        pass

    def imshow(self, **kws):
        image = self.h5.image
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial'], 'size':10})
        fig = plt.figure(figsize=(4.9, 4.9), dpi=150)
        ax = fig.add_axes((0.1, 0.1, 0.74, 0.85), xlabel="time (ns)", ylabel="pos (um)")
        cax = fig.add_axes((0.86, 0.1, 0.02, 0.85))

        vmin, vmax = dataimage.thresh_vlim(self.image, 0.01)
        locs = locals()
        [kws.update({k:locs[k]}) for k in ('vmin', 'vmax') if k not in kws]
        im = ax.imshow(self.image, extent=image.extent, aspect='auto', **kws)
        ax.set_title(self.h5.shot_id[:11] + " VISAR " + self.leg)
        ax.set_ylim(-1250, 1150)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("intensity (counts)")
        ax.yaxis.labelpad = -12
        cax.yaxis.labelpad = 2


    def plot_phase_roi(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmin, vmax = dataimage.thresh_vlim(self.image, 0.01)
        ax.imshow(self.image, aspect='equal', cmap='CMRmap', vmin=vmin, vmax=vmax)

    def plot_carrier_frequency(self):
        """Plot carrier frequency analysis of phase roi"""
        roi = self.phase_roi
        phase_slice = (slice(roi[2], roi[3]), slice(roi[0], roi[1]))
        # calculation
        S = self.image[phase_slice] # signal in phase_roi
        t_axis = self.image.x_axis[roi[0]:roi[1]] # [ns] time axis
        y_axis = self.image.y_axis[roi[2]:roi[3]] # [mic] spatial scale
        N = S.shape[0]//2

        s = np.fft.fft(S, axis=0) / N # fft
        s_abs = np.abs(s[:N,:])
        f_axis = np.arange(N) / (2 * N) # spatial frequency axis

        s_mean = np.log10(np.mean(s_abs, axis=1))
        i0 = np.argmax(s_mean[3:])
        f0 = f_axis[3+i0] # [px^-1] fringe carrier frequency (estimate)
        s0 = s_mean[3+i0]
        sys.stdout.write("{} VISAR-{} fringe period = {:.1f} px\n".format(
                         self.h5.shot_id[:11], self.leg, 1/f0))

        # plot calcs
        vlim_0 = dataimage.thresh_vlim(S, 0.01)
        vlim_1 = dataimage.thresh_vlim(np.log10(s_abs), (0.02, 0.005))
        tlim = (t_axis[0], t_axis[-1])
        ylim = (y_axis[0], y_axis[-1])
        flim = (0, 0.5) # [1/px]
        extent_0 = tlim + (0, S.shape[0]) # extent for signal
#        extent_0 = tlim + ylim # extent for signal
        extent_1 = tlim + flim # extent for fft

        # figure
        fig = plt.figure(figsize=(7,7), dpi=100)
        axs = []
        axs.append(fig.add_subplot(221, ylabel='[px]', title='signal'))
        axs.append(fig.add_subplot(222, sharey=axs[0], title='spatial lineout'))
        axs.append(fig.add_subplot(223, sharex=axs[0], title='log(fft(signal))',
                                   xlabel='time [ns]', ylabel="spatial frequency [px^-1]"))
        axs.append(fig.add_subplot(224, sharey=axs[2], xlabel='log10(power)', title='spectral lineout'))

        axs[0].imshow(S, extent=extent_0,
                      aspect='auto', vmin=vlim_0[0], vmax=vlim_0[1])
        axs[2].imshow(np.log10(s_abs), extent=extent_1,
                      aspect='auto', vmin=vlim_1[0], vmax=vlim_1[1])
        axs[1].plot(np.mean(S, axis=1), np.arange(S.shape[0]))
        axs[3].plot(s_mean, f_axis)
        axs[0].set_ylim(*extent_0[2:])
        
        axs[3].annotate("fringe period\n= {:.1f} px".format(1/f0),
                        (s0, f0), (0.95, 0.5), textcoords='axes fraction',
                        arrowprops=dict(width=1, headwidth=6, facecolor='k',
                        shrink=0.03), ha='right',)

        axs[3].axhline(f0*0.7, color='r', linestyle='dashed')
        axs[3].axhline(f0*1.4, color='r', linestyle='dashed')

        fig.tight_layout()
        fig.canvas.window().move(0,0)
        return fig

    def analyze_fringes(self, f0=None):
        """f0 = carrier frequency"""
        phase_roi, plot_roi = self.phase_roi, self.plot_roi
        roi_phase = (slice(phase_roi[2], phase_roi[3]),
                     slice(phase_roi[0], phase_roi[1]))
        roi_plot = (slice(plot_roi[0] - phase_roi[2], plot_roi[1] - phase_roi[2]),
                    slice(0, -1)) # subslice of phase_roi

        self.S = S = self.image[roi_phase] # signal

        if self.carrier_period is None:
            flim = (5, None)
        else:
            f0 = np.shape(S)[0] / self.carrier_period
            flim = (int(0.7*f0), int(1.4*f0))

        A, W, dphi2, U2 = ftm_analyze(S, flim=flim)

        self.A, self.W, self.dphi2, self.U2 = A, W, dphi2, U2
        
        f1 = mod_average(U2) # 2D
        f2 = mod_ave2(U2[roi_plot]) # 1D

        self.f1, self.f2 = f1, f2

    def plot_fringe_analysis(self):
        if not hasattr(self, "S"):
            self.analyze_fringes()
        S = self.S
        A, W, dphi2, U2 = self.A, self.W, self.dphi2, self.U2
        f1, f2 = self.f1, self.f2

        phase_roi, plot_roi = self.phase_roi, self.plot_roi
        roi_plot = (slice(plot_roi[0] - phase_roi[2], plot_roi[1] - phase_roi[2]),
                    slice(0, -1)) # subslice of phase_roi

    
        fig, axs = plt.subplots(2, 3, figsize=(12,8), sharex=True, sharey=True)
        axs[0,0].set_title("S")
        axs[1,0].set_title("A")
        axs[0,1].set_title("W")
        axs[1,1].set_title("dphi")
        axs[0,0].imshow(S,aspect='auto')
        axs[1,0].imshow(A,aspect='auto')
        axs[0,1].imshow(W,aspect='auto', cmap=get_cmap('cyc_mygbm'))
        axs[1,1].imshow(dphi2,aspect='auto')
        axs[0,2].imshow(U2,aspect='auto')
        axs[1,2].imshow(f1,aspect='auto')
    
        fig.tight_layout()
    
        plt.figure()
        plt.plot(self.vpf * np.mean(f1[roi_plot],0))
        plt.plot(self.vpf * (1+np.mean(f1[roi_plot],0)))
        plt.plot(self.vpf * f2)
        
    def print_summary(self):
        print("filename", self.fname)
        print("shot_id", self.h5.shot_id)
        print("location", self.h5.location)
        print("unit", self.h5.unit)
        print("identifier", self.h5.identifier)
        print("leg", self.leg)
        print("data_id", self.h5.data_id)
        print("data_taxon", self.h5.data_taxon)
        print("shape", self.h5.shape)
        print("min/max", np.min(self.h5.image), np.max(self.h5.image))
        print("---")


class VisarSet():
    """Represent set of NIF VISAR images"""
    LEGNAMES = ('A', 'B', 'SOP') # LEGS of NIF VISAR system
    DEFAULTS = dict( # default settings if no .ini file
        ND_A=0, vpf_A=1, carrier_A=0, band_A=np.array((0.7, 1.4)),
        ND_B=0, vpf_B=1, carrier_B=0, band_B=np.array((0.7, 1.4)),
        ND_SOP=0,
        roi_A=(50, -50, 300, -300),
        roi_B=(50, -50, 300, -300),
        N_regions=0, N_events_A=0, N_events_B=0, j0=0,
    )
    # data_id preference order for NIF VISAR images
    PREF_ORDER = ("TPS-CROSS-TIMING", "CROSS-TIMING", "RAW_STREAK_CAMERA")
        
    def __init__(self, handle, pathlist=None, legs=None):
        self.handle = handle
        self.legnames = legs if legs is not None else list(self.LEGNAMES)
        self.filenames = self.find_files(handle, pathlist)
        self.settings = self.parse_ini()
        self.leg_kws = self.parse_settings()
        self.legs = [VisarImage(**self.leg_kws[leg]) for leg in self.legnames]

    def __str__(self):
        return "<VisarSet: {}>".format(self.handle)

    def __len__(self):
        return len(self.legnames)

    def __iter__(self):
        for leg in self.legnames:
            yield leg

    def plot_carrier_frequency(self):
        for i, leg in enumerate(self.legnames):
            if leg == 'SOP':
                continue
            fig = self.legs[i].plot_carrier_frequency()
            fig.canvas.window().move(i*716, 0)
            fig.canvas.set_window_title("VISAR leg {}".format(leg))

    def find_files(self, handle, pathlist=None):
        # modify handle and pathlist
        if handle[-2] in '.-': # convert shorthand shotnum (e.g. 'N151213.1')
            handle = handle[:-2] + '-00' + handle[-1]

        pathlist = pathlist if pathlist is not None else []
        pathlist.extend(PATHLIST)
        
        # search pathlist for files containing shotnum
        filelist = []
        inilist = []
        for path in pathlist:
            filelist.extend(glob(osp.join(path, "*{}*.h5".format(handle))))
            inilist.extend(glob(osp.join(path, "*{}*.ini".format(handle))))

        # sort files; exclude files without "VISAR" in the filename
        filelist = [f for f in filelist if "VISAR" in f]
        sublists = []
        for leg in self.legnames:
            sublists.append([f for f in filelist if "-{}-".format(leg) in f])
        
        # construct dictionary with preferred files
        filedict = {}
        for i, leg in enumerate(self.legnames):
            filedict[leg] = self.choose_file(sublists[i])
        filedict['ini'] = self.choose_file(inilist, [osp.dirname(filedict['A'])])
        self.legnames = [leg for leg in self.legnames if filedict[leg] is not None]
        
        return filedict

    def choose_file(self, filelist, pref_order=None):
        """Return filename in filelist containing earliest string in pref_order"""
        if len(filelist) == 0:
            return None
        if len(filelist) > 1:
            pref_order = pref_order if pref_order is None else self.PREF_ORDER
            for p in pref_order:
                for f in filelist:
                    if p in f:
                        return f

        return filelist[0]

    def parse_ini(self):
        """Read ini file for visar setup and analysis settings."""
        cp = SafeConfigParser(defaults=self.DEFAULTS)
        if self.filenames['ini'] is not None:
            cp.read(self.filenames['ini'])
        return cp

    def parse_settings(self):
        """parse settings for VisarImage kws."""
        leg_kws = {leg: {'fname': self.filenames[leg]} for leg in self.legnames}

        if self.settings is None:
            return leg_kws

        cp = self.settings
        for leg in self.legnames:
            kws = leg_kws[leg]
            kws['leg'] = leg
            if leg == 'SOP':
                continue
            roi = cp.get('phase', 'roi_' + leg)
            if roi in (None, 'None', 0, '0'):
                kws['phase_roi'] = self.DEFAULTS['roi_' + leg]
                print(kws['phase_roi'])
            else:
                kws['phase_roi'] = np.array(roi.split(','), dtype=int)
            roi = cp.get('regions', 'roi_' + leg + '2')
            if roi in (None, 'None', 0, '0'):
                kws['plot_roi'] = np.array((0, -1))
            else:
                kws['plot_roi'] = np.array(roi.split(','), dtype=int)


            kws['carrier_period'] = cp.getfloat('phase', 'carrier_' + leg)
            kws['vpf'] = cp.getfloat('setup', 'vpf_' + leg)

        return leg_kws


    def plot_fringe_analysis(self):
        for i, leg in enumerate(self.legnames):
            if leg == 'SOP':
                continue
            self.legs[i].analyze_fringes()

        visA = self.legs[0]
        visB = self.legs[1]

        phase_roi, plot_roi = visA.phase_roi, visA.plot_roi
        tA = visA.h5.image.x_axis[phase_roi[0]:phase_roi[1]]
        roi_plot_A = (slice(plot_roi[0] - phase_roi[2], plot_roi[1] - phase_roi[2]),
                     slice(0, None)) 
        phase_roi, plot_roi = visB.phase_roi, visB.plot_roi
        tB = visB.h5.image.x_axis[phase_roi[0]:phase_roi[1]]
        roi_plot_B = (slice(plot_roi[0] - phase_roi[2], plot_roi[1] - phase_roi[2]),
                     slice(0, None)) 

        n = 1.55 # refractive index

        vpfA, vpfB = visA.vpf/n, visB.vpf/n
        v0A = vpfA * np.mean(visA.f1[roi_plot_A],0)+0.5
        v0B = vpfB * np.mean(visB.f1[roi_plot_B],0)

        
        fig = plt.figure()        
        ax = fig.add_subplot(111)
        
        for i in range(0,10):
            ax.plot(tA, v0A + i*vpfA, 'r')
        ax.plot(tB, v0B, 'b')
        ax.plot(tB, v0B+vpfB, 'b')
        ax.plot(tB, v0B+2*vpfB, 'b')    
        ax.plot(tB, v0B+3*vpfB, 'b')
        ax.plot(tB, v0B+4*vpfB, 'b')
        

#============================= MAIN ===========================================
def visar_summary():
    shotnum = 'N151214-003'
    vs = VisarSet(shotnum, legs=['A', 'B'])
    
    # leg A
    vs.legs[0].imshow(cmap='coolwarm')

    # leg B
    vs.legs[1].imshow(cmap='coolwarm')

    plt.show()


def visar_gui(fname=None):
    app = QtGui.QApplication(sys.argv)
    w1 = VisarWindow()
    w1._on_open(fname)
    w1.show()
    sys.exit(app.exec_())


def test_visar(shotnum='N151214-003'):
    """test VISAR analysis"""
    vs = VisarSet(shotnum, legs=['A', 'B'])

#    i = 1
#    visar_gui(fname=vs.filenames[vs.legnames[i]]) # check phase and region roi
#    vs.plot_carrier_frequency()
    
#    vs.legs[0].plot_fringe_analysis()
#    vs.legs[1].plot_fringe_analysis()
#    vs.plot_fringe_analysis()


#    folder = r"C:\usr\zchive\Data\NIF\TARDIS\VISAR"
#    fnameA = osp.join(folder, "TD_TC090-315_VISAR_STREAK-A-01-DB_SHOT_TPS-CROSS-TIMING_{}-999.h5".format(shotnum))
#    fnameB = osp.join(folder, "TD_TC090-315_VISAR_STREAK-B-01-DB_SHOT_TPS-CROSS-TIMING_{}-999.h5".format(shotnum))
#
#    dat = {
#        'N151129-002': dict(
#            vpfA=2.1857, vpfB=5.4604,
#            roi_phsA=(81, 1311, 458, 919,), roi_phsB=(),
#            roi_A1=(540,841), roiA2=(),
#            f0_A=17, f0_B=17, # carrier frequency (tied to choice of phase roi)
#            p0_A=27.1, p0_B=27, # carrier period (pixels)
#        ),
#    }.get(shotnum)
#    
#    visA = VisarImage(fnameA, vpf=dat['vpfA'],
#                      phase_roi=dat['roi_phsA'], plot_roi=dat['roi_A1'])
#    visA.print_summary()
#    visA.plot_carrier_frequency()

#    vs.legs[0].analyze_fringes(17)
#    visA.plot_fringe_analysis()
    vs.legs[0].h5.imshow(cmap='coolwarm')
    plt.xlabel("time (ns)")
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    visar_summary()
#    test_visar()
#    visar_gui()
#    si = ryxlsx.ShotInfo(XL_TARDIS_OVERVIEW, 'N140314',
#                         datasets=("setup", "target", "results"))
#    d = si.target
#    for k in sorted(d):
#        print k, d[k]

