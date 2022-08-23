#!/home/seis/anaconda2/bin/python

import numpy as np

def noise_floor(x):
    """
    Quick and dirty approximation of median envelope amplitude 
    aka "noise floor" which uses the half range of the 2nd to 
    98th percentile amplitudes.
    x: np array or list of numbers
    """
    xsort = np.sort(x)
    x2 = xsort[int(len(x)*0.02)]
    x98 = xsort[int(len(x)*0.98)]
    noisefloor = (x98-x2)/2.
    return noisefloor


def duration_exceed_RMS(x, ampthresh, RMSlen, dt):
    """
    Noise level metric that calcualtes the total duration of the RMS
    function above a theshold (ampthresh).
    RMSlen = window length used to calculate RMS (sec)
    x: np array or list of numbers
    st = an ObsPy stream
    dt = timeseries increment (sec)
    """
    from obspy.signal.util import smooth
    duration = 0
    if ( len(x) > int(RMSlen/dt) ):
        iRMSwinlen = int(RMSlen/dt)
        RMS = np.sqrt(smooth((x**2),iRMSwinlen))
        duration = ((RMS > ampthresh).sum())*dt
    else:
        duration = []
        duration = -1
#        print "Error duration_exceed_RMS: len(x)=" + str(len(x)*dt) + " must be greater than RMSlen=" + str(RMSlen)

    return duration


def count_peaks_stalta(x, y, sta, lta, mpd, mph, dt, ampthresh):
    """
    Counts peaks of a passed-through STA/LTA function (y) where 
       the absolute amplitude of the original timeseries (x) exceeds 
       an amplitude (ampthresh).
    x: np array or list of numbers of original time series
    y: np array or list of numbers of STA/LTA function of x
    sta: short term average used (sec)
    lta: long term average used (sec)
    mpd: minimum distance between peaks should be > sta+lta (sec)
    mph: minimum peak height of the stalta function
    ampthresh: absolute amplitude of the biggest local peak
    dt: timeseries increment (sec)
    """
    if ( len(y) > int(mpd/dt) and mpd > (sta+lta) ):
        #  Find peaks
        indexlocalmax = detect_peaks(y,mpd=int(mpd/dt),mph=mph)
        # Find the local maxima in x within sta+lta sec of peak
        peaksnr = np.zeros(len(indexlocalmax))
        peakamp = np.zeros(len(indexlocalmax))
        istalta = int((sta+lta)/dt)
        impd = int((mpd/dt))
        for j in range(0,len(indexlocalmax)):
            ipeak1 = max(0,indexlocalmax[j]-istalta)
            ipeak2 = min(indexlocalmax[j]+istalta,len(x))
            peaksnr[j] = max(y[ipeak1:ipeak2])
            peakamp[j] = max(abs(x[ipeak1:ipeak2]))
            # In case the x has zero mean but long period noise, use local demeaning
#            ilta1 = max(0,indexlocalmax[j]-impd-istalta)
#            ilta2 = min(indexlocalmax[j]-istalta,len(x))
#            peakamp[j] = max(abs(x[ipeak1:ipeak2]-np.mean(x[ilta1:ilta2])))
        peakcount = (peakamp>ampthresh).sum()
    else:
        peakcount = []
        peakcount = -1
#        print "Error count_peaks_stalta: len(x)=" + str(len(x)*dt) + " must be > sta+lta=" + str(sta+lta)
    return peakcount


def count_peaks_stalta_new(x0, y0, sta, lta, mpd, mph, dt, twin, ampthresh):
    """
    Counts peaks of a passed-through STA/LTA function (y) where 
       the absolute amplitude of the original timeseries (x) exceeds 
       an amplitude (ampthresh).
    x: np array or list of numbers of original time series
    y: np array or list of numbers of STA/LTA function of x
    sta: short term average used (sec)
    lta: long term average used (sec)
    mpd: minimum distance between peaks should be > sta+lta (sec)
    mph: minimum peak height of the stalta function
    twin: window length following trigger to measure amp (sec)
    ampthresh: absolute amplitude of the biggest local peak
    dt: timeseries increment (sec)
    """
    x = np.copy(x0)
    y = np.copy(y0)
    if ( len(y) > int(mpd/dt) and mpd > (sta+lta)*0 ):
        #  Get indicies of values breaching mph
        peakcount = 0
        igtmph = []
        y = np.asarray(y)
        y[ y > mph ] = mph
        y[ y < mph ] = 0
        ydiff = np.diff(y)
        ydiff[ ydiff < 0 ] = 0
        iydiff = np.nonzero(ydiff)[0].tolist()
        impd = int(mpd/dt)
        itwin = int(twin/dt)
        ista = int(sta/dt)
        # nuke any sta/lta > mph that doesn't have a max(amp(x)) > ampthresh
        if ( len(iydiff) > 0 ):
            x = abs(x)
            for i in range(0,len(iydiff)):
                i1 = iydiff[i] - ista
                i2 = min(i1+itwin,len(x))
                if ( max(x[i1:i2]) < ampthresh ):
                    iydiff[i] = 0
        iydiff = list(filter(None,iydiff))
        #  nuke any triggers within mpd sec of another trigger
        if ( len(iydiff) > 0 ):
           for i in range(0,len(iydiff)):
               for j in range(i+1,len(iydiff)):
                   if ( (iydiff[j] - iydiff[i]) < impd ):
                       iydiff[j] = 0
        iydiff = list(filter(None,iydiff))
        # count valid triggers
        for i in range(0,len(iydiff)):
            i1 = iydiff[i] - ista
            i2 = min(i1+itwin,len(x))
            if ( max(x[i1:i2]) > ampthresh ):
                peakcount = peakcount + 1
    else:
        peakcount = []
        peakcount = -1
#        print "Error count_peaks_stalta: len(x)=" + str(len(x)*dt) + " must be > sta+lta=" + str(sta+lta)
    return peakcount


def count_peaks_stalta_Elarms(xA0, xV0, xD0, y0, sta, lta, mpd, mph, dt, twin, ch0):
    """
    Counts peaks of a passed-through STA/LTA function (y) where 
       the absolute amplitude of the original timeseries (x) exceeds 
       an amplitude (ampthresh).
    x: np array or list of numbers of original time series
    y: np array or list of numbers of STA/LTA function of x
    sta: short term average used (sec)
    lta: long term average used (sec)
    mpd: minimum distance between peaks should be > sta+lta (sec)
    mph: minimum peak height of the stalta function
    twin: window length following trigger to measure amp (sec)
    ampthresh: absolute amplitude of the biggest local peak
    dt: timeseries increment (sec)
    """
    AorV = 0
    if ( ch0[1:2] == 'N' ):
        AorV = 1
    if ( ch0[1:2] == 'H' ):
        AorV = 2
    ampboxA = 0.000022
    ampboxV = 0.000000022
    # ln log amplitudes
#    ampAmin = 0.00082
#    ampVmin = 0.0000409
#    ampVmax = 0.2
#    ampDmin = 0.0000409
#    ampDmax = 0.331
    # log10 amplitudes
    ampAmin = 0.000031623
    ampVmin = 0.000000031623
    ampVmax = 10
    ampDmin = 0.000000031623
    ampDmax = 31.623

    xA = np.copy(xA0)
    xV = np.copy(xV0)
    xD = np.copy(xD0)
    y = np.copy(y0)
    if ( len(y) > int(mpd/dt) and mpd > (sta+lta)*0 ):
        #  Get indicies of values breaching mph
        peakcount = 0
        boxcount = 0
        igtmph = []
        iydiff = []
        y = np.asarray(y)
        y[ y > mph ] = mph
        y[ y < mph ] = 0
        ydiff = np.diff(y)
        ydiff[ ydiff < 0 ] = 0
        iydiff = np.nonzero(ydiff)[0].tolist()
        impd = int(mpd/dt)
        itwin = int(twin/dt)
        ista = int(sta/dt)
        ibox = int(0.1/dt)
        # nuke any sta/lta > mph that doesn't have a max(amp(x)) > ampthresh
        if ( len(iydiff) > 0 ):
            xA = abs(xA)
            xV = abs(xV)
            xD = abs(xD)
            for i in range(0,len(iydiff)):
                i1 = iydiff[i] - ista 
#                i1 = iydiff[i] - ista + ibox
                i2 = min(i1+itwin,len(xA))
                ib1 = min(iydiff[i] + 1 + ibox,len(xA)-1)
                ib2 = min(ib1 + ibox,len(xA))
                maxA = max(xA[i1:i2])
                maxV = max(xV[i1:i2])
                maxD = max(xD[i1:i2])
#                print ("TRIG : " + str(y0[i1+ista+1]) + " " + str(maxA/ampAmin) + " " + str(maxV/ampVmin) + " " + str(maxD/ampDmin)  + "  " + str(maxD/ampDmax)  )
                if ( maxA < ampAmin or maxV < ampVmin or maxV > ampVmax or maxD < ampDmin or maxD > ampDmax ):
                        iydiff[i] = 0
                else:
                    if ( AorV == 1 ):
                        boxA = max(xA[ib1:ib2]) - min(xA[ib1:ib2])
                        if ( boxA < ampboxA ):
#                            print ("BOXCAR: A " + str(boxA/100) + " " + str(boxA/ampboxA) + " " + str(y0[iydiff[i]])  )
                            iydiff[i] = 0
                            boxcount = boxcount + 1
                    if ( AorV == 2 ):
                        boxV = max(xV[ib1:ib2]) - min(xV[ib1:ib2])
                        if ( boxV < ampboxV ):
#                            print ("BOXCAR: V " + str(boxV/100) + " " + str(boxV/ampboxV) + " " + str(y0[iydiff[i]])   )
                            iydiff[i] = 0
                            boxcount = boxcount + 1
        iydiff = list(filter(None,iydiff))
        #  nuke any triggers within mpd sec of another trigger
        if ( len(iydiff) > 0 ):
           for i in range(0,len(iydiff)):
               for j in range(i+1,len(iydiff)):
                   if ( (iydiff[j] - iydiff[i]) < impd ):
                       iydiff[j] = 0
        iydiff = list(filter(None,iydiff))
        peakcount = len(iydiff)
    else:
        peakcount = []
        peakcount = -1
        boxcount = -1
#        print "Error count_peaks_stalta: len(x)=" + str(len(xA)*dt) + " must be > sta+lta=" + str(sta+lta)
    return [peakcount,boxcount]


def count_peaks_stalta_Elarms_times(xA0, xV0, xD0, y0, sta, lta, mpd, mph, dt, twin, ch0):
    """
    Counts peaks of a passed-through STA/LTA function (y) where 
       the absolute amplitude of the original timeseries (x) exceeds 
       an amplitude (ampthresh).
    x: np array or list of numbers of original time series
    y: np array or list of numbers of STA/LTA function of x
    sta: short term average used (sec)
    lta: long term average used (sec)
    mpd: minimum distance between peaks should be > sta+lta (sec)
    mph: minimum peak height of the stalta function
    twin: window length following trigger to measure amp (sec)
    ampthresh: absolute amplitude of the biggest local peak
    dt: timeseries increment (sec)
    """
    AorV = 0
    if ( ch0[1:2] == 'N' ):
        AorV = 1
    if ( ch0[1:2] == 'H' ):
        AorV = 2
    ampboxA = 0.000022
    ampboxV = 0.000000022
    # ln log amplitudes
#    ampAmin = 0.00082
#    ampVmin = 0.0000409
#    ampVmax = 0.2
#    ampDmin = 0.0000409
#    ampDmax = 0.331
    # log10 amplitudes
    ampAmin = 0.000031623
    ampVmin = 0.000000031623
    ampVmax = 10
    ampDmin = 0.000000031623
    ampDmax = 31.623

    xA = np.copy(xA0)
    xV = np.copy(xV0)
    xD = np.copy(xD0)
    y = np.copy(y0)
    iydiff = []
    if ( len(y) > int(mpd/dt) and mpd > (sta+lta)*0 ):
        #  Get indicies of values breaching mph
        peakcount = 0
        boxcount = 0
        igtmph = []
        iydiff = []
        y = np.asarray(y)
        y[ y > mph ] = mph
        y[ y < mph ] = 0
        ydiff = np.diff(y)
        ydiff[ ydiff < 0 ] = 0
        iydiff = np.nonzero(ydiff)[0].tolist()
        impd = int(mpd/dt)
        itwin = int(twin/dt)
        ista = int(sta/dt)
        ibox = int(0.1/dt)
        # nuke any sta/lta > mph that doesn't have a max(amp(x)) > ampthresh
        if ( len(iydiff) > 0 ):
            xA = abs(xA)
            xV = abs(xV)
            xD = abs(xD)
            for i in range(0,len(iydiff)):
                i1 = iydiff[i] - ista 
#                i1 = iydiff[i] - ista + ibox
                i2 = min(i1+itwin,len(xA))
                ib1 = min(iydiff[i] + 1 + ibox,len(xA)-1)
                ib2 = min(ib1 + ibox,len(xA))
                maxA = max(xA[i1:i2])
                maxV = max(xV[i1:i2])
                maxD = max(xD[i1:i2])
#                print ("TRIG : " + str(y0[i1+ista+1]) + " " + str(maxA/ampAmin) + " " + str(maxV/ampVmin) + " " + str(maxD/ampDmin)  + "  " + str(maxD/ampDmax)  )
                if ( maxA < ampAmin or maxV < ampVmin or maxV > ampVmax or maxD < ampDmin or maxD > ampDmax ):
                        iydiff[i] = 0
                else:
                    if ( AorV == 1 ):
                        boxA = max(xA[ib1:ib2]) - min(xA[ib1:ib2])
                        if ( boxA < ampboxA ):
#                            print ("BOXCAR: A " + str(boxA/100) + " " + str(boxA/ampboxA) + " " + str(y0[iydiff[i]])  )
                            iydiff[i] = 0
                            boxcount = boxcount + 1
                    if ( AorV == 2 ):
                        boxV = max(xV[ib1:ib2]) - min(xV[ib1:ib2])
                        if ( boxV < ampboxV ):
#                            print ("BOXCAR: V " + str(boxV/100) + " " + str(boxV/ampboxV) + " " + str(y0[iydiff[i]])   )
                            iydiff[i] = 0
                            boxcount = boxcount + 1
        iydiff = list(filter(None,iydiff))
        #  nuke any triggers within mpd sec of another trigger
        if ( len(iydiff) > 0 ):
           for i in range(0,len(iydiff)):
               for j in range(i+1,len(iydiff)):
                   if ( (iydiff[j] - iydiff[i]) < impd ):
                       iydiff[j] = 0
        iydiff = list(filter(None,iydiff))
        peakcount = len(iydiff)
    else:
        peakcount = []
        peakcount = -1
        boxcount = -1
#        print "Error count_peaks_stalta: len(x)=" + str(len(xA)*dt) + " must be > sta+lta=" + str(sta+lta)
    return [peakcount,boxcount,iydiff]



def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',kpsh=False, valley=False):
    """ Detect peaks in data based on their amplitude and other features.
    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 
    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])
    return ind


def get_power(trPower,inv,periodlist):
    """
    Calculates the PSD using ObsPy PPSD.
    trPower: input ObsPy trace that must be > 3600. sec long.
    inv:  ObsPy inventory.  Must be level = 'response'
    periodlist: a list of periods at which you want the power 
    returns: power in dB at the periods in periodlist
    """
    from obspy.signal import PPSD
    powers_at_periods = []
    if ( trPower.stats.delta * trPower.stats.npts > 3600. ):
        ppsd = PPSD(trPower.stats, metadata = inv)
        ppsd.add(trPower)
        if ( len(ppsd._binned_psds) > 0 ):
            psd_periods = ppsd._period_binning[2]
            psd_power = []
            psd_power = ppsd._binned_psds[0]
            for i in range(0,len(periodlist)):
                powers_at_periods.append(psd_power[np.argmin(abs(psd_periods-periodlist[i]))])
        else:
            for i in range(0,len(periodlist)):
                powers_at_periods.append(-1)

    return powers_at_periods


def count_triggers_FinDer(xA0, twin, ampthresh, dt):
    """
    Counts the number of instances of amplitude exceedance by convertnig input
    numpy array into a bit (trace of 0s or 1s) based on abs(amplitude) > mph.  
    xA0 = input numpy array of gain corrected acceleration amplitudes.
    twin = minimum time in sec. between counted instances of the peak amplitude 
           being exceed. 
    ampthresh = threshold above which amplitudes are considered, all others
                are zeroed.
    dt = sampling interval.
    Returns: integer count of number instances exceeding amplitude while being
             separated by at least twin seconds.
    """
    xA = abs(np.copy(xA0))
    xA[ xA < ampthresh ] = 0
    xA[ xA >= ampthresh ] = 1
    iXA = np.nonzero(xA)[0]
    # nuke all instances of xA = 1 that are within twin of a counted previous one
    if ( len(iXA) > 0 ):
        itwin = int(twin/dt)
#        iXB = iXA.copy()
#        for i in range(0,len(iXB)):
#            for j in range(i+1,len(iXB)):
#                if ( (iXB[j] - iXB[i]) < itwin and iXB[j] > 0 ):
#                    iXB[j] = 0
#        iXB = np.nonzero(iXB)[0]
        i = 0
        while i < len(iXA):
            j1 = iXA[i]
            j2 = iXA[i] + itwin -1
            if ( iXA[i] > 0 ):
                iXA[ (iXA > j1) & (iXA < j2) ] = 0
            i = i + 1
        iXA = np.nonzero(iXA)[0]
    return len(iXA)


