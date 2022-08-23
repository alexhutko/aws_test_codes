#!/home/ahutko/miniconda3/envs/aws/bin/python

#
# Calculate some station metrics - noise floor, power at different frequencies, Nspikes...
# This is the object oriented version.

import os
import sys
from obspy import read
import numpy as np
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.signal.util import smooth
import datetime
from datetime import date
import time
import timeit
from obspy.signal.trigger import z_detect
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.filter import envelope 
from obspy.signal.filter import bandpass
from obspy.signal import PPSD
from get_data_metadata import *
from noise_metrics import *
from plot_station import *
from obspy import Stream

try:
    from squacapi_client.models.write_only_measurement_serializer \
    import WriteOnlyMeasurementSerializer

    from squacapi_client.pnsn_utilities \
    import get_client, make_channel_map, make_metric_map, perform_bulk_create

    no_squacapi = False
except Exception as e:
    print("Info: squacapi_client not available, cannot use --squac option")
    no_squacapi = True

# LOCAL, 0.05GING, PRODUCTION
#ENVIRONMENT = '0.05GING'
#HOST = 'https://squacapi.pnsn.org'
#HOSTstaging = 'https://staging-squacapi.pnsn.org'
HOSTproduction = 'https://squacapi.pnsn.org'

try:
    USER = os.environ['SQUACAPI_USER']
    PASSWORD = os.environ['SQUACAPI_PASSWD']
#    HOST = os.environ[ENVIRONMENT + '_HOST']
except KeyError:
    sys.exit("Requires ENV vars SQUACAPI_USER, SQUACAPI_PASSWD and HOST")

# create API client and retrieve foreign keys from the SQUAC database
#squac_client_staging = get_client(USER, PASSWORD, HOSTstaging)
squac_client_production = get_client(USER, PASSWORD, HOSTproduction)
#measurements_staging = []
measurements_production = []
print("SQUAC client: ", squac_client_production)
#squac_channels_staging = squac_client_staging.v1_0_nslc_channels_list()
#squac_metrics_staging = squac_client_staging.v1_0_measurement_metrics_list()
#channel_map_staging = make_channel_map(squac_channels_staging)
#metric_map_staging = make_metric_map(squac_metrics_staging)
#squac_channels_production = squac_client_production.v1_0_nslc_channels_list()
squac_metrics_production = squac_client_production.v1_0_measurement_metrics_list()
#channel_map_production = make_channel_map(squac_channels_production)
metric_map_production = make_metric_map(squac_metrics_production)
#nsquac_staging = 0
nsquac_production = 0

#----- Set up variables
chanfile = "chanfile.PNSN.test"
datacenter = "IRIS"
sta = 0.05
lta = 5.0
padding = 120
duration = 3600
starttime = datetime.datetime(2022,8,23,17,0,0)
endtime = starttime + datetime.timedelta(0,duration)
starttimesquac = starttime
endtimesquac = endtime
print("Starttime: " + str(starttime))
freqBP1 = 0.06
freqBP2 = 0.075
freqBP3 = 15.0
freqBP4 = 17.0
freqHP = 3.0
mpd = 10.0
twin = 4.0

RMSlen = 5.0
GainOrFullResp = "Gain"
#PSDperiods = [100,90,80,70,60,50,40,30,20,15,12,10,9,8,7,6,5,4,3,2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01 ]
PSDperiods = [50, 30, 20, 12, 10, 5, 1, 0.2, 0.1, 0.05, 0.02 ]

if ( "esp" in GainOrFullResp or "ESP" in GainOrFullResp ):
    ifullresp = 1
else: 
    ifullresp = 0
iplot = 0   #---- careful: of order 25sec/channel for 1 hr long traces
try:
    dbuser = os.environ['POSTGRES_USER']
    dbpass = os.environ['POSTGRES_PASSWD']
except:
    dbuser = None
    dbpass = None

#----- Set up the FDSN client
try:
    client = Client(datacenter,timeout=800,debug=True)
except:
    print("Failed client connection")
    exit()

#----- Download waveforms using channel list
tbuffer = sta + lta
Time1padded = starttime - datetime.timedelta(0,tbuffer+120)
stAll,inv = download_waveforms_metadata_fdsn_bulk(chanfile,Time1padded,duration+tbuffer+(2*padding),client)

#----- Get the time slices we want for analysis.
Time1 = starttime - datetime.timedelta(0,tbuffer)
Time2 = starttime + datetime.timedelta(0,duration)
lenrequested = (Time2-Time1).seconds

#-------- open file for writing trigger times
try:
    suffix = chanfile.split('.')[1]
except:
    suffix = 'temp'
outfile = open('LOGFILES/log.TRIGGERS.2022.8.23.17.' + suffix,"w")
outfileBB = open('LOGFILES/log.TRIGGERS_BB.2022.8.23.17.' + suffix,"w")

#----- Make channel map dictionary for squac_ids
f = open('channels_squacids_west_coast')
channel_map_production = {}
for line in f.readlines():
    sncl = line.split()[0]
    try:
        squacid = int(line.split()[1])
        channel_map_production[sncl] = int(squacid)
    except:
        pass
f.close()

#----- Make dictionary of metric names from this code to squac names

squac_metric_ids = {}
squac_metric_ids['rawmin'] = 82  # hourly_min
squac_metric_ids['rawmax'] = 83  # hourly_max
squac_metric_ids['rawmean'] = 84  # hourly_mean
squac_metric_ids['rawrange'] = 85  # hourly_range
squac_metric_ids['finder_2cm_hp'] = 86 # acc_gt_2.0
squac_metric_ids['RMSduration_0p07cm'] = 87  # rms__bp_above_.07     TYPO in squac name!!!!
squac_metric_ids['RMSduration_0p07cmHP'] = 88  # rms_above_.07
squac_metric_ids['snr20_0p34cm'] = 89  # acc_bp_spikes_gt_.34
squac_metric_ids['snr20_0p34cmHP'] = 90  # acc_spikes_gt_.34
squac_metric_ids['NTrigElarmSAlex'] = 91  # approximate_epic_triggers
squac_metric_ids['NTrigElarmSAlexBB'] = 92  # approximate_epic_bp_triggers
squac_metric_ids['NoiseFloorAcc'] = 93  # hourly_noise_floor_bp_acc
squac_metric_ids['accmax'] = 94  # hourly_max_bp_acc

#----- Loop over stations.  If a station has gaps, a stream of that sncl will 
#-----      be made from all the traces.

itrace = 0
T0 = timeit.default_timer()
Tsum = 0
for i in range(0,len(stAll)):
    ntr = 1
    iduplicate = 0
    if ( stAll[i].stats.location == "" ):
        sncl = stAll[i].stats.network + "." + stAll[i].stats.station + ".--." + stAll[i].stats.channel
        scnl = stAll[i].stats.station + "." + stAll[i].stats.channel + "." + stAll[i].stats.network + ".--"
    else:
        sncl = stAll[i].stats.network + "." + stAll[i].stats.station + "." + stAll[i].stats.location + "." + stAll[i].stats.channel
        scnl = stAll[i].stats.station + "." + stAll[i].stats.channel + "." + stAll[i].stats.network + "." + stAll[i].stats.location
#    sncl_id = sncl_list.index(sncl) + 1   #--- added 
    isnclfail = 0
#    try:
#        sncl_id = sncl_id_list[sncl_list.index(sncl)]   #--- added Thu Aug  9 17:54:41 PDT 2018
#    except:
#        isnclfail = 1
    datasrc_id = 1
#    print("SNCL " + sncl + "  sncl_id " + str(sncl_id))
    for j in range(0,len(stAll)):
        if ( stAll[j].stats.location == "" ):
            sncl2 = stAll[j].stats.network + "." + stAll[j].stats.station + ".--." + stAll[j].stats.channel
        else:
            sncl2 = stAll[j].stats.network + "." + stAll[j].stats.station + "." + stAll[j].stats.location + "." + stAll[j].stats.channel
        if ( sncl == sncl2 and i != j ):
            if ( j < i ):
                iduplicate = 1
                break
            else:
                ntr = ntr + 1
    if ( iduplicate == 0 and isnclfail == 0 ):
        itrace = itrace + 1
        tr = []
        trRaw = []
        trVelTrig = []
        trDisHighPass = []
        trVelHighPass = []
        trAccHighPass = []
        trAccFilt = []
        trVelFilt = []
        trDisFilt = []
        trAccFilt10Hz = []
        trVelFilt10Hz = []
        trDisFilt10Hz = []
        trAccFilt20Hz = []
        trVelFilt20Hz = []
        trDisFilt20Hz = []
        trAccFilt25Hz = []
        trVelFilt25Hz = []
        trDisFilt25Hz = []
        stalta = []
        dataraw = []
        dataVelTrig = []
        dataDisHighPass = []
        dataVelHighPass = []
        dataAccHighPass = []
        dataAccFilt = []
        dataAccFilt10Hz = []
        dataAccFilt20Hz = []
        dataAccFilt25Hz = []
        dataVelFilt = []
        dataVelFilt10Hz = []
        dataVelFilt20Hz = []
        dataVelFilt25Hz = []
        dataDisFilt = []
        dataDisFilt10Hz = []
        dataDisFilt20Hz = []
        dataDisFilt25Hz = []
        datastalta = []
        nptstotal = 0
        #inv = download_metadata_fdsn(stAll[i].stats.network, stAll[i].stats.station, stAll[i].stats.location, stAll[i].stats.channel, Time1, client)
        dt = stAll[i].stats.delta
        freqBB3 = 0.9 * (1/(2*dt))
        freqBB4 = 0.99 * (1/(2*dt))

        powers = np.zeros(len(PSDperiods))
        segmentlong = 0
        segmentshort = 9e6
        for j in range(i,i+ntr):
            trRaw = stAll[j].copy()
            trRaw = trRaw.slice(UTCDateTime(Time1),UTCDateTime(Time2))
            npts = trRaw.stats.npts
            nptstotal = nptstotal + npts
            trPadded = stAll[j].copy()
            trPadded = trPadded.detrend(type='polynomial',order=3)
            if ( dt*npts > segmentlong ):
                segmentlong = dt*npts
            if ( dt*npts < segmentshort ):
                segmentshort = dt*npts
            trVelTrig = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,0,freqHP,0,0,inv).copy()
            trDisHighPass = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Dis",0,0,freqBP2,0,0,inv).copy()
            trVelHighPass = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,0,freqBP2,0,0,inv).copy()
            trAccHighPass = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Acc",0,0,freqBP2,0,0,inv).copy()
            trAccFilt = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Acc",0,freqBP1,freqBP2,freqBP3,freqBP4,inv).copy()
            trVelFilt = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,freqBP1,freqBP2,freqBP3,freqBP4,inv).copy()
            trDisFilt = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Dis",0,freqBP1,freqBP2,freqBP3,freqBP4,inv).copy()
#            trAccFilt10Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Acc",0,freqBP1,freqBP2,10.0,freqBP4,inv).copy()
#            trVelFilt10Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,freqBP1,freqBP2,10.0,freqBP4,inv).copy()
#            trDisFilt10Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Dis",0,freqBP1,freqBP2,10.0,freqBP4,inv).copy()
#            trAccFilt20Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Acc",0,freqBP1,freqBP2,20.0,freqBP4,inv).copy()
#            trVelFilt20Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,freqBP1,freqBP2,20.0,freqBP4,inv).copy()
#            trDisFilt20Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Dis",0,freqBP1,freqBP2,20.0,freqBP4,inv).copy()
#            trAccFilt25Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Acc",0,freqBP1,freqBP2,25.0,freqBP4,inv).copy()
#            trVelFilt25Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Vel",0,freqBP1,freqBP2,25.0,freqBP4,inv).copy()
#            trDisFilt25Hz = raw_trace_to_ground_motion_filtered_pruned(trPadded,Time1,Time2,"Dis",0,freqBP1,freqBP2,25.0,freqBP4,inv).copy()
            stalta = trVelTrig.copy()
            dataraw = np.concatenate((dataraw,trRaw.data),axis=0)
            dataDisHighPass = np.concatenate((dataDisHighPass,trVelHighPass.data),axis=0)
            dataVelHighPass = np.concatenate((dataVelHighPass,trVelHighPass.data),axis=0)
            dataAccHighPass = np.concatenate((dataAccHighPass,trAccHighPass.data),axis=0)
            dataAccFilt = np.concatenate((dataAccFilt,trAccFilt.data),axis=0)
            dataVelFilt = np.concatenate((dataVelFilt,trVelFilt.data),axis=0)
            dataDisFilt = np.concatenate((dataDisFilt,trDisFilt.data),axis=0)
#            dataAccFilt10Hz = np.concatenate((dataAccFilt10Hz,trAccFilt10Hz.data),axis=0)
#            dataVelFilt10Hz = np.concatenate((dataVelFilt10Hz,trVelFilt10Hz.data),axis=0)
#            dataDisFilt10Hz = np.concatenate((dataDisFilt10Hz,trDisFilt10Hz.data),axis=0)
#            dataAccFilt20Hz = np.concatenate((dataAccFilt20Hz,trAccFilt20Hz.data),axis=0)
#            dataVelFilt20Hz = np.concatenate((dataVelFilt20Hz,trVelFilt20Hz.data),axis=0)
#            dataDisFilt20Hz = np.concatenate((dataDisFilt20Hz,trDisFilt20Hz.data),axis=0)
#            dataAccFilt25Hz = np.concatenate((dataAccFilt25Hz,trAccFilt25Hz.data),axis=0)
#            dataVelFilt25Hz = np.concatenate((dataVelFilt25Hz,trVelFilt25Hz.data),axis=0)
#            dataDisFilt25Hz = np.concatenate((dataDisFilt25Hz,trDisFilt25Hz.data),axis=0)

            if ( npts*dt > tbuffer ):
#                stalta = z_detect(trVelHighPass,int(0.05/dt))                      #--- of order 0.1 sec/trace 
                stalta = classic_sta_lta(trVelTrig,int(sta/dt),int(lta/dt))  #--- of order 0.4 sec/trace
                datastalta = np.concatenate((datastalta,stalta),axis=0)
            else:
                datastalta = np.concatenate((datastalta,np.zeros(trVelTrig.stats.npts)),axis=0)
#            if ( dt*npts > 3600.05 ):
#                powers = get_power(trRaw,inv,PSDperiods)

            T2 = timeit.default_timer()
#
#  Need to update the plots with the appropriate traces.  8/24/2018
#
            if ( iplot == 1 ):
                if ( j == i ):
                    strawplot = Stream(traces=[trRaw])
                    stVelHighPassplot = Stream(traces=[trVelHighPass])
                    stAccFiltplot = Stream(traces=[trAccFilt])
                    stVelFiltplot = Stream(traces=[trVelFilt])
                else:
                    strawplot += Stream(traces=[trRaw])
                    stVelHighPassplot += Stream(traces=[trVelHighPass])
                    stAccFiltplot += Stream(traces=[trAccFilt])
                    stVelFiltplot += Stream(traces=[trVelFilt])
                if ( j == i + ntr -1 ):
                    VelHighPasslabel = "cm/s " + ">" + str(freqHP) + " Hz"
                    AccFiltlabel = "cm/s^2 " + str(freqBP2) + "-" + str(freqBP3) + " Hz"
                    VelFiltlabel = "cm/s " + str(freqBP2) + "-" + str(freqBP3) + " Hz"
                    for ij in range(0,len(stVelHighPassplot)):
                        stVelHighPassplot[ij].data = stVelHighPassplot[ij].data*100
                        stAccFiltplot[ij].data = stAccFiltplot[ij].data*100
                        stVelFiltplot[ij].data = stVelFiltplot[ij].data*100
                    make_station_figure(strawplot,stVelHighPassplot,stAccFiltplot,stVelFiltplot,datastalta,"Raw",VelHighPasslabel,AccFiltlabel,VelFiltlabel,"Sta/Lta")  #--- of order 25sec/channel for 1 hr long traces
            T3 = timeit.default_timer()

#        pow50sec = powers[0]
#        pow30sec = powers[1]
#        pow20sec = powers[2]
#        pow12sec = powers[3]
#        pow10sec = powers[4]
        power_5sec = powers[5]
        power_1Hz = powers[6]
        power_5Hz = powers[7]
#        pow10Hz = powers[8]
#        pow20Hz = powers[9]
#        pow50Hz = powers[10]

        rawrange = float(np.ptp(dataraw))
        rawmean = np.mean(dataraw)
#        rawrms = np.sqrt(np.vdot(dataraw, dataraw)/dataraw.size)
        rawrms = np.sqrt(np.mean((dataraw-rawmean)**2))
        rawmin = min(dataraw)
        rawmax = max(dataraw)
        accmax = max(abs(dataAccFilt))*100.
        accmaxHP = max(abs(dataAccHighPass))*100.
        velmax = max(abs(dataVelFilt))*100.
        pctavailable = 100.*(nptstotal/(lenrequested/dt))
        ngaps = ntr - 1

        NoiseFloorAcc = noise_floor(dataAccFilt)*100.
        NoiseFloorAccHP = noise_floor(dataAccHighPass)*100.
        NoiseFloorVel = noise_floor(dataVelFilt)*100.

#        snr10_0p01cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,10,dt,twin,0.0001)
#        snr20_0p05cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.0005)  #--- of order 0.02 sec/trace each
#        snr20_0p082cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.00082)
#        snr20_0p1cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.001)
#        snr20_0p17cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.0017)  #--- 2017 ShakeAlert stat. acceptance thresh.
        snr20_0p34cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.0034)  #--- 2018 ShakeAlert stat. acceptance thresh.
#        snr20_1cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.01)
#        snr20_3cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.03)
#        snr20_5cm = count_peaks_stalta_new(dataAccFilt,datastalta,sta,lta,mpd,20,dt,twin,0.05)

        TrigTimes = []
        TrigTimesBB = []
        NTrigElarmSAlex,NboxcarsAlex,TrigTimes = count_peaks_stalta_Elarms_times(dataAccHighPass,dataVelHighPass,dataDisHighPass,datastalta,sta,lta,mpd,20,dt,twin,stAll[i].stats.channel) #--- Elarms threshold
        NTrigElarmSAlexBB,NboxcarsAlexBB,TrigTimesBB = count_peaks_stalta_Elarms_times(dataAccFilt,dataVelFilt,dataDisFilt,datastalta,sta,lta,mpd,20,dt,twin,stAll[i].stats.channel)
 
        for ijk in range(0,len(TrigTimes)):
            ishift = TrigTimes[ijk]
            trigtime = trVelHighPass.stats.starttime + (ishift*dt)
#            print(sncl,trigtime)
            outfile.write(sncl + " " + str(trigtime) + "\n" )
        for ijk in range(0,len(TrigTimesBB)):
            ishift = TrigTimesBB[ijk]
            trigtime = trVelFilt.stats.starttime + (ishift*dt)
            outfileBB.write(sncl + " " + str(trigtime) + "\n" )

#        NTrigElarmSAlexBB10Hz,NboxcarsAlexBB10Hz = count_peaks_stalta_Elarms(dataAccFilt10Hz,dataVelFilt10Hz,dataDisFilt10Hz,datastalta,sta,lta,mpd,20,dt,twin,stAll[i].stats.channel)
#        NTrigElarmSAlexBB20Hz,NboxcarsAlexBB20Hz = count_peaks_stalta_Elarms(dataAccFilt20Hz,dataVelFilt20Hz,dataDisFilt20Hz,datastalta,sta,lta,mpd,20,dt,twin,stAll[i].stats.channel)
#        NTrigElarmSAlexBB25Hz,NboxcarsAlexBB25Hz = count_peaks_stalta_Elarms(dataAccFilt25Hz,dataVelFilt25Hz,dataDisFilt25Hz,datastalta,sta,lta,mpd,20,dt,twin,stAll[i].stats.channel)

#        snr20_0p05cmBB10Hz = count_peaks_stalta_new(dataAccFilt10Hz,datastalta,sta,lta,mpd,20,dt,twin,0.0005)
#        snr20_0p05cmBB20Hz = count_peaks_stalta_new(dataAccFilt20Hz,datastalta,sta,lta,mpd,20,dt,twin,0.0005) 
#        snr20_0p05cmBB25Hz = count_peaks_stalta_new(dataAccFilt25Hz,datastalta,sta,lta,mpd,20,dt,twin,0.0005)

#        snr10_0p01cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,10,dt,twin,0.0001)
#        snr20_0p05cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.0005)   #--- of order 0.02 sec/trace each
#        snr20_0p082cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.00082) #--- Elarms threshold
#        snr20_0p1cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.001)
#        snr20_0p17cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.0017)   #--- 2017 ShakeAlert stat. acceptance thresh.
        snr20_0p34cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.0034)   #--- 2018 ShakeAlert stat. acceptance thresh.
#        snr20_1cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.01)
#        snr20_3cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.03)
#        snr20_5cmHP = count_peaks_stalta_new(dataAccHighPass,datastalta,sta,lta,mpd,20,dt,twin,0.05)

#        RMSduration_0p01cm = duration_exceed_RMS(dataAccFilt,0.0001,RMSlen,dt)  #--- of order 0.05 sec/trace each
#        RMSduration_0p035cm = duration_exceed_RMS(dataAccFilt,0.00035,RMSlen,dt)  #--- 2017 ShakeAlert stat. acceptance thresh.
        RMSduration_0p07cm = duration_exceed_RMS(dataAccFilt,0.0007,RMSlen,dt)  #--- 2018 ShakeAlert stat. acceptance thresh.
        RMSduration_0p07cmHP = duration_exceed_RMS(dataAccHighPass,0.0007,RMSlen,dt)  #--- 2018 ShakeAlert stat. acceptance thresh.
#        RMSduration_0p1cm = duration_exceed_RMS(dataAccFilt,0.001,RMSlen,dt)
#        RMSduration_1cm = duration_exceed_RMS(dataAccFilt,0.01,RMSlen,dt)

        finder_1cm_15Hz = count_triggers_FinDer(dataAccFilt,30., 0.01, dt)
        finder_1cm_hp = count_triggers_FinDer(dataAccHighPass,30., 0.01, dt)
        finder_2cm_15Hz = count_triggers_FinDer(dataAccFilt,30., 0.02, dt)
        finder_2cm_hp = count_triggers_FinDer(dataAccHighPass,30., 0.02, dt)

        #----- all metrics
#        numstring = "%d %d %d %d %d %d %d %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.6g %.6g %.4g %d    %.4g %.4g %.4g %.4g %.4g " % ( snr10_0p01cm, snr20_0p05cm, snr20_0p17cm, snr20_0p34cm, snr20_1cm, snr20_3cm, snr20_5cm, pow50sec, pow30sec, pow20sec, pow10sec, power_5sec, power_1Hz, power_5Hz, pow10Hz, pow20Hz, pow50Hz, NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps, RMSduration_0p01cm, RMSduration_0p035cm, RMSduration_0p07cm, RMSduration_0p1cm, RMSduration_1cm )

#        metriclist = [ "snr10_0p01cm", "snr20_0p05cm", "snr20_0p17cm", "snr20_0p34cm", "snr20_1cm", "snr20_3cm", "snr20_5cm", "pow50sec", "pow30sec", "pow20sec", "pow12sec", "pow10sec", "power_5sec", "power_1Hz", "power_5Hz", "pow10Hz", "pow20Hz", "pow50Hz", "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps", "RMSduration_0p01cm", "RMSduration_0p035cm", "RMSduration_0p07cm", "RMSduration_0p1cm", "RMSduration_1cm", "accmax", "velmax", "snr20_0p1cm", "snr20_0p082cm", "snr10_0p01cmHP", "snr20_0p05cmHP", "snr20_0p17cmHP", "snr20_0p34cmHP", "snr20_1cmHP", "snr20_3cmHP", "snr20_5cmHP", "snr20_0p082cmHP", "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB", "NTrigElarmSAlexBB20Hz", "NboxcarsAlexBB20Hz", "NTrigElarmSAlexBB25Hz", "NboxcarsAlexBB25Hz", "snr20_0p05cmBB20Hz", "snr20_0p05cmBB25Hz", "NTrigElarmSAlexBB10Hz", "NboxcarsAlexBB10Hz", "snr20_0p05cmBB10Hz", "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp" ]
#        valuelist = [ snr10_0p01cm, snr20_0p05cm, snr20_0p17cm, snr20_0p34cm, snr20_1cm, snr20_3cm, snr20_5cm, pow50sec, pow30sec, pow20sec, pow12sec, pow10sec, power_5sec, power_1Hz, power_5Hz, pow10Hz, pow20Hz, pow50Hz, NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps, RMSduration_0p01cm, RMSduration_0p035cm, RMSduration_0p07cm, RMSduration_0p1cm, RMSduration_1cm, accmax, velmax, snr20_0p1cm, snr20_0p082cm, snr10_0p01cmHP, snr20_0p05cmHP, snr20_0p17cmHP, snr20_0p34cmHP, snr20_1cmHP, snr20_3cmHP, snr20_5cmHP, snr20_0p082cmHP, NTrigElarmSAlex, NboxcarsAlex, NTrigElarmSAlexBB, NboxcarsAlexBB, NTrigElarmSAlexBB20Hz, NboxcarsAlexBB20Hz, NTrigElarmSAlexBB25Hz, NboxcarsAlexBB25Hz, snr20_0p05cmBB20Hz, snr20_0p05cmBB25Hz, NTrigElarmSAlexBB10Hz,  NboxcarsAlexBB10Hz ,  snr20_0p05cmBB10Hz, finder_1cm_15Hz, finder_1cm_hp, finder_2cm_15Hz, finder_2cm_hp ]

#        metrictuple = ( "snr10_0p01cm", "snr20_0p05cm", "snr20_0p17cm", "snr20_0p34cm", "snr20_1cm", "snr20_3cm", "snr20_5cm", "pow50sec", "pow30sec", "pow20sec", "pow12sec", "pow10sec", "power_5sec", "power_1Hz", "power_5Hz", "pow10Hz", "pow20Hz", "pow50Hz", "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps", "RMSduration_0p01cm", "RMSduration_0p035cm", "RMSduration_0p07cm", "RMSduration_0p1cm", "RMSduration_1cm", "accmax", "velmax", "snr20_0p1cm", "snr20_0p082cm", "snr10_0p01cmHP", "snr20_0p05cmHP", "snr20_0p17cmHP", "snr20_0p34cmHP", "snr20_1cmHP", "snr20_3cmHP", "snr20_5cmHP", "snr20_0p082cmHP", "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB", "NTrigElarmSAlexBB20Hz", "NboxcarsAlexBB20Hz", "NTrigElarmSAlexBB25Hz", "NboxcarsAlexBB25Hz", "snr20_0p05cmBB20Hz", "snr20_0p05cmBB25Hz", "NTrigElarmSAlexBB10Hz", "NboxcarsAlexBB10Hz", "snr20_0p05cmBB10Hz", "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp"  )

        #----- abbreviated metrics list
#        numstring = "%d %d %d %d %d %d %d %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.6g %.6g %.4g %d    %.4g %.4g %.4g %.4g %.4g " % ( snr10_0p01cm, snr20_0p05cm, snr20_0p17cm, snr20_0p34cm, snr20_1cm, snr20_3cm, snr20_5cm, power_5sec, power_1Hz, power_5Hz, NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps, RMSduration_0p01cm, RMSduration_0p035cm, RMSduration_0p07cm, RMSduration_0p1cm, RMSduration_1cm )

#        metriclist = [ "snr10_0p01cm", "snr20_0p05cm", "snr20_0p17cm", "snr20_0p34cm", "snr20_1cm", "snr20_3cm", "snr20_5cm", "power_5sec", "power_1Hz", "power_5Hz", "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps", "RMSduration_0p01cm", "RMSduration_0p035cm", "RMSduration_0p07cm", "RMSduration_0p1cm", "RMSduration_1cm", "accmax", "velmax", "snr20_0p1cm", "snr20_0p082cm", "snr10_0p01cmHP", "snr20_0p05cmHP", "snr20_0p17cmHP", "snr20_0p34cmHP", "snr20_1cmHP", "snr20_3cmHP", "snr20_5cmHP", "snr20_0p082cmHP", "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB", "NTrigElarmSAlexBB20Hz", "NboxcarsAlexBB20Hz", "NTrigElarmSAlexBB25Hz", "NboxcarsAlexBB25Hz", "snr20_0p05cmBB20Hz", "snr20_0p05cmBB25Hz", "NTrigElarmSAlexBB10Hz", "NboxcarsAlexBB10Hz", "snr20_0p05cmBB10Hz", "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp", "RMSduration_0p07cmHP" ]
#        valuelist = [ snr10_0p01cm, snr20_0p05cm, snr20_0p17cm, snr20_0p34cm, snr20_1cm, snr20_3cm, snr20_5cm, power_5sec, power_1Hz, power_5Hz, NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps, RMSduration_0p01cm, RMSduration_0p035cm, RMSduration_0p07cm, RMSduration_0p1cm, RMSduration_1cm, accmax, velmax, snr20_0p1cm, snr20_0p082cm, snr10_0p01cmHP, snr20_0p05cmHP, snr20_0p17cmHP, snr20_0p34cmHP, snr20_1cmHP, snr20_3cmHP, snr20_5cmHP, snr20_0p082cmHP, NTrigElarmSAlex, NboxcarsAlex, NTrigElarmSAlexBB, NboxcarsAlexBB, NTrigElarmSAlexBB20Hz, NboxcarsAlexBB20Hz, NTrigElarmSAlexBB25Hz, NboxcarsAlexBB25Hz, snr20_0p05cmBB20Hz, snr20_0p05cmBB25Hz, NTrigElarmSAlexBB10Hz,  NboxcarsAlexBB10Hz ,  snr20_0p05cmBB10Hz, finder_1cm_15Hz, finder_1cm_hp, finder_2cm_15Hz, finder_2cm_hp, RMSduration_0p07cmHP ]

#        metrictuple = ( "snr10_0p01cm", "snr20_0p05cm", "snr20_0p17cm", "snr20_0p34cm", "snr20_1cm", "snr20_3cm", "snr20_5cm", "power_5sec", "power_1Hz", "power_5Hz", "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps", "RMSduration_0p01cm", "RMSduration_0p035cm", "RMSduration_0p07cm", "RMSduration_0p1cm", "RMSduration_1cm", "accmax", "velmax", "snr20_0p1cm", "snr20_0p082cm", "snr10_0p01cmHP", "snr20_0p05cmHP", "snr20_0p17cmHP", "snr20_0p34cmHP", "snr20_1cmHP", "snr20_3cmHP", "snr20_5cmHP", "snr20_0p082cmHP", "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB", "NTrigElarmSAlexBB20Hz", "NboxcarsAlexBB20Hz", "NTrigElarmSAlexBB25Hz", "NboxcarsAlexBB25Hz", "snr20_0p05cmBB20Hz", "snr20_0p05cmBB25Hz", "NTrigElarmSAlexBB10Hz", "NboxcarsAlexBB10Hz", "snr20_0p05cmBB10Hz", "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp", "RMSduration_0p07cmHP"  )


        #----- very (squac) abbreviated metrics list
        numstring = "%d %d %d %d %d %d %d %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.4g %.6g %.6g %.4g %d    %.4g %.4g %.4g %.4g %.4g " % ( -1, -1, -1, snr20_0p34cm, -1, -1, -1, power_5sec, power_1Hz, power_5Hz, NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps, -1, -1, RMSduration_0p07cm, accmaxHP, NoiseFloorAccHP )

        metriclist = [  "snr20_0p34cm",  "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps",  "RMSduration_0p07cm",  "accmax", "velmax",  "snr20_0p34cmHP",  "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB", "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp", "RMSduration_0p07cmHP", "accmaxHP", "NoiseFloorAccHP" ]
        valuelist = [ snr20_0p34cm,  NoiseFloorAcc, NoiseFloorVel, rawrange, rawmean, rawrms, rawmin, rawmax, segmentshort, segmentlong, pctavailable, ngaps,  RMSduration_0p07cm, accmax, velmax,  snr20_0p34cmHP,  NTrigElarmSAlex, NboxcarsAlex, NTrigElarmSAlexBB, NboxcarsAlexBB,  finder_1cm_15Hz, finder_1cm_hp, finder_2cm_15Hz, finder_2cm_hp, RMSduration_0p07cmHP, accmaxHP, NoiseFloorAccHP ]

        metrictuple = (  "snr20_0p34cm", "NoiseFloorAcc", "NoiseFloorVel", "rawrange", "rawmean", "rawrms", "rawmin", "rawmax", "segmentshort", "segmentlong", "pctavailable", "ngaps",  "RMSduration_0p07cm",  "accmax", "velmax",  "snr20_0p34cmHP",  "NTrigElarmSAlex", "NboxcarsAlex", "NTrigElarmSAlexBB", "NboxcarsAlexBB",  "finder_1cm_15Hz", "finder_1cm_hp", "finder_2cm_15Hz", "finder_2cm_hp", "RMSduration_0p07cmHP", "accmaxHP", "NoiseFloorAccHP"  )


        T1 = timeit.default_timer()
        T2 = timeit.default_timer()

#        Tsum = Tsum + T2 - T1

        print( sncl + " " + numstring + "          ------" + str(itrace) + " " + str(T2-T1) + "  " + str(T2-T0) + " " + str(accmax) + " " + str(velmax) + " " + str(NTrigElarmSAlex) + " " + str(NboxcarsAlex) + " " + str(NTrigElarmSAlexBB) + " " + str(NboxcarsAlexBB) )

        '''
        if ( sncl in channel_map_staging ):
            nsquac_staging = nsquac_staging + 1
            upload_to_squac = [ rawmin, rawmax, rawmean, rawrange, finder_2cm_hp, RMSduration_0p07cm, RMSduration_0p07cmHP, snr20_0p34cm, snr20_0p34cmHP, NTrigElarmSAlex, NTrigElarmSAlexBB, NoiseFloorAcc, accmax, pctavailable, ngaps, segmentshort, segmentlong ]
            upload_to_squac_ids = [ 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98 ]
            for isquac in range(0,len(upload_to_squac)):
                value = upload_to_squac[isquac]
                metric_id = upload_to_squac_ids[isquac]

                measurement_staging = WriteOnlyMeasurementSerializer(
                              metric=metric_id,
                              channel=channel_map_staging[sncl],
                              value=value,
                              starttime=starttimesquac,
                              endtime=endtimesquac)
                measurements_staging.append(measurement_staging)
#            print(nsquac_staging,sncl,' ---- ',measurements)
        '''

        if ( sncl in channel_map_production ):
            nsquac_production = nsquac_production + 1
            upload_to_squac = [ rawmin, rawmax, rawmean, rawrange, finder_2cm_hp, RMSduration_0p07cm, RMSduration_0p07cmHP, snr20_0p34cm, snr20_0p34cmHP, NTrigElarmSAlex, NTrigElarmSAlexBB, NoiseFloorAcc, accmax, pctavailable, ngaps, segmentshort, segmentlong, accmaxHP, NoiseFloorAccHP ]
            upload_to_squac_ids = [ 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 109, 110 ]
            measurements_production = []
            for isquac in range(0,len(upload_to_squac)):
                value = upload_to_squac[isquac]
                metric_id = upload_to_squac_ids[isquac]

                measurement_production = WriteOnlyMeasurementSerializer(
                              metric=metric_id,
                              channel=channel_map_production[sncl],
                              value=value,
                              starttime=starttimesquac,
                              endtime=endtimesquac)
                measurements_production.append(measurement_production)
#            print(nsquac_production,sncl,' ---- ',measurements)
            T1 = timeit.default_timer()
  #####            response, errors = perform_bulk_create(measurements_production, squac_client_production, chunk=500)
            T2 = timeit.default_timer()
            Tsum = Tsum + (T2 - T1)
  #####            print("ERRORS: ", errors," N channels uploaded to SQUAC production: ",nsquac_production,"  in ",(T2-T1)," sec", Tsum, Tsum/nsquac_production )


#print('')
#print('measurements=')
#print(measurements)

#T1 = timeit.default_timer()
#response, errors = perform_bulk_create(measurements_staging, squac_client_staging)
#T2 = timeit.default_timer()
#print("ERRORS: ", errors," N channels uploaded to SQUAC staging: ",nsquac_staging,"  in ",(T2-T1)," sec")
##print(" N channels uploaded to SQUAC staging: ",nsquac_staging,"  in ",(T2-T1)," sec")

#T1 = timeit.default_timer()
#response, errors = perform_bulk_create(measurements_production, squac_client_production, chunk=500)
#T2 = timeit.default_timer()
#print("ERRORS: ", errors," N channels uploaded to SQUAC production: ",nsquac_production,"  in ",(T2-T1)," sec")
#print(" N channels uploaded to SQUAC production: ",nsquac_production,"  in ",(T2-T1)," sec")

outfile.close()
outfileBB.close()

