#!/bin/csh

cd /home/seis/STATION_METRICS

set datacenter = 'NCEDC'
set list = "NCEDC.4"

#set chanfile = "ShakeAlertList.$list"
set chanfile = "chanfile.NCEDC.4"
set duration = 3600

set year = `date --date="2 hours ago" -u +%-Y`
set month = `date --date="2 hours ago" -u +%-m`
set day = `date --date="2 hours ago" -u +%-d`
set hour = `date --date="2 hours ago" -u +%-H`
set minute = 0
set second = 0
set sta = 0.05
set lta = 5.0
set padding = 120

set fBP1 = 0.06
set fBP2 = 0.075
set fBP3 = 15.0
set fBP4 = 17.0
set fHP = 3.0

set gainorfullresp = 'Gain'
#set gainorfullresp = 'Resp'
set mpd = 10.0
set rmslen = 5.0

#set iplot = 1
set iplot = 0

if ( -e status_running.$list ) then
   echo "Currently running: " $list $year.$month.$day.$hour
   echo "Skipped: " $list $year.$month.$day.$hour >> skipped_hours
else

set datenow = `date +"%Y.%m.%d.%H.%M.%S" -u`
set datenow2 = `date`
echo "Currently running: " $year.$month.$day.$hour $list $datenow $datenow2 > status_running.$list

if ( -e log.METRICS.$year.$month.$day.$hour.$list ) /bin/rm log.METRICS.$year.$month.$day.$hour.$list

/bin/cp calculate_station_metrics.py calculate_station_metrics.$list.py
sed -i 's/DATACENTER/'$datacenter'/g' calculate_station_metrics.$list.py
sed -i 's/CHANFILE/'$chanfile'/g' calculate_station_metrics.$list.py
sed -i 's/DURATION/'$duration'/g' calculate_station_metrics.$list.py
sed -i 's/YEAR/'$year'/g' calculate_station_metrics.$list.py
sed -i 's/MONTH/'$month'/g' calculate_station_metrics.$list.py
sed -i 's/DAY/'$day'/g' calculate_station_metrics.$list.py
sed -i 's/HOUR/'$hour'/g' calculate_station_metrics.$list.py
sed -i 's/MINUTE/'$minute'/g' calculate_station_metrics.$list.py
sed -i 's/SECOND/'$second'/g' calculate_station_metrics.$list.py
sed -i 's/STA/'$sta'/g' calculate_station_metrics.$list.py
sed -i 's/LTA/'$lta'/g' calculate_station_metrics.$list.py
sed -i 's/PADDING/'$padding'/g' calculate_station_metrics.$list.py
sed -i 's/fBP1/'$fBP1'/g' calculate_station_metrics.$list.py
sed -i 's/fBP2/'$fBP2'/g' calculate_station_metrics.$list.py
sed -i 's/fBP3/'$fBP3'/g' calculate_station_metrics.$list.py
sed -i 's/fBP4/'$fBP4'/g' calculate_station_metrics.$list.py
sed -i 's/fHP/'$fHP'/g' calculate_station_metrics.$list.py
sed -i 's/GAINORFULLRESP/'$gainorfullresp'/g' calculate_station_metrics.$list.py
sed -i 's/MPD/'$mpd'/g' calculate_station_metrics.$list.py
sed -i 's/RMSLEN/'$rmslen'/g' calculate_station_metrics.$list.py
sed -i 's/IPLOT/'$iplot'/g' calculate_station_metrics.$list.py

./calculate_station_metrics.$list.py > log.METRICS.$year.$month.$day.$hour.$list

/bin/rm status_running.$list
/bin/mv log.METRICS.$year.$month.$day.$hour.$list LOGFILES/

endif


