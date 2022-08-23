#!/home/ahutko/miniconda3/envs/aws/bin/python

import os
import datetime
from pytz import timezone

# Import modules
try:
    from squacapi_client.models.write_only_measurement_serializer \
    import WriteOnlyMeasurementSerializer
    from squacapi_client.models.write_only_group_serializer \
    import WriteOnlyGroupSerializer
    from squacapi_client.pnsn_utilities \
    import get_client, make_channel_map, make_metric_map, perform_bulk_create
    no_squacapi = False
except Exception as e:
    print("Info: squacapi_client not available, cannot use --squac option")
    no_squacapi = True

USER = os.environ['SQUACAPI_USER']
PASSWORD = os.environ['SQUACAPI_PASSWD']
HOST = 'https://squacapi.pnsn.org'

# Iniate the client
squac_client = get_client(USER, PASSWORD, HOST )

# Example: retrieve all metrics
metrics = squac_client.v1_0_measurement_metrics_list()
metric_map = make_metric_map(metrics)
print(metric_map)

# Example: retrieve all channels.  Note can take up to 30s.
channels = squac_client.v1_0_nslc_channels_list(network='cc')
for channel in channels:
       print(channel)

