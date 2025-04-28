#! /usr/bin/env python3

import OctopusCluster.SimpleCluster as Cluster

jobids=Cluster.submitCommands(
     Cluster.Cluster.octopus_cloud2,
     [
         ['ls', '-l'],
         ['id'],
     ],
singularity_image='/mnt/rclsfserv005/local/dipole_fitting.sif',
     delete_successful_logs=False,
     debug=True,
)
print(jobids)
