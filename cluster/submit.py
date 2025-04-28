#! /usr/bin/env python3

import OctopusCluster.SimpleCluster as Cluster

jobids=Cluster.submitCommands(
     Cluster.Cluster.octopus_cloud2,
     [
        ['/mnt/rclsfserv005/users/tfq96423/dipole_fitting/run_test_application.sh',
        '/opt/matlabruntime/R2024b',
        '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/sims/sims_hinterer.tif',
        '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/sims/params_hinterer.csv',
        '/mnt/rclsfserv005/users/tfq96423/dipole_fitting/results/results_hinterer_test_11-20.csv',
        'hinterer',
        '988',
        '11',
        '20'],

     ],
singularity_image='/mnt/rclsfserv005/local/dipole_fitting.sif',
     delete_successful_logs=True,
     debug=True,
)

print('Submitted')
