import subprocess

cmd='/mnt/rclsfserv005/users/tfq96423/dipole_fitting/test_application /opt/matlabruntime/R2024b /mnt/rclsfserv005/users/tfq96423/dipole_fitting/sims/sims_hinterer.tif /mnt/rclsfserv005/users/tfq96423/dipole_fitting/sims/params_hinterer.csv /mnt/rclsfserv005/users/tfq96423/dipole_fitting/results/results_hinterer_test.csv hinterer 988'

subprocess.call(cmd, shell=True)
