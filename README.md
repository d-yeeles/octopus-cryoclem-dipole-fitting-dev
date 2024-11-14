diPOLE_python3_loop.py uses Mortensen's fixed dipole fitting code diPOLE.py, available from 'Supplementary Software' at the bottom of https://www.nature.com/articles/ncomms9621
I've translated it into python3 to work, and modified Estimate() so that it returns the position, orientation, and covariance matrix for each spot.
That code is intended for a single spot on a single frame, so here we just loop over every spot in a frame, and then over every frame.

One complication is that to do this we need to be able to identify spots in a frame. For this, we rely on thunderSTORM. diPOLE_python3_loop.py first runs thunderSTORM on a directory of frames. It then takes those localisations, extracts an NxN pixel region centred on them, and then applies Mortensen to that patch. This returns some adjustment which we add to the localisations given in the thunderSTORM results table.

Currently testing on:
/mnt/cryosil_ro/AttoDRY800/Developments_Autumn_2024/2024-10-09_ASIL240923C05/StormData/StormData_1/Results
and using the STORM parameters from the protocol.txt file in there.

ISSUES

1) Getting ImageJ to integrate properly with python is hard. I've tried to use pyimagej but it seems a bit dodgy, thought it might just be that I'm unfamiliar with ImageJ. It runs thunderSTORM fine, but I would prefer it to run headless if possible. I couldn't get that to work. I also tried to get it to run automatically at the end to generate the visualisation from the results table, but have been hitting a wall with that.

2) It's super slow. Like 15 seconds per spot / 5 minutes per frame. So this would never work for our usual 10,000 frame stack. The log-likelihood maximisation done in diPOLE_python3.py is done by fmin_powell(). That goes through about 200-300 function evaluations, each taking about 0.3 seconds, so ~5 seconds comes from that. Even with just a single evaluation, diPOLE_python3.py's Estimate() function still takes about 5 seconds, but I'm not sure which part is taking all that time.

3) I'm not using the right parameters for the experimental setup, or for the initial thunderSTORM run. So need to correct that and change this so it can read off metadata from the images or whatever. And Mortensen assumes EMCCD but I think we use CMOS?

4) There's a qick script to extract frames from the tif stack into a directory of frames. There's probably a better way to work with the tif stack directly.
