matlab_abp
==========

MATLAB tools for feature extraction from Arterial Blood Pressure waveforms

This repo contains the following third-party dependencies (under the external
 directory):

* [Cardiac Output Estimation from Arterial Blood Pressure Waveforms](http://www.physionet.org/physiotools/cardiac-output/). 
Licensed under the terms of the [GNU General Public License (GPL)](http://www.fsf.org/copyleft/gpl.html)

`matlab_abp` depends on these other repos:
 
 * [matlab_misc](https://github.com/germangh/matlab_misc)
 * [matlab_mperl](https://github.com/germangh/matlab_misc)
 * [matlab_io](https://github.com/germangh/matlab_io)
 
To install, download `matlab_abp` and the three repos above. Then add them 
(including directories) to your MATLAB path and run the following to get 
some usage information:

    import abp.*;
	help abp_features

