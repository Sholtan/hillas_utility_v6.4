# hillas_utility
Utility for data processing TAIGA-IACT

First you need to compile the program using the Makefile. Just enter the 'make' in the Makefile directory

Then change the param_file according to the file you want to process. You need to select a telescope, a calibration file, pixel coordinates, paths to the source file, an out folder and a folder with the corresponding pointing files.

The launch is carried out from the console: ./iact_oper file_name_param_name

The program will display basic information about the run. Please note the delay between the start time and the selected pointing files. If everything is OK, press any key.

If several runs were written in the param_file, after the completion of the first one, the program will display information about the next run, you should again press any key.

It will contain a table with event parameters in csv format and text files with positions and pixel amplitudes of individual events.
