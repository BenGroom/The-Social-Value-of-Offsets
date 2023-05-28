Data and code for the paper The Social Value of Offsets, Nature 2023.

How calculate the SVO of your own project?
In excell file Value of an offset, set your parameters in sheet 1 (begin and end of project, risk, discount rate, ...)
Sheet 2 shows the value of the offset (cell N8) and the relative value expressed as a proportion of the SCC (cell O8). 
You can also value a project with a time-variable absorbed carbon stock (change column M) or a different temperature path (change column C).

Sheet 2, Column W calculates table SM5.
Sheet 3 shows the value of the project if the exact impact response function of Joos/Geoffroy (shown in Fig1) is used rather than an stepfunction.
Sheet 4 contains the calculations in the Table of SM1 (literature review).
The following sheets calculate the biomass example from SM9.

The Matlab file Table_SVO2.m calculates Table 1 (the excell sheet gives the same results). 

The Matlab file IRFproject.m makes Fig 1, Figure SM2 and Figures SM3a-e. 

The other files are input for the matlab files.

Ben and Frank

