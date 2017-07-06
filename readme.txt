To run:

[root_squared_error] = design_waveform(design_type, globals)

design_type = 0: no adaption
design_type = 1: MMSE waveform design
design_type = 2: RMB waveform design

globals is a structure with the fields, it must be called globals:

globals.ASNR (in dB)
globals.NXtheta (Ns - we used 250)
globals.L (forces it to max(NR,NT) if RMB)
globals.NR 
globals.NT 
globals.K (a multiple of 3 to plot properly)
globals.maxits (max iterations for cost function evaluation, we used 40)
globals.Ntargets (Number of targets, only works for 1 or 2)
globals.J (i.e., J in citation [4] - we used 10)

'.m' files to be made available upon first request after 24/09/2018

 

 




