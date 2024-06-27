README for the MATLAB scripts
by Thomas Bissinger

Hi!
This is a collection of all the MATLAB files used to analyze data for the
2023 paper on the mobile XY model.

You can access all files here, though they are not commented and not all names may be
intelligible. As a personal recommendation, most files that could be useful to you are 
in the folders readfunctions and ScanningSnaphots, possibly also StructureFactor

- Many different plots can be created from the files in 2020_PaperPlots and 2023_ThesisPlots.
There are also some redundant versions present (they all just use basic MATLAB plotting routines,
and so they can be quite quickly recreated). You may browse through them when looking for 
a specific plot, however you may have a hard time finding an accurate plot. Note that there
are also many subfolders.

- auxiliary_funcs: These are functions implemented to handle the Nelson-Fisher law and some plotting
specifics. 

- FitFuncs: Routines to handle data fitting

- FourierTransform: Routines for Fourier transforms

- readfunctions: Routines for reading simulation data 
  This may be the most interesting folder for your purposes, the files in it can be used to gather
  data, extract information and collect it in .mat-Files. There are examples of files like that in
  the folder with the example calculation
  
- ScanningSnaphots: Files that can extract data from a snapshot.
  These may also be interesting to you. 

- StructureFactor: Can calculate structure factors and two point correlation functions
  Can be used for properties such as g(r), but also spin-spin correlations like <s(0) s(r)>.
  Also treats bond order parameter calculation (though possibly with implementation errors)
