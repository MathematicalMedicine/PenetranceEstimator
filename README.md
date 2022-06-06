Penetrance Estimator
====================

This is a very simple app for calculating and displaying the approximately ascertainment-corrected penetrance estimate (hereafter _f-tilde_) and the corresponding uncorrected penetrance estimate (hereafter _f-tilde-star_) given variables _alpha_, _f_, _s_, and _k_, and showing distributions of those calculations as a function of a number of simulated families. In essence, this serves as a pair of "interactive figures".

For details on what is going on here, please see the following paper:
> FIXME: paper not yet published


Requirements
------------

* Python 3.7 or higher
* PyQt5
* seaborn (and thus matplotlib)
* pandas


Operation
---------

1. Start the program by running `python3 PenEstApp`.
2. Enter in values for presented variables as desired.
3. Press the "Refresh Plot" button.

Statistics may be saved to a plain text file using the "Save statistics to file" button. There is also a self-explanatory "Save plot image to file" button.


Current Potential Gotchas
-------------------------

* _k_ for Figure 2 takes a comma-separated list of values. This is not presently properly validated, so incorrect input will crash the program.
* As _s_ increases (particularly as it goes above 5), refreshes slow down considerably because of the calculation time required. We intend to add a warning dialog to indicate this.
* The default (on launch) number of sim replicates for Figure 1 is set at 100. On older computers this may take a while to finish (minutes or more).


Acknowledgements
----------------

[Matthew Parker](https://github.com/mparker2) for the "sinaplot" code.
