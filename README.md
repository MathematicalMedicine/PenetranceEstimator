betaApp
=======

This is a very simple app for calculating and displaying betaHat and betaHatStar given variables _alpha_, _beta_, _s_, and _k_.


Requirements
------------

* Python 3.7 or higher
* PyQt5
* matplotlib


Operation
---------

1. Start the program by running `python3 BetaHatApp`.
2. Select which model for calculating betaHat and betaHatStar is desired.
3. Enter in values for _alpha_, _beta_, and _s_ as desired. If multiple possible values (for plotting along the X-axis) are desired, this can be done by providing a comma-separated list. Only one variable may have multiple values at any given time.
4. Press the "Refresh Plot" button.

Statistics may be saved to a plain text file using the "Save statistics to file" button. There is also a self-explanatory "Save plot image to file" button.


Current Potential Gotchas
-------------------------

* You may specify a comma-separated list of values for _alpha_ OR _beta_, NOT both at the same time. Attempting to do both will crash the program.
* The backend code is capable of using multiple values for _s_ as well (hypothetically; it's untested), but the GUI does not presently allow for entering multiple values. This will likely be changed later.
* _k_ is not presently editable.
* _alpha_ and _beta_ values must be between 0 and 1. Attempts to use values outside of this range will not work properly.
