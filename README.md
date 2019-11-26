## Multi Source Self Calibration
Multi-Source Self-Calibration (MSSC) is a direction-dependent, calibration technique which provides an additional step to standard phase referencing. MSSC uses multiple faint sources detected within the primary beam and combines them together. The combined response of many sources across the field-of-view is more than sufficient to allow phase corrections to be derived. Each source have their models divided into the visibilities which results in multiple point sources. These are stacked in the uv plane to increase the S/N, permitting self-calibration to become feasible. It is worth noting that this process only applies to wide-field VLBI data sets that detect and image multiple sources within one epoch.  Recent improvements in the capabilities of VLBI correlators is ensuring that wide-field VLBI is a reality and as a result there will be an increased number of experiments which can utilise MSSC.

This new version is designed to be used with [CASA](https://casa.nrao.edu) and is under active development.

### Installation requirements
* [CASA](https://casa.nrao.edu)

### How to run
This is under development and will have an input file etc. For the moment you need to edit lines 20-25 in `MSSC_CASA.py` to set the inputs.

The code can then be run from within the CASA ipython environment using `execfile('MSSC_CASA.py')` or run directly in the command line via `casa -c MSSC_CASA.py`.

Until the new formulation is written into a paper, please reference [Radcliffe et al. (2016)](http://adsabs.harvard.edu/abs/2016A%26A...587A..85R) if you use this code. 
