# BLACKRAY
Raytracing code in non-Kerr spacetime with a convolver to obtain a full relativistic reflection spectrum using the XILLVER FITS file.

The current version calculates iron line profiles and full reflection spectra of infinitesimally-thin Novikov-Thorne disks in the Johannsen spacetime with deformation parameters $\alpha_{13}$, $\alpha_{22}$, $\alpha_{52}$, and $\epsilon_3$ (higher order deformation parameters are assumed to vanish). The Kerr solution is recovered for $\alpha_{13} = \alpha_{22} = \alpha_{52} = \epsilon_3 = 0$.

To run the code, we provide the Python script run.py. Before running the Python script, it is necessary to specify the location of the XILLVER table (line 9 in run.py) as well as all the parameters of the model (lines 10-23 in run.py). The parameters rstep and pstep regulate the resolution of the image of the observer and, in turn, the number of photons to fire to the disk. The default values are rstep = $1.008$ and pstep = $2\pi/720$, which should work for most situations. In the case you want to calculate spectra for the next generation of X-ray missions (e.g., eXTP or Athena), you may have to increase the resolution (e.g., rstep = $1.0001$ and pstep = $2\pi/3600$).

To run the Python script, the command line is

python3 run.py

run.py generates two outputs: an output for a relativistically broadened iron line profile and another output for a relativistically broadened full reflection spectrum.

NOTE: In the Johannsen spacetime, equatorial circular orbits are always vertically stable and the ISCO radius (which is assumed to set the inner edge of the accretion disk in BLACKRAY) is determined by the stability along the radial direction. In the current version of the code, the ISCO radius is thus calculated checking the stability along the radial direction only. However, there are non-Kerr spacetimes in which the ISCO radius is determined by the stability of the orbit along the vertical direction. If you change the metric and the new metric can have equatorial circular orbits that are vertically unstable, it is necessary to modify even the subroutine to determine the ISCO radius.
