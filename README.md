# beam_alignment

This is a repository for some of the code used to analyze and correct for beam alignment problems in particle accelerators. 
All these files are meant to be compiled using ROOT v5

<====================================================================>

The file:
eta_modification_grid_original.C

is used to calculate the modification of detectors due to beam alignment effects

The file:
/boost_angle/simulate_boost_rotation_distortion.C

is used to simulate the distortion of the emitted particle angle distributions due to beam alignment effects.

The file:
/boost_angle/scan_boost_angle.C

is used to scan over multiple possible beam orientations and computes the Lorentz boost and rotation necessary in order to bring the beams to a colinear reference frame

The folder:
/boost_angle/derivation/

is the latex files to produce the pdf detailing the derivation of beam orientaion boost and rotation dependence.