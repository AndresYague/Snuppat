SNUPPAT v1.0.0 Release

This is SNUPPAT, an S-process NUcleosynthesis Post-Processing code for ATon.
This code was written by Andrés Yagüe López (and.yague@gmail.com) in collaboration
with Dr. Paolo Venura, Dr. Aníbal García-Hernández and Dr. Maria Lugaro.

To compile the code a Makefile is included with the source, although it may be
compiled by hand. Also, a bash script was created to run the code containing a
rudimentary warning.

The "data" directory referred by some parts of the code is too big to upload
to GitHub. It contains the network and reaction lists. For further
information, please contact the author.

In order to run the code, an input file named "parameters.in" must be written.
The input file must have a specific layout. A template is detailed below.

-------------------TEMPLATE  BEGIN-----------------
1.D-5      # Relative accuracy
1.D-24     # Yscale
100        # ExtraNum
0.10       # Envelope overshooting parameter
0.002      # PDCZ overshooting parameter
advective  # overshooting mode (advective | diffusive)
0.018      # Metallicity
5          # MixFreq
5          # Minimum integDt
0          # Minimum pulsIntegDt
40         # Write frequency
--------------------TEMPLATE  END------------------

The meaning of each row is:

-Relative accuracy: The permitted error relative to the value.

-Yscale: The threshold dividing a meaningful quantity from zero.

-ExtraNum: Number of extra c13pocket shells.

-Envelope overshooting parameter: Overshooting parameter at the bottom of the
convective envelope.

-PDCZ overshooting parameter: Overshooting parameter at the bottom of the
PDCZ.

-Overshooting mode: Select between advective and diffusive overshooting

-Metallicity: Desired metallicity.

-MixFreq: Divide integration timestep in mixFreq steps to perform mixing (both
convective and extramixing).

-Minimum integDt: The minium dt in years to consider an integration. It may be
higher if the model or sum of models have inherently a greater timestep but
never lower.

-Minimum pulsIntegDt: Same as minimum integDt but for a pulse.

-Write frequency: The number of models to skip between calls to the output
writing function. Only "integrated" models count as skipped.

There are additional "Readme.txt" files in the "AtonOutput" and "input"
directories.
