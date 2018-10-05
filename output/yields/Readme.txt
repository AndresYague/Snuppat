The *.datYields files are the yield evolutions for the Snuppat simulations.
Each model is prefaced by a line describing the total mass and the total time
since the beginning of the TP-AGB phase, with the format:

"# Mass: [Total Mass] MSun. Time: [Time since beginning of TP-AGB] ky"

This header is then followed by a three column table with the format:
"[Element|Isotope name] [Proton number] [Yield in solar masses]"

For example, neutrons and nitrogen are both specified by the letter "n",
with the only discriminant between these in an "elements" file being the
proton number (0 and 7, respectively).

The initial_compAGS2009.dat file contains the initial composition used in the
simulations in mass number (mass fraction / isotope mass). This abundance is
then scaled to the selected metallicity for each specific model.

Additionally, SNUPPAT does not currently accurately predict HBB
nucleosynthesis. Therefore, for masses above 3.5 or 4 MSun, only heavier
than iron elemental yields should be trusted.
