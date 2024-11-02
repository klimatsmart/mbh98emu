# MBH98emu: a hockey stick emulation

MBH98emu is an accurate emulation of the iconic and much-debated Mann, Bradley & Hughes 1998 (MBH98) "hockey stick" temperature reconstruction.

## Status

In testing, the emulation of the Northern Hemisphere temperature reconstruction (the "hockey stick") is within a millionth of a degree of the original. The "old" version of the reconstruction, shown in Figure 7 of MBH98, is reproduced exactly (the archived data has been rounded to four decimal places). The accuracy might differ between computers since the calculations involve a lot of single-precision floating-point arithmetic. Regardless, this emulation should be much more accurate than those found in the literature. Reasons for this include corrections to the proxy networks and differences in how the proxy data is processed.

Principal component selection rules are outside the scope of this project, as the exact rules used in MBH98 are not known. Selections from MBH98 are hard-coded in `mbh98emu/config/`. This was handled similarly in published emulations.

## Usage

Download the MBH98emu repository and run the Python script `mbh98emu/script/mbh98emu.py`. The results are saved in `mbh98emu/reconstruction/` and `mbh98emu/validation/`. Some intermediate results are saved in other directories.

Corresponding files for the "old" version are located in `mbh98emu/old/`.

A run may take around 10 minutes. Use the `--fast` flag for faster but less accurate emulation.

## Background reading

Mann ME, Bradley RS, Hughes MK. 1998. Global-scale temperature patterns and climate forcing over the past six centuries. Nature. 392: 779–787

McIntyre S, McKitrick R. 2005. The M&M critique of the MBH98 Northern Hemisphere climate index: update and implications. Energy Environ. 16: 69–100

McIntyre S, McKitrick R. 2005. Hockey sticks, principal components, and spurious significance. Geophys Res Lett. 32: L03710

Wahl ER, Ammann C. 2007. Robustness of the Mann, Bradley, Hughes reconstruction of Northern Hemisphere surface temperatures: Examination of criticisms based on the nature and processing of proxy climate evidence. Clim Change. 85: 33–69
