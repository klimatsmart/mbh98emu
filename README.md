# MBH98emu: a hockey stick emulation

MBH98emu is an accurate emulation of the iconic and much-debated Mann, Bradley & Hughes 1998 (MBH98) "hockey stick" temperature reconstruction.

## Status

In testing, the Northern Hemisphere temperature reconstruction (the "hockey stick") has small errors in the last (seventh) decimal place. Reconstructed temperature principal components also have small errors in the last (eighth) decimal place.

This accuracy depends on IEEE-754 compliance and may not apply to all systems. The emulation should still be orders of magnitude more accurate than emulations found in the literature, for various reasons.

The scope of MBH98emu does not include principal component selection rules, as the exact rules used in MBH98 are not known. Selections from MBH98 are hard-coded in `mbh98emu/config/`. This was handled similarly in published emulations.

## Usage

Download the MBH98emu repository and run the Python script `mbh98emu/script/mbh98emu.py`. The results are saved in `mbh98emu/reconstruction/` and `mbh98emu/validation/`. Some intermediate results are saved in other directories.

A run may take around 10 minutes. Use the `--fast` flag for faster but less accurate emulation.

## Background reading

Mann ME, Bradley RS, Hughes MK. 1998. Global-scale temperature patterns and climate forcing over the past six centuries. Nature. 392: 779–787

McIntyre S, McKitrick R. 2005. The M&M critique of the MBH98 Northern Hemisphere climate index: update and implications. Energy Environ. 16: 69–100

McIntyre S, McKitrick R. 2005. Hockey sticks, principal components, and spurious significance. Geophys Res Lett. 32: L03710

Wahl ER, Ammann C. 2007. Robustness of the Mann, Bradley, Hughes reconstruction of Northern Hemisphere surface temperatures: Examination of criticisms based on the nature and processing of proxy climate evidence. Clim Change. 85: 33–69
