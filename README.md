# MBH98emu

## Description

MBH98emu is an accurate emulation of the Mann, Bradley & Hughes 1998 (MBH98) "hockey stick" temperature reconstruction.

## Current status

The emulated Northern Hemisphere reconstruction should be within 0.00001 °C of the original on most modern systems. That is several orders of magnitude closer than published emulations.

## Usage

Run the Python script mbh98emu/script/mbh98emu.py. The results are saved in mbh98emu/reconstruction/ and mbh98emu/validation/.

Tested on Linux with Python 3.9.2 and the following packages installed:
<pre>
requests  2.25.1
NumPy     1.19.5
SciPy     1.6.0
pandas    1.1.5
</pre>

## Background reading

Mann ME, Bradley RS, Hughes MK. 1998. Global-scale temperature patterns and climate forcing over the past six centuries. Nature. 392: 779–787

McIntyre S, McKitrick R. 2005. The M&M critique of the MBH98 Northern Hemisphere climate index: update and implications. Energy Environ. 16: 69–100

McIntyre S, McKitrick R. 2005. Hockey sticks, principal components, and spurious significance. Geophys Res Lett. 32: L03710

Wahl ER, Ammann C. 2007. Robustness of the Mann, Bradley, Hughes reconstruction of Northern Hemisphere surface temperatures: Examination of criticisms based on the nature and processing of proxy climate evidence. Clim Change. 85: 33–69
