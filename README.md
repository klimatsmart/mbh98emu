# MBH98emu: a hockey stick emulation

MBH98emu is an accurate emulation of the iconic and much-debated Mann, Bradley & Hughes 1998 (MBH98) "hockey stick" temperature reconstruction.

## Status

In testing, the emulation of the Northern Hemisphere temperature reconstruction (the "hockey stick") is well within a millionth of a degree of the original. The small difference is likely due to roundoff error, and the exact numbers might differ between computers.

An "old" version of the reconstruction, shown in Figure 7 of MBH98, is emulated as well. The archived reconstruction is only given to four decimal places and is reproduced exactly.

As this project demonstrates, MBH98 used several tree-ring principal components the authors have claimed were excluded according to objective selection criteria. It is unclear if they were retained by mistake or if some undisclosed rule was applied. The emulation skips that part of the algorithm and reads the PC selections from configuration files instead.

## Usage

Run the Python script `mbh98emu/script/mbh98emu.py` to emulate the main MBH98 reconstruction. The results are saved in `mbh98emu/reconstruction/` and `mbh98emu/validation/`. Some intermediate results are saved in other directories.

Corresponding files for the "old" reconstruction are located in `mbh98emu/old/`.

A run may take around 10 minutes. Use the `--fast` flag for faster but less accurate emulation.

## Data sources

No instrumental or proxy data is hosted in this repository. Instead, the script downloads the necessary files `anomalies-new` and `mbh98.tar` to `mbh98emu/downloads/`.

This step can be carried out manually if the downloads fail. Links are given below.

### Instrumental data

The instrumental data is available from PSU and NOAA:

<https://www.meteo.psu.edu/holocene/public_html/shared/research/MANNETAL98/INSTRUMENTAL/anomalies-new><br>
<https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/mann1998/corrigendum2004/instrumental/anomalies-new.txt>

### Proxy data

The PSU and NOAA archives also contain proxy data, but this dataset was prepared after MBH98 was published and is missing several records the authors mistakenly believed they had not used.

Luckily, a copy of the original MBH98 proxy data was circulated online following a data breach at the Climatic Research Unit in 2009. The file is available here:

<https://sealevel.info/FOIA/2009/FOIA/documents/mbh98-osborn/mbh98.tar><br>
<http://junksciencearchive.com/FOIA/documents/mbh98-osborn/mbh98.tar>

## Background reading

Mann ME, Bradley RS, Hughes MK. 1998. Global-scale temperature patterns and climate forcing over the past six centuries. Nature. 392: 779–787

Mann ME, Bradley RS, Hughes MK. 2004. Correction: Corrigendum: Global-scale temperature patterns and climate forcing over the past six centuries. Nature. 430: 105

McIntyre S, McKitrick R. 2005. The M&M critique of the MBH98 Northern Hemisphere climate index: update and implications. Energy Environ. 16: 69–100

McIntyre S, McKitrick R. 2005. Hockey sticks, principal components, and spurious significance. Geophys Res Lett. 32: L03710

Wahl ER, Ammann C. 2007. Robustness of the Mann, Bradley, Hughes reconstruction of Northern Hemisphere surface temperatures: Examination of criticisms based on the nature and processing of proxy climate evidence. Clim Change. 85: 33–69
