# mbh98emu

## Description

An accurate emulation of the iconic "hockey stick" temperature reconstruction from the 1998 paper "Global-scale temperature patterns and climate forcing over the past six centuries" (MBH98) by Mann et al.

## Current status

Most of the algorithm has been implemented, the main exception being that the selections of empirical orthogonal functions for the calibration steps are hardcoded.

The emulation is within 0.00001 degrees of the original reconstruction.

## Usage

Start an Octave session, enter the Code directory, run the script full_algorithm.m and wait a few minutes. The results are saved in the directory Data/Emulation.

## System requirements

The code has been tested on a computer with
<pre>
OS:        Ubuntu 20.04
Software:  Octave version 5.2
Other:     IEEE 754 floating-point standard
</pre>
