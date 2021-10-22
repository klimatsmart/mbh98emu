# mbh98emu

## Description

The aim of this project is to accurately emulate the iconic "hockey stick" temperature reconstruction that appeared in the 1998 paper "Global-scale temperature patterns and climate forcing over the past six centuries" (MBH98) by Mann et al.

## Current status

Most of the algorithm has been implemented, the main exception being that the selections of empirical orthogonal functions for the calibration steps are hardcoded.

The emulation is within 0.00001 degrees of the original reconstruction.

## Usage

Start an Octave session, enter the Code directory, run the script full_algorithm.m and wait a few minutes.

## System requirements

The code has been tested on a computer with

OS:        Ubuntu 20.04  
Software:  Octave version 5.2  
Other:     IEEE 754 floating-point standard
