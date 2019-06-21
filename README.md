# CC-PBSA

## Required software

+ Python 3.7:
  + python3-pymol
  + pandas
  + numpy
+ CONCOORD
+ GROMACS 2019.2
+ GroPBS (zip packages for Mac/Linux in the repository)

## Installation

The installation of some of the components is not trivial and differs depending on the OS.
The following sections will be about some problems that might occur.

### Linux

#### python3-pymol

This module does not seem to exist in in PyPI, so you have to get it from somewhere else.
On Ubuntu the software exists in the apt repository, so to install that just use:

`sudo apt install python3-pymol`

No other software component should cause any trouble.

### Mac OS

#### Python

I found it easiest to just get PyMOL and use its the own python interpreter.
Additional modules can then be install by just using `pip`.

#### CONCOORD

By itself CONCOORD is basically just a .exe file, which should be able to run after download.
However it depends on a number of dynamic libraries with an absolute path.
To find out where it looks for them, run:

`otool -L /path/to/dist.exe`
