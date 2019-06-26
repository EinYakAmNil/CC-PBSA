# CC-PBSA
A fast tool to compute mutational free energy.

## Required software

+ Python 3.7:
  + python3-pymol
  + pandas
  + numpy
+ CONCOORD
+ GROMACS 2019.2
+ GroPBS (instructions in progress)

## Installation

The installation of some of the components is not trivial and differs depending on the OS.
The following sections will be about some problems that might occur.
The installation of GROMACS will not be covered here though, since it has a well [documented site](http://manual.gromacs.org/2019.2/install-guide/index.html) already.

### Linux

#### python3-pymol

This module does not seem to exist in in PyPI, so you have to get it from somewhere else.
On Ubuntu the software exists in the [apt repository](https://packages.ubuntu.com/disco/python3-pymol), so to install that use:

`sudo add-apt-repository deb http://de.archive.ubuntu.com/ubuntu eoan main universe`

`sudo apt update`

`sudo apt install python3-pymol`

No other software component should cause any trouble.

### Mac OS

#### Python

I found it easiest to just get [PyMOL](https://pymol.org/2/) and use its the own python interpreter.
Additional modules (like numpy or pandas) can then be install by just using `/path/to/PyMOL/bin/pip install module`.

#### CONCOORD

By itself [CONCOORD](https://www3.mpibpc.mpg.de/groups/de_groot/concoord/concoord.html) is basically just a .exe file, which should be able to run after download.
However it depends on a number of dynamic libraries with an absolute path.
To find out where it looks for them, run:

`otool -L /path/to/dist.exe`

If the specified path does not exist, but you have the library somewhere else,
use:

`install_name_tool -change /old/path /new/path /path/to/dist.exe`

If all the packages are installed. Run:

`python3 setup.py install --user`
