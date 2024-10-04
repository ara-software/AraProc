# AraProc

## Installation

You are **strongly encouraged** to install AraProc into a 
python virtual environment built on an already existing ARA software stack.
The easiest way to do this is to grab the ARA software stack from cvmfs.
You can also build your own via the cvmfs scripts (https://github.com/ara-software/cvmfs).

You are also strongly encouraged to install it in "editable"
mode so that you can make changes to the analysis classes,
which you are very likely to do.

The instructions below show you the steps to the recommended installation method.

```

# activate the python in cvmfs
source /cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/setup.sh

# build a virtual environment
python3 -m venv /path/to/araproc_venv

# then activate your new virtual environment
source /path/to/araproc_venv/bin/activate

# then get, and install, AraProc
(araproc_venv) git clone git@github.com:clark2668/AraProc.git
(araproc_venv) cd AraProc
(araproc_venv) pip install -e .

```

To do work in the future, you will always need to start up your environment,
and then your virtual env. So you will always start by doing:
```
source /cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/setup.sh
source /path/to/araproc_venv/bin/activate
```
The former is necessary so that you get all the *other* ARA depenencies
(e.g. ROOT, GSL, etc.).
The latter actually turns on your venv.

### Some Warnings

One issue we encountered early were memory leaks.
In particular, AraRoot and FFTTools like to return *pointers* to ROOT objects.
E.g. [the AraRoot getGraphFromRFChan function](https://github.com/ara-software/AraRoot/blob/master/AraEvent/UsefulIcrrStationEvent.cxx#L51C1-L52C1).
This is dangerous, because pyroot doesn't know how to correctly cleanup
pointers; see for example, [this article](https://github.com/root-project/root/issues/11397)
about this on the ROOT forums.
The solution is that you need to explicitly give python control of the objects memory *when it is created*. For example:
```
wave = usefulPtr.getGraphFromRFChan(0)
ROOT.SetOwnership(wave, True)
```