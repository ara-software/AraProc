# AraProc

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