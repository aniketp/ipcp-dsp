# Hardware assisted Cache Prefetching
Hardware prefetchers play an important role in hiding long data access latency 
in various multiprocessor systems. Significant number of prefetchers have been proposed using state-of-the-art techniques which independently augment 
the existing instruction pointer based and cache address based prefetching.

In this paper, we attempt to categorically evaluate the submissions of the _Third Data 
Prefetching Championship_ ([DPC3](http://dpc3.compas.cs.stonybrook.edu/)). Additionally, we also propose __IPCP++__, an enhanced version of Instruction Pointer Classifier based prefetcher (IPCP) reinforced with IP-Delta based sequence predictor and IP-based stride prefetcher (Sangam++). The confluence of these IP based prefetching techniques achieves a speedup of _9.39%_ over no prefetching when averaged over 20 single-thread 6XX SPEC CPU 2017 traces.

We also present a dynamic degree stream prefetcher. The dynamic stream prefetcher (DSP) takes into account factors like prefetch accuracy, cache pollution and strength of the stream in a given window of 1000 cycles to determine the prefetch degree dynamically.

This research work is done as a part of the Course Project for
Advanced Computer Architecture **(CS622A)**, Fall Semester - 2019, 
instructed by [Prof. Mainak Chaudhury](https://www.cse.iitk.ac.in/users/mainakc/).

The simulation results for all the prefetchers are compiled [here](./results).
And the project paper can be found [here](./Project-Paper.pdf).

### Group Members (G16)

| __Name__ | __Email__ | __Roll__ |
|-------------|------------|------------|
| Aditya Rohan | [raditya@iitk.ac.in](mailto:raditya@iitk.ac.in) | 160053 |
| Aniket Pandey | [aniketp@iitk.ac.in](mailto:aniketp@iitk.ac.in) | 160113 |


## Brief Build Instructions
* Install `champsim` repository
``` bash
 > cd scripts
 > ./install-champsim.sh
```

* Install DPC3 traces (SPEC 2017) <br>

Either install the traces from within champsim repository in the project or
move already existing traces to `champsim/dpc3_traces`.

* Build CPU executables and run simulation on all prefetchers for all available
traces.
``` golang
 > go run main.go
```

**NOTE** If you already have built the CPUs using above command and introduced
new traces for which you want to re-run the simulation, use
``` golang
 > go run main.go skip-cpu
```
This skips building the CPU executables again.

## Acknowledgement
* The Third Data Prefetching Championship ([online](https://dpc3.compas.cs.stonybrook.edu/))

* Sangam: A multi-component core cache prefetcher ([paper](https://dpc3.compas.cs.stonybrook.edu/pdfs/Sangam.pdf))

* Bouquet of instruction pointers: Instruction pointer classifier based hardware prefetching ([paper](https://dpc3.compas.cs.stonybrook.edu/pdfs/Bouquet.pdf))