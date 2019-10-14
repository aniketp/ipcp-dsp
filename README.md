# Data Prefetching Techniques
Course Project for Advanced Computer Architecture (CS622)

### Brief Build Instructions
* Install `champsim`
``` bash
 > cd scripts
 > ./install-champsim.sh
```

* Install DPC3 traces <br>

Either install the traces from within champsim directory in the project or
move already existing traces to `champsim/dpc3_traces`.

* Build CPU and Run simlution
``` golang
 > go run main.go
```

**NOTE** If you already have the CPUs built using above command and introduced
new traces for which you want to run the simulation, use
``` golang
 > go run main.go skip-cpu
```
This skips building the CPU executables again.