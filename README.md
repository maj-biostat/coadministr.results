# coadministr.results

Runs coadministr simulation, has functionality to process the results and generate report.
Depends on [coadministr.stanc](https://github.com/maj-biostat/coadministr.stanc) and [coadministr](https://github.com/maj-biostat/coadministr) r packages.

## Getting started

### Running simulation

Edit `R/run_sim.R` as necessary, e.g. change number of simulations, scenarios etc, then

```
$ cd coadministr.results
$ Rscript R/run_sim.R > log.txt 2>&1 &
$ tail -f log.txt 
e15935c3 Started (scenario 1)
2021-04-08 16:25:57 rep 27
b3b11366 Started (scenario 1)
2021-04-08 16:25:58 rep 28
4fcb562d Started (scenario 1)
5b9847e0 Started (scenario 1)
...
```



