# PSSJulia

Recap of tiers

1. `g4simple`: Geant4 simulated events

2. `mcstp`: Geant4 events in the format requested by SSD / siggen

3. `mcpss`: output of SSD / siggen (waveforms)

4. `mcraw`: tier1-like format including DAQ and electronics effects    


## g4simple -> mcraw

Main script: `g4_to_mcraw.jl`

Run via `julia g4_to_mcraw.jl`

Launches the whole chain `g4 -> mcstp -> mcpss -> mcraw`.

Each intermediate stage may be saved if requested. Currently output of each stage is saved in `cache/`.

Inputs: detector geometry (example: `data/V05266A.json`), g4simple file (example: `/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/IC160A/Th228/uncollimated/top_source_holder/hdf5`)

The script calls three functions contained in the individual scripts `g4_to_mcstp.jl`, `mcstp_to_mcpss.jl`, `mcpss_to_mcraw.jl`. These scripts can also be launched separatly by uncommenting the `main()` function in each and running the directly with `julia code.jl`, defining the inputs in the code; or by commenting out stages in the main script.

## g4simple -> mcstp

```console
julia g4_to_mcstp.jl
```

Reads: g4simple hdf5 (here, a file in `/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/IC160A/Th228/uncollimated/top_source_holder/hdf5`)

Writes: `cache/filename_mcstp.h5` (filename same as input base name)    

## mcstp -> mcpss

```console
julia mcstp_to_mcpss.jl
```

Reads: `data/V05266A.json` (the first time, then reads cached `h5`), `cache/V05266A_mcstp.h5`

Writes: `cache/V05266A_mcpss.h5`, `cache/V05266A.h5f`

## mcpss -> mcraw

```console
julia mcpss_to_mcraw.jl
```

Reads: `cache/V05266A_mcpss.h5`

Writes: `cache/V05266A_mcraw.h5`

