# PSSJulia

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

