# PSSJulia

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

## Example output

### mcstp -> mcpss

```console
[redchuk@lfs1 PSSJulia]$ julia mcstp_to_mcpss.jl 
[ Info: Loading packages...
[ Info: Reading detector h5f
[ Info: I/O of charge drift model not yet supported. Loading default: ADLChargeDriftModel
Reading MC events from HDF5.
2767819 hits before clustering
 17.953064 seconds (39.31 M allocations: 4.119 GiB, 8.40% gc time)
294636 hits after clustering
[ Info: Group by detector...
[ Info: Adding fano noise...
[ Info: Simulating waveforms...
[ Info: Detector has 2 contact(s)
[ Info: Table has 800 physics events (2671 single charge depositions).
Progress: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:07
[ Info: Generating waveforms...
Extending tail -> 2000 baseline samples, wf length 8000
[ Info: Saving table...
Done
```

### mcpss -> mcraw

```console
[redchuk@lfs1 PSSJulia]$ julia mcpss_to_mcraw.jl 
[ Info: Loading packages...
[ Info: Reading mcpss
cache/V05266A_mcpss.h5
[ Info: Reading MC truth
cache/V05266A_mcpss.h5
[ Info: Saving table
Done
```
