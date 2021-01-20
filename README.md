# PSSJulia

## mcraw -> mcpss

```console
julia mcraw_to_mcpss.jl
```

Reads: `data/V05266A.json`, `cache/V05266A_mcraw.h5`

Writes: `cache/V05266A_mcpss.h5`

## mcpss -> t1pss

```console
julia mcpss_to_t1pss.jl
```

Doesn't do anything yet, but breaks down at around waveform number 68, giving an error "trying to access Array at negative index" at lin 194

```julia
stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq_nsamples-1]);
```

Line 201 contains a variable `curdir` because testing with Atom-Juno it would write into home directory otherwise and I don't know how to change that :)
