# Convert mcpss HDF5 format to t1pss
# 19.01.2021 Mariia Redchuk mariia.redchuk@pd.infn.it
# Author: Lukas Hauertmann

@info "Loading packages..."

using DSP # Julia's Digital Signal Processor (DSP) Package
using RadiationDetectorSignals # #master branch - By Oliver Schulz
# using RadiationDetectorDSP # By Oliver Schulz
# Most of this notebook should be part of RadiationDetectorDSP.jl.
# We will have to work on this....

using Unitful
using Distributions
using StatsBase # Histograms
using Random

using ArraysOfArrays, Tables, TypedTables, Unitful
using RadiationDetectorSignals
using LegendDataTypes, LegendHDF5IO
using LegendHDF5IO: readdata, writedata
import HDF5

using Plots

## Fix basic parameters

T = Float32

germanium_ionization_energy = T(2.95)u"eV"

# DAQ

# Let's say DAQ always stores 5000 samples. So the final waveform length is 5000.
# Has to be smaller than total waveform length (total_waveform_length in mcraw_to_mcpss.jl)
daq_nsamples = 4000; # samples per waveform

daq_baseline_length = 400;
daq_sampling_rate = 250u"MHz"
daq_Δt = uconvert(u"ns", inv(daq_sampling_rate));

daq_type = UInt16
daq_max_e = 10_000u"keV" # == typemax(UInt16)
daq_offset = 2_000u"keV";
c = typemax(daq_type) / ((daq_max_e-daq_offset)/uconvert(u"keV", germanium_ionization_energy))

# PreAmp (pa)
# pa_τ_decay = T(50)u"μs"; # decay constant of the preamp
# pa_τ_rise = T(20)u"ns"; # has something to do with the bandwidth of the preamp I believe
pa_τ_decay = T(5)*u"μs"
pa_τ_rise = T(2)*u"ns"

noise_σ = uconvert(u"eV", T(3)u"keV") / germanium_ionization_energy


##

function read_mcpss(filename)
    @info "Reading mcpss"
    println(filename)

    HDF5.h5open(filename) do input
        Table(
            channel = readdata(input, "mcpss/mcpss/channel"),
            ievt = readdata(input, "mcpss/mcpss/ievt"),
            waveforms = readdata(input, "mcpss/mcpss/waveform")
        )
    end
end

function read_mctruth(filename)
    @info "Reading MC truth"
    println(filename)

    HDF5.h5open(filename) do input
        Table(
            detno = readdata(input, "mcpss/mctruth/detno"),
            edep = readdata(input, "mcpss/mctruth/edep"),
            ievt = readdata(input, "mcpss/mctruth/ievt"),
            pos = readdata(input, "mcpss/mctruth/pos"),
            thit = readdata(input, "mcpss/mctruth/thit")
        )
    end
end
##

function dspjl_differentiator_filter(gain::Real)
    # The waveform is in units of induced charge (integrated form).
    # But we want to transform it into its differential form (current).
    # We will use a BiQuad filter for this: https://en.wikipedia.org/wiki/Digital_biquad_filter
    # The following functions to create filters are in general in `RadiationDetectorDSP` .jl(#dev branch)
    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end

function differentiate_wf(wf)
    diff_biquad_filter = dspjl_differentiator_filter(T(1)) # We want a gain of 1 here
    filter_output = filt(diff_biquad_filter, wf.value)
    # !!! We have to implement a method do parse a RDWaveform to `filt`
    wf_diff = RDWaveform(wf.time, filter_output)
    wf_diff
end


##

# simulate a charge sensitive amplifier (`CSA`)
# For this, we need a `RC` filter and an `integrator` filter (the inverser of the `dspjl_differentiator_filter`)
# This filters can all be made by BiQuad filter (with different paramters)
# These BiQuad filters can also be combined together nicely, see `dspjl_simple_csa_response_filter` where to filter are multiplied.

function dspjl_rc_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(α), T(0), T(0), T(α - 1), T(0))
end

function dspjl_integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end

function dspjl_simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter
    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    dspjl_rc_filter(T(τ_rise)) * dspjl_integrator_cr_filter(T(gain), T(τ_decay))
end

function simulate_csa(wf)
    # Here, the parameters τ_rise and τ_decay have to be given in units of samples,
    # because the `filt`-function does not know the Δt between the samples.

    csa_filter = dspjl_simple_csa_response_filter(
        pa_τ_rise / step(wf.time),
        uconvert(u"ns", pa_τ_decay) / step(wf.time))

    pa_wf = RDWaveform(wf.time, filt(csa_filter, wf.value))
    pa_wf
end

##

function simulate_noise(wf)
    # I am no expert here. I don't know at which point one should introduce noise.
    # Also, different noise could be added at different stages. This really depends on the electronics.
    # I will just add some Gaussian Noise (σ of 3 keV defined on top)
    # lets generate 1000 random samples from this normal distribution
    gaussian_noise_dist = Normal(T(0), T(noise_σ)) #  Normal() -> Distributions.jjl
    samples = rand(gaussian_noise_dist, 1000)
    # h = fit(Histogram, samples, nbins = 50) # -> StatsBase.jl
    # plot(h)

    # Now, lets add this Gaussian noise to other waveform (here, after the filters (but might be also added before))
    wf_noise = RDWaveform(wf.time, wf.value .+ rand!(gaussian_noise_dist, similar(wf.value)))
    wf_noise
end

##

function daq_baseline(wf)
    o = c * uconvert(u"eV", daq_offset) / germanium_ionization_energy

    # invert the pulse if needed
    sign = wf.value[end] < 0 ? -1 : 1

    daq_buffer_wv = RDWaveform(wf.time, daq_type.(round.(sign * wf.value .* c .+ o, digits = 0)))
    daq_buffer_wv
end

function daq_online_filter(values::AbstractVector, offset::Int, window_lengths::NTuple{3, Int}, threshold)
    wl = window_lengths
    r1 = offset+wl[1]+wl[2]:offset+wl[1]+wl[2]+wl[3]
    r2 = offset+wl[1]+wl[2]+wl[3]:offset+sum(window_lengths)
    r = mean(values[r2]) - mean(values[r1])
    r, r >= threshold
end

function daq_trigger(wf)
    daq_trigger_window_lengths = (250,250,250)
    daq_trigger_window_length = sum(daq_trigger_window_lengths)
    daq_trigger_threshold = noise_σ * 10 * c

    online_filter_output = zeros(T, length(wf.value) - daq_trigger_window_length)
    t0_idx = 0
    trig = false

    # while(not trig)
    #     online_filter_output[i], trig = daq_online_filter(wf.value, i-1, daq_trigger_window_lengths, daq_trigger_threshold)
    #     t0_idx = i

    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = daq_online_filter(wf.value, i-1, daq_trigger_window_lengths, daq_trigger_threshold)
        if trig && t0_idx == 0
            t0_idx = i
        end
    end

    ts = range(T(0)u"ns", step = daq_Δt, length = daq_nsamples)
    # in case it didn't trigger
    if(t0_idx == 0)
        stored_waveform = RDWaveform(ts, wf.value[1:daq_nsamples]) # just to return something
        online_energy = 0u"keV" # flag meaning didn't trigger
    else
        online_energy = uconvert(u"keV", maximum(online_filter_output) * germanium_ionization_energy / c)
        iStart = t0_idx-daq_baseline_length
        stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq_nsamples-1]);
    end

    stored_waveform, online_energy
end

##

function main()
    ### Read mcpss events
    det_name = "V05266A"
    mcpss_file = "cache/$(det_name)_mcpss.h5"
    mcpss = read_mcpss(mcpss_file)
    
    ### Create arrays to be filled with results, and online energy
    idx_end = size(mcpss.waveforms,1)
    wf_array = Array{RDWaveform}(undef, idx_end)
    temp = 1.0u"keV"; K = typeof(temp);
    online_energy = Array{K}(undef, idx_end)
    baseline = Array{T}(undef, idx_end)
    baseline_rms = Array{T}(undef, idx_end)
   
    @info "Processing waveforms..."
    ### loop over each wf and process it
    for i = 1:idx_end
        ### Differentiate
        wf = differentiate_wf(mcpss.waveforms[i])
        wf_array[i] = wf
    
        ### 1. PreAmp Simulation
    
        #### 1.1. Simulate a charge sensitive amplifier (`CSA`)
        wf_array[i] = simulate_csa(wf_array[i])
    
        # plot_wf = plot(
        #     begin
        #         plot(mcpss.waveforms[i], label = "true")
        #         plot!(wf_array[i], label = "CSA Output")
        #     end,
        #     plot(wf_array[i], xlims = (2000-1000, 2000+1000)),
        #     layout = (2, 1)
        # )
        # png(plot_wf, "step03-preamp-decay1_wf$i.png")
    
    
        #### 1.2. Noise
        wf_array[i] = simulate_noise(wf_array[i])
    
    
        ### 2. DAQ Simulation
    
        # Now we want to simulate, what the DAQ records
        # The waveform will go rough the DAQ buffer on which the DAQ filter runs
        # in order to trigger in case of an event and calculate some parameters
        # like online energy
        # In order to do so. The waveform should be longer the than the DAQ Buffer (here, 5000 samples)
        # We took care of this at the beginning
    
        #### 2.1. DAQ units and baseline
        wf_array[i] = daq_baseline(wf_array[i])
    
    
        #### 2.2. Trigger method
        # if online energy is zero, means didn't trigger
        wf_array[i], online_energy[i] = daq_trigger(wf_array[i])
    
        # plot_wf = plot(wf_array[i], label = "Online Energy = $(online_energy[i])", title = "raw-data-like-waveform")
        # png(plot_wf, "step06-daq-trigger1_wf$i.png")
    
        baseline[i], baseline_rms[i] = mean_and_std(wf_array[i].value[1:daq_baseline_length])
    
    end
    
    wf_final = ArrayOfRDWaveforms(wf_array)
    wf_final = ArrayOfRDWaveforms((wf_final.time, VectorOfSimilarVectors(wf_final.value)))
    
    # plot_wf = plot(wf_array)
    # png(plot_wf, "step02-mcpss-wf-current.png")
    # png(plot_wf, "step04-preamp-noise1.png")
    # png(plot_wf, "step05-daq-baseline1.png")
    # println(online_energy[5])
    # println(baseline[5])
    
    mctruth = read_mctruth(mcpss_file)
    
    ##
    mcraw = Table(
        baseline = baseline,
        channel = mcpss.channel,
        energy = online_energy,
        ievt = mcpss.ievt,
        numtraces = ones(length(baseline)), # number of triggered detectors (1 for HADES)
        packet_id = zeros(length(baseline)), # means to packet losses
        timestamp = getindex.(mctruth.thit, 1), # frist MC truth hit time of each event?
        tracelist = VectorOfVectors([[1] for idx in 1:length(baseline)]), # lists of ADCs that triggered, 1 for HADES all the time
        waveform = wf_final,
        wf_max = maximum.(wf_final.value), # ?
        wf_std = std.(wf_final.value[1:daq_baseline_length]) # ?
    )
    
    
    ##
    
    #
    # electruth = Table(
    #     # fill in noise and DAQ parameters
    # )
    
    #
    
    @info "Saving table..."
    out_filename = "cache/$(det_name)_mcraw.h5"
    println("-> $out_filename")
    HDF5.h5open(out_filename, "w") do f
        LegendDataTypes.writedata(f, "raw", mcraw)
    end
    println("Done")


end

##

main()
