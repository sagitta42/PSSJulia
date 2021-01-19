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
daq_nsamples = 5000; # samples per waveform

daq_baseline_length = 1000;
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
gaussian_noise_dist = Normal(T(0), T(noise_σ)) #  Normal() -> Distributions.jjl

##

function read_mcpss(filename)
    @info "Reading mcpss"
    println(filename)

    HDF5.h5open(filename) do input
        Table(
            channel = readdata(input, "mcpss/raw/channel"),
            ievt = readdata(input, "mcpss/raw/ievt"),
            waveforms = readdata(input, "mcpss/raw/waveform")
        )
    end
end

function read_mctruth(filename)
    @info "Reading mcpss"
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
    daq_buffer_wv = RDWaveform(wf.time, daq_type.(round.(wf.value .* c .+ o, digits = 0)))
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
    daq_trigger_threshold = gaussian_noise_dist.σ * 10 * c

    online_filter_output = zeros(T, length(wf.value) - daq_trigger_window_length)
    t0_idx = 0

    for i in eachindex(online_filter_output)
        online_filter_output[i], trig = daq_online_filter(wf.value, i-1, daq_trigger_window_lengths, daq_trigger_threshold)
        if trig && t0_idx == 0
            t0_idx = i
        end
    end

    online_energy = uconvert(u"keV", maximum(online_filter_output) * germanium_ionization_energy / c)

    ts = range(T(0)u"ns", step = daq_Δt, length = daq_nsamples)
    iStart = t0_idx-daq_baseline_length
    stored_waveform = RDWaveform(ts, wf.value[iStart:iStart+daq_nsamples-1]);

    stored_waveform, online_energy
end

##

curdir = "/home/sagitta/_legend/pss/"

function main()
    ### Read mcpss events
    det_name = "V05266A"
    mcpss_file = "$curdir/cache/$(det_name)_mcpss.h5"
    mcpss = read_mcpss(mcpss_file)

    # select one RDwaveform for test
    wfpss = mcpss.waveforms[1]

    # plot_wf = plot(wf)
    # png(plot_wf, "step01-mcpss-decay.png")

    ### Differentiate
    wf = differentiate_wf(wfpss)

    # plot_wf = plot(wf)
    # png(plot_wf, "step02-mcpss-wf-current.png")

    ### 1. PreAmp Simulation

    #### 1.1. Simulate a charge sensitive amplifier (`CSA`)
    wf = simulate_csa(wf)

    # plot_wf = plot(
    #     begin
    #         plot(wfpss, label = "true")
    #         plot!(wf, label = "CSA Output")
    #     end,
    #     plot(wf, xlims = (1700-200, 1700+200)),
    #     layout = (2, 1)
    # )
    # png(plot_wf, "step03-preamp-decay1.png")


    #### 1.2. Noise
    wf = simulate_noise(wf)

    # plot_wf = plot(wf)
    # png(plot_wf, "step04-preamp-noise1.png")

    ### 2. DAQ Simulation

    # Now we want to simulate, what the DAQ records
    # The waveform will go rough the DAQ buffer on which the DAQ filter runs
    # in order to trigger in case of an event and calculate some parameters
    # like online energy
    # In order to do so. The waveform should be longer the than the DAQ Buffer (here, 5000 samples)
    # We took care of this at the beginning

    #### 2.1. DAQ units and baseline
    wf = daq_baseline(wf)

    # plot_wf = plot(wf)
    # png(plot_wf, "step05-daq-baseline1.png")

    #### 2.2. Trigger method
    wf, online_energy = daq_trigger(wf)

    # plot_wf = plot(wf, label = "Online Energy = $(online_energy)", title = "raw-data-like-waveform")
    # png(plot_wf, "step06-daq-trigger1.png")

    baseline, baseline_rms = mean_and_std(wf.value[1:daq_baseline_length])

    mctruth = read_mctruth(mcpss_file)

    t1pss_raw = Table(
        baseline = baseline,
        channel = mcpss.channel,
        energy = online_energy,
        ievt = mcpss.ievt,
        # numtraces = ?
        # packet_id = array of 1 (0?)
        # timestamp = mctruth.thit first element of each
        # tracelist = connected with numtraces
        # waveform = wf # but should be ArrayOfRDWaveforms
        # wf_max = ?
        # wf_min = ?
    )

    electruth = Table(
        # fill in noise and DAQ parameters
    )





end

##

main()
