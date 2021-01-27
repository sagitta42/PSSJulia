include("g4_to_mcstp.jl")
include("mcstp_to_mcpss.jl")
include("mcpss_to_mcraw.jl")

function main()

    det_name = "V05266A"
    g4_dir = "/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/IC160A/Th228/uncollimated/top_source_holder/hdf5/"
    mc_filename = "raw-IC160A-Th228-uncollimated-top-run0002-source_holder-bi-hdf5-01-test"
    processed_dir = "cache/"

    @info "----- g4simple -> mcstp"
    mcstp = g4_to_mcstp(g4_dir, processed_dir, mc_filename, ".hdf5", ".h5", save=false)

    @info "----- mcstp -> mcpss"
    mc_events = prepare_mcstp(mcstp)
    mcpss, mctruth = mcstp_to_mcpss(det_name, mc_events)
#    mcpss, mctruth = mcstp_to_mcpss(det_name, mc_events, mc_name=mc_filename) # if name given, will save the file

    @info "----- mcpss -> mcraw"
    mcpss_to_mcraw(mcpss, mctruth, mc_filename) # saves the final output
end


main()
