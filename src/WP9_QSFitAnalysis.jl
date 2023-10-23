module WP9_QSFitAnalysis

using Dierckx, DataStructures, Statistics, Dates
using MyAstroUtils, Gnuplot
using GModelFit, QSFit

export run_analysis, extract_catalogs

include("WP9Type1AGN.jl")
include("WP9Type2AGN.jl")
include("WP9Type1AGNSim.jl")
include("WP9Type2AGNSim.jl")


function getFilename(qsfit_source)
    fname = qsfit_source.name
    fname = replace(fname, " " => "_")
    fname = replace(fname, "(" => "")
    fname = replace(fname, ")" => "")
    return fname
end


function logGrid(x, y, e)
    xlog = exp10.(range(log10(minimum(x)), stop=log10(maximum(x)), length=length(x))) 
    ylog =      Spline1D(x, y)(xlog)
    elog = abs.(Spline1D(x, e)(xlog))
    return xlog, ylog, elog
end


function run_analysis(catalog_fits, spectra_fits, sim_dir, output)
    # Load catalog and spectra
    catalog = fits2df(catalog_fits)
    catalog = catalog[catalog.Z_OBS .> 0.9, :]

    f = FITS(spectra_fits)
    λ = read(f[2]) .* 1. # Read in wavelength array (same for all spectra)
    spectra = read(f[3]) # Read in spectra
    close(f)

    # Create output folders
    isdir(output)  ||  mkdir(output)

    # Loop through sources in catalog
    Threads.@threads for ii in 1:nrow(catalog)
        source_id = catalog[ii, :SOURCE_ID]
        z = catalog[ii, :Z_OBS]
        @info Threads.threadid() ii source_id z

        # Prepare incident spectrum
        y = 1.e+17 .* spectra[:, catalog[ii, :SED_INDEX] + 1]
        x, y, e = logGrid(λ, y, fill(1., length(y)))
        e = abs.(y .* 0.05)
        jj = findall(12000 .< x .< 19000)  # find wavelength samples in NISP red grism coverage
        x = x[jj]
        y = y[jj]
        e = e[jj]
        
        source = QSFit.Source("WP9 Incident $(source_id)", z)
        add_spec!(source, Spectrum(x, y, e))
        # Apply recipe depending on spectrum SOURCE_ID
        if source_id < 100000
            job = QSFit.Job{WP9Type1AGN}()
        else
            job = QSFit.Job{WP9Type2AGN}()
        end
        try
            filename = joinpath(output, getFilename(source))
            if !isfile(filename * ".json.gz")
                res = QSFit.run(source, job)
                # GFitViewer.save_html(ViewerData(source, res, comps=true), filename=filename * ".html")
                GModelFit.serialize(filename * ".json", res.bestfit, res.fitstats, res.pspec.data, compress=true)
            else
                @info "$filename already exists, skipping..."
            end
        catch err
            @warn "Error on source_id (incident): $source_id"
            println(err)
        end


        # Prepare simulated spectrum
        spec = fits2df(joinpath(sim_dir, string(source_id) * ".fits"))
        x, y, e = logGrid(spec.Wavelength .* 1., spec.Flux, spec.Uncertainty)
        jj = findall(12000 .< x .< 19000)  # find wavelength samples in NISP red grism coverage
        x = x[jj]
        y = y[jj]
        e = e[jj]
        source = QSFit.Source("WP9 Simulated $(source_id)", z)
        add_spec!(source, Spectrum(x, y, e, resolution=790))
        # Apply recipe depending on spectrum SOURCE_ID
        if source_id < 100000
            job = QSFit.Job{WP9Type1AGNSim}()
        else
            job = QSFit.Job{WP9Type2AGNSim}()
        end
        try
            filename = joinpath(output, getFilename(source))
            if !isfile(filename * ".json.gz")
                res = QSFit.run(source, job)
                # GFitViewer.save_html(ViewerData(source, res, comps=true), filename=filename * ".html")
                GModelFit.serialize(filename * ".json", res.bestfit, res.fitstats, compress=true)
            else
                @info "$filename already exists, skipping..."
            end
        catch err
            @warn "Error on source_id (simulated): $source_id"
            println(err)
        end
    end
end


function extract_catalogs(json_dir)
    # Read from zip file
    alldf = DataFrame()
    linenames = Vector{Symbol}()
    unknames  = Vector{Symbol}()

    for filename in readdir(json_dir, join=true)
        isnothing(match(r".*json\.gz", filename))  &&  continue
        @info "Reading $filename ..."
        res = GModelFit.deserialize(filename)

        if isnothing(findfirst("Simulated", filename))
            m = match(r".*WP9_Incident_(\d+)\.json\.gz", filename)
        else
            m = match(r".*WP9_Simulated_(\d+)\.json\.gz", filename)
        end
        df = DataFrame(:SOURCE_ID => Meta.parse(m.captures[1]))
        df[!, :sim] .= !isnothing(findfirst("Simulated", filename))
        df[!, :t2] .= (df[1, :SOURCE_ID] >= 100000)
        df[!, :file] .= filename

        for cname in keys(res[1])
            comp = res[1][cname]

            if !(cname in [:qso_cont, :balmer, :ironuv, :ironoptbr, :ironoptna, :galaxy, :Continuum, :main, :Iron, :BroadLines, :NarrowLines, :VeryBroadLines, :UnkLines])
                if isnothing(findfirst("unk", String(cname)))
                    push!(linenames, cname)
                    df[!, Symbol(cname, "_addcomp")] .= 0
                else
                    push!(unknames , cname)
                end
            end

            for (pname, par) in comp
                df[!, Symbol(cname, "_", pname)] .= par.actual
                df[!, Symbol(cname, "_", pname, "_err")] .= NaN
                if !par.fixed              &&
                    isnothing(par.patch)
                    df[!, Symbol(cname, "_", pname, "_err")] .= par.unc
                end
            end
        end

        df[!, :redchisq] .= res[2].fitstat
        df[!, :ndata]    .= res[2].ndata
        df[!, :dof]      .= res[2].dof
        append!(alldf, df, cols=:union)
    end

    linenames = unique(sort(linenames))
    unknames  = unique(sort(unknames))

    # Replace missing and nothing (which can not be stored in a FITS file)
    # with NaN, empty string or -9999 (depending on type)
    nn = names(alldf)
    for icol in 1:ncol(alldf)
        t = eltype(skipmissing(alldf[:, icol]))
        (t == Bool)  &&  continue
        i = findall(ismissing.(alldf[:, icol]))
        if t == Nothing
            nn[icol] = ""
        elseif t == String
            alldf[i, icol] .= ""
        elseif t == Int64
            alldf[i, icol] .= -9999
        else
            alldf[i, icol] .= NaN
        end

        if t == Union{Nothing, Float64}
            i = findall(isnothing.(alldf[:, icol]))
            if length(i) > 0
                alldf[i, icol] .= NaN
            end
            alldf[!, icol] .= float.(alldf[:, icol])
        end
    end
    select!(alldf, nn[nn .!= ""])
    disallowmissing!(alldf)

    #=
    Associate unknown lines
    Caveats:
    - only normalization is considered, FWHM and Voff are not modified
    - also normalization error is not considered since in a few cases normalization is patched, hence uncertainties are NaN
    =#
    for i in 1:nrow(alldf)
        for l in linenames
            for u in unknames
                c = alldf[i, Symbol(l, "_center")]
                w = alldf[i, Symbol(l, "_fwhm")]
                o = alldf[i, Symbol(l, "_voff")]
                if isfinite(c)  &&  isfinite(w)  && isfinite(o)
                    c *= (1 - o / 3e5)
                    w = w/3e5 * c
                    c2 = alldf[i, Symbol(u, "_center")]
                    w2 = alldf[i, Symbol(u, "_fwhm")]
                    w2 = w2/3e5 * c2
                    isfinite(w2)  &&  (w = max(w, w2))
                    if isfinite(c2)
                        if abs(c - c2) < w / 2
                            @assert isfinite(alldf[i, Symbol(l, "_norm")])
                            @assert isfinite(alldf[i, Symbol(u, "_norm")])
                            if  isfinite(alldf[i, Symbol(u, "_norm_err")])  &&
                                (alldf[i, Symbol(u, "_norm")] >
                                 alldf[i, Symbol(u, "_norm_err")])
                                @info "SOURCE_ID=$(alldf[i, :SOURCE_ID]): associating $l and $u"
                                alldf[i, Symbol(l, "_addcomp")]  += 1
                                alldf[i, Symbol(l, "_norm")]     += alldf[i, Symbol(u, "_norm")]
                                alldf[i, Symbol(l, "_norm_err")] += alldf[i, Symbol(u, "_norm_err")]
                            end
                        end
                    end
                end
            end
        end
    end


    # Sort column names
    nn = names(alldf)
    ii = findall((nn .== "file")     .|
                 (nn .== "SOURCE_ID").|
                 (nn .== "sim")      .|
                 (nn .== "t2")       .|
                 (nn .== "redchisq") .|
                 (nn .== "ndata")    .|
                 (nn .== "dof"))
    deleteat!(nn, ii)
    sort!(nn)
    nn = ["SOURCE_ID", "file", "sim", "t2", nn..., "redchisq", "ndata", "dof"]
    select!(alldf, nn)

    # Add continuum normalization at specific wavelengths
    for l in [1350, 3000, 5100]
        alldf[!, Symbol(:Cont, l)] .= alldf.qso_cont_norm .* (l ./ alldf.qso_cont_x0).^alldf.qso_cont_alpha
    end

    # Save FITS
    filenames = "QSFit_Catalog_" .* ["Incident", "Simulated"] .* "_" .* string(now())[1:end-7] .* ".fits"
    f = FITS(filenames[1], "w")
    write(f, alldf[findall(.!alldf.sim), :])
    close(f)
    f = FITS(filenames[2], "w")
    write(f, alldf[findall(  alldf.sim), :])
    close(f)

    println()
    println("Catalogs written in:")
    println(join(filenames, '\n'))
 end


end

