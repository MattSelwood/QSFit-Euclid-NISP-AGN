import QSFit: Options, add_qso_continuum!, add_patch_functs!, Job, JobState,
    EmLineComponent, SpecLineLorentz, SpecLineGauss, SpecLineVoigt

abstract type WP9Type2AGN <: DefaultRecipe end

function QSFit.Options(::Type{T}) where T <: WP9Type2AGN
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)
    out[:skip_lines] = Symbol[]
    out[:host_template] = Dict(:library=>"swire", :template=>"Ell5")
    out[:host_template_ref_wavelength] = 5500. # A
    out[:use_host_template] = true
    out[:host_template_range] = [4000., 7000.]
    out[:use_balmer] = false
    out[:use_ironuv] = false
    out[:use_ironopt] = false
    #out[:use_lorentzian_profiles] = false
    out[:n_unk] = 4
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50,
                       6565 .+ [-1,1] .* 150,
                       5008 .+ [-1,1] .* 25]  # Angstrom
    out[:unk_maxoffset_from_guess] = 1e3      # km/s

    out[:line_broadening] = true
    out[:norm_integrated] = true
    out[:line_profiles] = :voigt

    lines = OrderedDict{Symbol, QSFit.EmLineDescription}()
    lines[:Lya         ] = StdEmLine(:Lya       ,  :narrow)
    lines[:NIV_1483    ] = CustomEmLine(1483.32 ,  :narrow)
    lines[:NIV_1487    ] = CustomEmLine(1486.5  ,  :narrow)
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  ,  :narrow)
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :narrow)
    lines[:NeIV_2424   ] = CustomEmLine(2424.0  ,  :narrow) 
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :narrow)
    lines[:NeIII_3342  ] = CustomEmLine(3342.18 ,  :narrow)
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow)
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow)
    lines[:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow)
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow)
    lines[:Hd          ] = StdEmLine(:Hd        ,  :narrow)
    lines[:Hg          ] = StdEmLine(:Hg        ,  :narrow)
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow)
    lines[:Hb          ] = StdEmLine(:Hb        ,  :narrow)
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow)
    lines[:OIII_4959_bw] = StdEmLine(:OIII_4959 ,  :narrow)
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow)
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow)
    lines[:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow)
    lines[:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow)
    lines[:NII_6549    ] = StdEmLine(:NII_6549  ,  :narrow)
    lines[:Ha          ] = StdEmLine(:Ha        ,  :narrow)
    lines[:NII_6583    ] = StdEmLine(:NII_6583  ,  :narrow)
    lines[:SII_6716    ] = StdEmLine(:SII_6716  ,  :narrow)
    lines[:SII_6731    ] = StdEmLine(:SII_6731  ,  :narrow)
    lines[:ArIV_7171   ] = CustomEmLine(7170.7  ,  :narrow)
    lines[:OII_7325    ] = CustomEmLine(7325    ,  :narrow)  # Blend of 4 lines 
    lines[:Pa10        ] = CustomEmLine(9016    ,  :narrow)  # Are Paschen series narrow only in T2 ?
    lines[:SIII_9069   ] = CustomEmLine(9068.62 ,  :narrow) 
    lines[:Pa9         ] = CustomEmLine(9228.97 ,  :narrow)  
    lines[:SIII_9531   ] = CustomEmLine(9530.62 ,  :narrow)  
    lines[:Pa8         ] = CustomEmLine(9545.93 ,  :narrow)
    lines[:CI_9850     ] = CustomEmLine(9850    ,  :narrow)  # Blend of 3 lines 
    lines[:Pad         ] = CustomEmLine(10049   ,  :narrow)  # Paschen δ
    lines[:Pag         ] = CustomEmLine(10938   ,  :narrow)  # Paschen γ
    lines[:OI_11287    ] = CustomEmLine(11287   ,  :narrow)
    lines[:PII_11886   ] = CustomEmLine(11886   ,  :narrow)
    lines[:Pab         ] = StdEmLine(:Pab       ,  :narrow)  # Paschen β

    out[:lines] = lines
    return out
end


function QSFit.add_qso_continuum!(::Type{T}, job::JobState) where T <: WP9Type2AGN
    λ = coords(domain(job.model))

    comp = QSFit.powerlaw(median(λ))
    comp.alpha.val = -1.8

    job.model[:qso_cont] = comp
    push!(job.model[:Continuum].list, :qso_cont)
    GModelFit.update!(job.model)
end

function QSFit.EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Val{:narrow}) where T <: WP9Type2AGN
    lc = QSFit.EmLineComponent(supertype(T), job, λ, Val(:narrow)) # invoke parent recipe
    lc.comp.fwhm.low  = 10
    lc.comp.fwhm.high = 1000
    lc.comp.voff.high = 500
    return lc
end

function QSFit.EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Val{:unknown}) where T <: WP9Type2AGN
    lc = QSFit.EmLineComponent(supertype(T), job, λ, Val(:unknown)) # invoke parent recipe
    lc.comp.norm.val = 0.
    lc.comp.center.fixed = false
    lc.comp.center.low = 0
    lc.comp.center.high = Inf
    lc.comp.fwhm.val  = 500
    lc.comp.fwhm.low  = 10
    lc.comp.fwhm.high = 2e3
    lc.comp.voff.fixed = true
    return lc
end

function QSFit.add_patch_functs!(::Type{T}, job::JobState) where T <: WP9Type2AGN
    # Patch parameters
    if haskey(job.model, :OIII_4959)  &&  haskey(job.model, :OIII_5007)
        # job.model[:OIII_4959].norm.patch = @λ m -> m[:OIII_5007].norm / 3
        job.model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(job.model, :OIII_5007)  &&  haskey(job.model, :OIII_5007_bw)
        job.model[:OIII_5007_bw].voff.patch = @λ (m, v) -> v + m[:OIII_5007].voff
        job.model[:OIII_5007_bw].fwhm.patch = @λ (m, v) -> v + m[:OIII_5007].fwhm
        #job.model[:OIII_5007_bw].norm.patch = @λ (m, v) -> v + m[:OIII_5007].norm
    end
    if haskey(job.model, :OI_6300)  &&  haskey(job.model, :OI_6364)
        # job.model[:OI_6300].norm.patch = @λ m -> m[:OI_6364].norm / 3
        job.model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(job.model, :NII_6549)  &&  haskey(job.model, :NII_6583)
        # job.model[:NII_6549].norm.patch = @λ m -> m[:NII_6583].norm / 3
        job.model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(job.model, :SII_6716)  &&  haskey(job.model, :SII_6731)
        # job.model[:SII_6716].norm.patch = @λ m -> m[:SII_6731].norm / 3
        job.model[:SII_6716].voff.patch = :SII_6731
    end
end
