import QSFit: Options, add_qso_continuum!, add_patch_functs!, Job, JobState,
    EmLineComponent, SpecLineLorentz, SpecLineGauss, SpecLineVoigt

abstract type WP9Type1AGN <: DefaultRecipe end

function QSFit.Options(::Type{T}) where T <: WP9Type1AGN
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)
    out[:skip_lines] = Symbol[]

    out[:use_host_template] = false
    out[:host_template_ref_wavelength] = 5500. # A
    out[:use_balmer] = true
    out[:use_ironuv] = true;      out[:ironuv_fwhm]    = 3000.
    out[:use_ironopt] = true;     out[:ironoptbr_fwhm] = 3000.;  out[:ironoptna_fwhm] =  500.
    out[:use_lorentzian_profiles] = false
    out[:n_unk] = 5
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50,
                       6565 .+ [-1,1] .* 150,
                       5008 .+ [-1,1] .* 25] # Angstrom
    out[:unk_maxoffset_from_guess] = 1e3     # km/s

    out[:line_broadening] = true
    out[:iron_broadening] = true
    out[:norm_integrated] = true
    out[:line_profiles] = :gauss

    out[:lines] = OrderedDict{Symbol, QSFit.EmLineDescription}()
    out[:lines][:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :narrow, :broad )
    out[:lines][:CII_2326    ] = StdEmLine(:CII_2326  ,  :broad  )
    out[:lines][:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :narrow, :broad)
    out[:lines][:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow )
    out[:lines][:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow )
    out[:lines][:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow )
    out[:lines][:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow )
    out[:lines][:Hd          ] = StdEmLine(:Hd        ,  :broad  )
    out[:lines][:Hg          ] = StdEmLine(:Hg        ,  :broad  )
    out[:lines][:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow )
    out[:lines][:HeII_4686   ] = StdEmLine(:HeII_4686 ,  :broad  )
    out[:lines][:Hb          ] = StdEmLine(:Hb        ,  :narrow, :broad)
    out[:lines][:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow )
    out[:lines][:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow )
    out[:lines][:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow )
    out[:lines][:HeI_5876    ] = StdEmLine(:HeI_5876  ,  :broad  )
    out[:lines][:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow )
    out[:lines][:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow )
    out[:lines][:NII_6549    ] = StdEmLine(:NII_6549  ,  :narrow )
    out[:lines][:Ha          ] = StdEmLine(:Ha        ,  :narrow, :broad, :verybroad)
    out[:lines][:NII_6583    ] = StdEmLine(:NII_6583  ,  :narrow )
    out[:lines][:SII_6716    ] = StdEmLine(:SII_6716  ,  :narrow )
    out[:lines][:SII_6731    ] = StdEmLine(:SII_6731  ,  :narrow )

    return out
end


function QSFit.add_qso_continuum!(::Type{T}, job::JobState) where T <: WP9Type1AGN
    λ = coords(domain(job.model))

    comp = QSFit.powerlaw(median(λ))
    comp.alpha.val = -1.8

    job.model[:qso_cont] = comp
    push!(job.model[:Continuum].list, :qso_cont)
    GModelFit.update!(job.model)
end
