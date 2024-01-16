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

    lines = OrderedDict{Symbol, QSFit.EmLineDescription}()
    lines[:OVI_1034    ] = CustomEmLine(1033.82 ,  :narrow, :broad)
    lines[:Lya         ] = StdEmLine(:Lya       ,  :narrow, :broad)
    lines[:HeII_1215   ] = StdEmLine(:HeII_1215 ,  :broad )
    lines[:NV_1241     ] = StdEmLine(:NV_1241   ,  :broad )
    # lines[:OI_1306     ] = StdEmLine(:OI_1306   ,  :broad  )
    # lines[:CII_1335    ] = StdEmLine(:CII_1335  ,  :broad  )
    # lines[:SiIV_1400   ] = StdEmLine(:SiIV_1400 ,  :broad  )
    lines[:NIV_1483    ] = CustomEmLine(1483.32 ,  :narrow )
    lines[:NIV_1487    ] = CustomEmLine(1486.5  ,  :narrow )
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  ,  :narrow, :broad)
    lines[:HeII_1640   ] = StdEmLine(:HeII_1640 ,  :broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  :broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  :broad  )
    # lines[:AlIII_1858  ] = StdEmLine(:AlIII_1858,  :broad  )
    lines[:NIII_1750   ] = CustomEmLine(1750    ,  :broad  )
    lines[:SiIII_1883  ] = CustomEmLine(1882.71 ,  :broad  )
    lines[:SiIII_1892  ] = CustomEmLine(1892.03 ,  :broad  )
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :broad  )
    lines[:CII_2326    ] = StdEmLine(:CII_2326  ,  :broad  )
    lines[:NeIV_2424   ] = CustomEmLine(2424.0  ,  :narrow  )
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :narrow, :broad)
    lines[:NeIII_3342  ] = CustomEmLine(3342.18 ,  :narrow )
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow )
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow )
    lines[:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow )
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow )
    lines[:Hd          ] = StdEmLine(:Hd        ,  :broad  )
    lines[:Hg          ] = StdEmLine(:Hg        ,  :broad  )
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow )
    lines[:HeII_4686   ] = StdEmLine(:HeII_4686 ,  :broad  )
    lines[:Hb          ] = StdEmLine(:Hb        ,  :narrow, :broad)
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow )
    lines[:OIII_4959_bw] = StdEmLine(:OIII_4959 ,  :narrow )
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:HeI_5876    ] = StdEmLine(:HeI_5876  ,  :broad  )
    lines[:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow )
    lines[:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow )
    lines[:NII_6549    ] = StdEmLine(:NII_6549  ,  :narrow )
    lines[:Ha          ] = StdEmLine(:Ha        ,  :narrow, :broad, :verybroad)
    lines[:NII_6583    ] = StdEmLine(:NII_6583  ,  :narrow )
    lines[:SII_6716    ] = StdEmLine(:SII_6716  ,  :narrow )
    lines[:SII_6731    ] = StdEmLine(:SII_6731  ,  :narrow )
    lines[:ArIV_7171   ] = CustomEmLine(7170.7  ,  :narrow )
    lines[:OII_7325    ] = CustomEmLine(7325    ,  :narrow )  # Blend of 4 lines 
    lines[:Pa10        ] = CustomEmLine(9016    ,  :broad  )  # Are Paschen series broad + narrow or? 
    lines[:SIII_9069   ] = CustomEmLine(9068.62 ,  :narrow ) 
    lines[:Pa9         ] = CustomEmLine(9228.97 ,  :broad  )  
    lines[:SIII_9531   ] = CustomEmLine(9530.62 ,  :narrow )  
    lines[:Pa8         ] = CustomEmLine(9545.93 ,  :broad  )
    lines[:CI_9850     ] = CustomEmLine(9850    ,  :narrow )  # Blend of 3 lines 
    lines[:Pad         ] = CustomEmLine(10049   ,  :broad  )  # Paschen δ
    lines[:HeII_10126  ] = CustomEmLine(10126   ,  :broad  ) 
    lines[:HeI_10830   ] = CustomEmLine(10830.3 ,  :broad  )
    lines[:Pag         ] = CustomEmLine(10938   ,  :broad  )  # Paschen γ
    lines[:OI_11287    ] = CustomEmLine(11287   ,  :narrow )
    lines[:PII_11886   ] = CustomEmLine(11886   ,  :narrow )
    lines[:Pab         ] = StdEmLine(:Pab       ,  :narrow, :broad)  # Paschen β

    out[:lines] = lines
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
