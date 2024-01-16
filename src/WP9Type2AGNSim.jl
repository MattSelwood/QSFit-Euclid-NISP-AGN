abstract type WP9Type2AGNSim <: WP9Type2AGN end

function QSFit.Options(::Type{T}) where T <: WP9Type2AGNSim
    out = QSFit.Options(supertype(T))
    out[:n_unk] = 2
    out[:line_profiles] = :gauss

    lines = OrderedDict{Symbol, QSFit.EmLineDescription}()
    lines[:Lya         ] = StdEmLine(:Lya       ,  :narrow)
    lines[:NIV_1483    ] = CustomEmLine(1483.32 ,  :narrow)
    lines[:NIV_1487    ] = CustomEmLine(1486.5  ,  :narrow)
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  , :narrow)
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 , :narrow)
    lines[:NeIV_2424   ] = CustomEmLine(2424.0  ,  :narrow) 
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 , :narrow)
    lines[:NeIII_3342  ] = CustomEmLine(3342.18 ,  :narrow)
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow)
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow)
    lines[:OII_3727    ] = StdEmLine(:OII_3727  , :narrow)
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869, :narrow)
    lines[:Hd          ] = StdEmLine(:Hd        ,  :narrow)
    lines[:Hg          ] = StdEmLine(:Hg        ,  :narrow)
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow)
    lines[:Hb          ] = StdEmLine(:Hb        ,  :narrow)
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 , :narrow)
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 , :narrow)
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 , :narrow)
    lines[:OI_6300     ] = StdEmLine(:OI_6300   , :narrow)
    lines[:OI_6364     ] = StdEmLine(:OI_6364   , :narrow)
    lines[:Ha          ] = StdEmLine(:Ha        , :broad )
    lines[:SII_6716    ] = StdEmLine(:SII_6716  , :narrow)
    lines[:SII_6731    ] = StdEmLine(:SII_6731  , :narrow)
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

function QSFit.EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Val{:broad}) where T <: WP9Type2AGNSim
    lc = QSFit.EmLineComponent(supertype(T), job, λ, Val(:broad)) # invoke parent recipe
    lc.comp.fwhm.val = 2000
    lc.comp.fwhm.high = 3500
    return lc
end

