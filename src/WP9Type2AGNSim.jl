abstract type WP9Type2AGNSim <: WP9Type2AGN end

function QSFit.Options(::Type{T}) where T <: WP9Type2AGNSim
    out = QSFit.Options(supertype(T))
    out[:n_unk] = 2
    out[:line_profiles] = :gauss

    out[:lines] = OrderedDict{Symbol, QSFit.EmLineDescription}()
    out[:lines][:CIV_1549  ] = StdEmLine(:CIV_1549  , :narrow)
    out[:lines][:CIII_1909 ] = StdEmLine(:CIII_1909 , :narrow)
    out[:lines][:MgII_2798 ] = StdEmLine(:MgII_2798 , :narrow)
    out[:lines][:NeV_3426  ] = StdEmLine(:NeV_3426  , :narrow)
    out[:lines][:OII_3727  ] = StdEmLine(:OII_3727  , :narrow)
    out[:lines][:NeIII_3869] = StdEmLine(:NeIII_3869, :narrow)
    out[:lines][:Hg        ] = StdEmLine(:Hg        , :narrow)
    out[:lines][:Hb        ] = StdEmLine(:Hb        , :narrow)
    out[:lines][:OIII_4959 ] = StdEmLine(:OIII_4959 , :narrow)
    out[:lines][:OIII_5007 ] = StdEmLine(:OIII_5007 , :narrow)
    out[:lines][:OIII_5007_bw]=StdEmLine(:OIII_5007 , :narrow)
    out[:lines][:OI_6300   ] = StdEmLine(:OI_6300   , :narrow)
    out[:lines][:OI_6364   ] = StdEmLine(:OI_6364   , :narrow)
    out[:lines][:Ha        ] = StdEmLine(:Ha        , :broad )
    out[:lines][:SII_6716  ] = StdEmLine(:SII_6716  , :narrow)
    out[:lines][:SII_6731  ] = StdEmLine(:SII_6731  , :narrow)
    return out
end

function QSFit.EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Val{:broad}) where T <: WP9Type2AGNSim
    lc = QSFit.EmLineComponent(supertype(T), job, λ, Val(:broad)) # invoke parent recipe
    lc.comp.fwhm.val = 2000
    lc.comp.fwhm.high = 3500
    return lc
end

