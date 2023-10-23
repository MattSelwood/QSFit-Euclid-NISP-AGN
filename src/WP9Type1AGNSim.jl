abstract type WP9Type1AGNSim <: WP9Type1AGN end

function QSFit.Options(::Type{T}) where T <: WP9Type1AGNSim
    out = QSFit.Options(supertype(T))
    out[:n_unk] = 2

    for (name, line) in out[:lines]
        if length(line.types) > 1
            out[:lines][name] = StdEmLine(line.tid, :broad)
        end
    end

    out[:lines] = OrderedDict{Symbol, QSFit.EmLineDescription}()
    out[:lines][:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :broad  )
    out[:lines][:CII_2326    ] = StdEmLine(:CII_2326  ,  :broad  )
    out[:lines][:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :broad  )
    out[:lines][:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow )
    out[:lines][:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow )
    out[:lines][:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow )
    out[:lines][:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow )
    out[:lines][:Hd          ] = StdEmLine(:Hd        ,  :broad  )
    out[:lines][:Hg          ] = StdEmLine(:Hg        ,  :broad  )
    out[:lines][:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow )
    out[:lines][:HeII_4686   ] = StdEmLine(:HeII_4686 ,  :broad  )
    out[:lines][:Hb          ] = StdEmLine(:Hb        ,  :broad  )
    out[:lines][:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow )
    out[:lines][:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow )
    out[:lines][:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow )
    out[:lines][:HeI_5876    ] = StdEmLine(:HeI_5876  ,  :broad  )
    out[:lines][:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow )
    out[:lines][:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow )
    out[:lines][:Ha          ] = StdEmLine(:Ha        ,  :broad  )
    out[:lines][:SII_6716    ] = StdEmLine(:SII_6716  ,  :narrow )
    out[:lines][:SII_6731    ] = StdEmLine(:SII_6731  ,  :narrow )

    return out
end
