abstract type WP9Type1AGNSim <: WP9Type1AGN end

function QSFit.Options(::Type{T}) where T <: WP9Type1AGNSim
    out = QSFit.Options(supertype(T))
    out[:n_unk] = 2

    for (name, line) in out[:lines]
        if length(line.types) > 1
            out[:lines][name] = StdEmLine(line.tid, :broad)
        end
    end

    lines = OrderedDict{Symbol, QSFit.EmLineDescription}()
    lines[:OVI_1034    ] = CustomEmLine(1033.82 ,  :broad  )
    lines[:Lya         ] = StdEmLine(:Lya       ,  :broad  )
    lines[:HeII_1215   ] = StdEmLine(:HeII_1215 ,  :broad  )
    lines[:NV_1241     ] = StdEmLine(:NV_1241   ,  :broad  )
    # lines[:OI_1306     ] = StdEmLine(:OI_1306   ,  :broad  )
    # lines[:CII_1335    ] = StdEmLine(:CII_1335  ,  :broad  )
    # lines[:SiIV_1400   ] = StdEmLine(:SiIV_1400 ,  :broad  )
    lines[:NIV_1483    ] = CustomEmLine(1483.32 ,  :narrow )
    lines[:NIV_1487    ] = CustomEmLine(1486.5  ,  :narrow )
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  ,  :broad  )
    lines[:HeII_1640   ] = StdEmLine(:HeII_1640 ,  :broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  :broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  :broad  )
    # lines[:AlIII_1858  ] = StdEmLine(:AlIII_1858,  :broad  )
    lines[:NIII_1750   ] = CustomEmLine(1750    ,  :broad  )
    lines[:SiIII_1883  ] = CustomEmLine(1882.71 ,  :broad  )
    lines[:SiIII_1892  ] = CustomEmLine(1892.03 ,  :broad  )
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :broad  )
    lines[:CII_2326    ] = StdEmLine(:CII_2326  ,  :broad  )
    lines[:NeIV_2424   ] = CustomEmLine(2424.0  ,  :narrow )
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :broad  )
    lines[:NeIII_3342  ] = CustomEmLine(3342.18 ,  :narrow )
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow )
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow )
    lines[:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow )
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow )
    lines[:Hd          ] = StdEmLine(:Hd        ,  :broad  )
    lines[:Hg          ] = StdEmLine(:Hg        ,  :broad  )
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow )
    lines[:HeII_4686   ] = StdEmLine(:HeII_4686 ,  :broad  )
    lines[:Hb          ] = StdEmLine(:Hb        ,  :broad  )
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow )
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:HeI_5876    ] = StdEmLine(:HeI_5876  ,  :broad  )
    lines[:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow )
    lines[:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow )
    lines[:Ha          ] = StdEmLine(:Ha        ,  :broad  )
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
    lines[:Pab         ] = StdEmLine(:Pab       ,  :broad  )  # Paschen β

    out[:lines] = lines
    return out
end
