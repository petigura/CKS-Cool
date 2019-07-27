import ckscool.io


def tab_star():
    df = ckscool.io.load_table('planets-cuts1+iso',cache=2)
    df['cks_sprov'] = df.cks_sprov.str.replace('smemp','emp').\
                      str.replace('smsyn','syn')

    df['rm_sb2'] = df.rm_sb2.fillna(-1)

    df = df.groupby('id_koi',as_index=False).nth(0)
    df = df.sort_values(by='id_koi')
    lines = []
    for i, row in df.iterrows():
        s = r""
        s+="{id_koi:0.0f} & "
        s+="{m17_kmag:0.1f} & "
        s+="{gaia2_sparallax:0.2f} & "
        s+="{cks_steff:0.0f} & "
        s+="{cks_smet:0.2f} & "
        s+="{cks_svsini:0.1f} & "
        s+="{gdir_srad:0.2f} & "
        s+="{giso_smass:0.2f} & "
        s+="{giso_srad:0.2f} & "
        s+="{giso_srho:0.2f} & "
        # s+="{giso_slogage:0.2f} & "

        s+="{giso2_sparallax:0.2f} & "
        s+="{cks_sprov:s} & "
        s+="{rm_sb2:0.0f} "
        # s+=r"{gaia2_gflux_ratio:0.2f} & " 
        # s+=r"{fur17_rcorr_avg:.3f} \\"
        s+=r" \\"
        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)

    return lines

def tab_planet():
    df = ckscool.io.load_table('planets-cuts1+iso',cache=2)
    df = df.sort_values(by='id_koicand')
    lines = []
    for i, row in df.iterrows():

        s = r""
        s+=r"{id_koicand:s} & "
        s+=r"{koi_period:0.1f} & "
        s+=r"{dr25_ror:.3f}  & "
        s+=r"{gdir_prad:.2f} & "  
        s+=r"{giso_sma:.3f} & "
        s+=r"{giso_sinc:.1f} \\  "

        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)

    return lines


def tab_star_csv():
    df = ckscool.io.load_table('planets-cuts1+iso',cache=2)
    df['cks_sprov'] = df.cks_sprov.str.replace('smemp','emp').\
                      str.replace('smsyn','syn')

    df['rm_sb2'] = df.rm_sb2.fillna(-1)
    df = df.groupby('id_koi',as_index=False).nth(0)
    df = df.sort_values(by='id_koi')

    s = r""
    s+="{id_koi:0.0f} "
    s+="{m17_kmag:0.1f} {m17_kmag_err:0.1f} "
    s+="{gaia2_sparallax:0.2f} {gaia2_sparallax_err:0.2f} "
    s+="{cks_steff:0.0f} {cks_steff_err:0.0f} "
    s+="{cks_smet:0.2f} {cks_smet_err:0.2f} "
    s+="{cks_svsini:0.1f} {cks_svsini_err:0.1f} "
    s+="{gdir_srad:0.2f} {gdir_srad_err1:0.2f} {gdir_srad_err2:0.2f} "
    s+="{giso_smass:0.2f} {giso_smass_err1:0.2f} {giso_smass_err2:0.2f} "
    s+="{giso_srad:0.2f} {giso_srad_err1:0.2f} {giso_srad_err2:0.2f} "
    s+="{giso_srho:0.2f} {giso_srho_err1:0.2f} {giso_srho_err2:0.2f} "
    s+="{giso2_sparallax:0.2f} {giso2_sparallax_err1:0.2f} {giso2_sparallax_err2:0.2f} "
    s+="{cks_sprov:s} "
    s+="{rm_sb2:0.0f}" # last line no space

    cols = s.split(' ')
    header = [c.split(':')[0][1:] for c in cols]
    header = ",".join(header)
    lines = []
    lines.append(header)
    for i, row in df.iterrows():
        fmt = s.replace(' ',',')
        l = fmt.format(**row)
        lines.append(l)

    return lines

def tab_planet_csv():
    df = ckscool.io.load_table('planets-cuts2+iso',cache=2)
    df = df.sort_values(by='id_koicand')

    s = r""
    s+="{id_koicand:s} "
    s+="{koi_period:0.6f} {koi_period_err1:0.6f} {koi_period_err2:0.6f} "
    s+="{dr25_ror:0.6f} {dr25_ror_err1:0.6f} {dr25_ror_err2:0.6f} "
    s+="{gdir_prad:0.3f} {gdir_prad_err1:0.3f} {gdir_prad_err2:0.3f} "
    s+="{giso_sma:0.3f} {giso_sma_err1:0.3f} {giso_sma_err2:0.3f} "
    s+="{giso_sinc:0.3f} {giso_sinc_err1:0.3f} {giso_sinc_err2:0.3f}"

    cols = s.split(' ')
    header = [c.split(':')[0][1:] for c in cols]
    header = ",".join(header)
    lines = []
    lines.append(header)
    for i, row in df.iterrows():
        fmt = s.replace(' ',',')
        l = fmt.format(**row)
        lines.append(l)

    return lines

def tab_planet_full_csv():
    df = ckscool.io.load_table('planets-cuts2+iso',cache=2)
    df = df.sort_values(by='id_koicand')
    return df.to_csv(index=False).split('\n')

def tab_field_full_csv():
    df = ckscool.io.load_table('field-cuts')
    return df.to_csv(index=False).split('\n')
