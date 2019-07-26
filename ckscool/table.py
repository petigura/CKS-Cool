import ckscool.io

def tab_planet():
    df = ckscool.io.load_table('ckscool-planets',cache=2)
    df = df.sort_values(by='id_koicand')
    lines = []
    for i, row in df.iterrows():

        s = r""
        s+=r"{id_koicand:s} & "
        s+=r"{koi_period:0.1f} & "
        s+=r"{koi_ror:.3f}  & "
        s+=r"{gdir_prad:.2f} & "  
        s+=r"{giso_sma:.3f} & "
        s+=r"{giso_sinc:.0f} \\  "

        s = s.format(**row)
        s = s.replace('nan','\\nodata')
        lines.append(s)

    return lines

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

def tab_star_machine():
    df = ckscool.io.load_table('cksgaia-planets')
    df = df.groupby('id_starname', as_index=False).nth(0)
    df = df.sort_values(by='id_starname')

    cols = ['id_starname',
            'cks_steff', 'cks_steff_err1', 'cks_steff_err2',
            'cks_smet', 'cks_smet_err1', 'cks_smet_err2',
            'm17_kmag', 'm17_kmag_err',
            'gaia2_sparallax', 'gaia2_sparallax_err',
            'gdir_srad', 'gdir_srad_err1', 'gdir_srad_err2',
            'giso_smass', 'giso_smass_err1', 'giso_smass_err2',
            'giso_slogage', 'giso_slogage_err1', 'giso_slogage_err2',
            'gaia2_gflux_ratio', 'fur17_rcorr_avg']

    df = df[cols]

    lines = []
    head = ",".join(cols)
    lines.append(head)
    for i, row in df.iterrows():
        l = "{id_starname:s},\
{cks_steff:.0f}, {cks_steff_err1:.0f}, {cks_steff_err2:.0f},\
{cks_smet:.2f},{cks_smet_err1:.2f},{cks_smet_err2:.2f},\
{m17_kmag:.3f},{m17_kmag_err:.3f},\
{gaia2_sparallax:.3f},{gaia2_sparallax_err:.3f},\
{gdir_srad:.3f},{gdir_srad_err1:.3f},{gdir_srad_err2:.3f},\
{giso_smass:.3f},{giso_smass_err1:.3f},{giso_smass_err2:.3f},\
{giso_slogage:.2f},{giso_slogage_err1:.2f},{giso_slogage_err2:.2f},\
{gaia2_gflux_ratio:.3f},{fur17_rcorr_avg:.4f}".format(**row)

        lines.append(l)

    return lines

def tab_planet_machine():
    df = ckscool.io.load_table('cksgaia-planets')
    df = df.sort_values(by='id_koicand')

    cols = ['id_koicand',
            'koi_period', 'koi_period_err1', 'koi_period_err2',
            'koi_ror', 'koi_ror_err1', 'koi_ror_err2',
            'gdir_prad', 'gdir_prad_err1', 'gdir_prad_err2',
            'giso_sma', 'giso_sma_err1', 'giso_sma_err2',
            'giso_insol', 'giso_insol_err1', 'giso_insol_err2']

    lines = []
    head = ",".join(cols)
    lines.append(head)
    for i, row in df.iterrows():
        l = "{id_koicand:s},\
{koi_period:.9f}, {koi_period_err1:.9f}, {koi_period_err2:.9f},\
{koi_ror:.6f}, {koi_ror_err1:.6f}, {koi_ror_err2:.6f},\
{gdir_prad:.3f},{gdir_prad_err1:.3f},{gdir_prad_err2:.3f},\
{giso_sma:.5f},{giso_sma_err1:.5f},{giso_sma_err2:.5f},\
{giso_insol:.1f},{giso_insol_err1:.1f},{giso_insol_err2:.1f}".format(**row)

        lines.append(l)

    return lines
