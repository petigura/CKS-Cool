# Code from EAP's K2-24 paper, use as template
'''
def tab_rv():
    df = ktwo24.io.load_table('rv')
    lines = []
    for i, row in df.iterrows():
        line = r""
        line+=r"{time:.6f} & {mnvel:.2f} & {errvel:.2f} \\"
        line = line.format(**row)
        lines.append(line)
    return lines

def tab_transit_times_predict():
    df = ktwo24.io.load_table('times-predict',cache=1)
    df = df.reset_index(drop=True)
    df['s_planet'] = df.i_planet.astype(str).str.replace('1','b').\
                     str.replace('2','c')
    df['date'] = pd.Series(Time(df.tc+bjd0,format='jd').iso).str.slice(stop=10)
    df['tc_err'] = 0.5 * (df.tc_err1 - df.tc_err2)
    df = df[~df.date.str.contains('2024')]
    lines = []
    for i, row in df.iterrows():
        line = r""
#        line+=r"{s_planet:s} & {i_epoch:.0f} & {date:s} & ${{{tc:.4f}}}^{{+{tc_err1:.4f}}}_{{{tc_err2:.4f}}}$ \\"
        line+=r"{s_planet:s} & {i_epoch:.0f} & {date:s} & {tc:.4f} & {tc_err:.4f} \\"
        line = line.format(**row)
        lines.append(line)
    return lines

'''

