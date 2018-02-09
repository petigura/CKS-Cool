# Code from EAP's K2-24 paper use as a template


def val_stat(return_dict=False):
    """
    Statistics of sample
    """
    d = OrderedDict()
    # load up data from p16
    

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d

    return lines

