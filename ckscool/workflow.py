import os
from collections import OrderedDict
from matplotlib import pylab as plt

class Workflow(object):
    def __init__(self):
        self.plot = OrderedDict()
        self.table = OrderedDict()
        self.val = OrderedDict()
        self.csv = OrderedDict()
        self.build_dir = "./build"
        self.paper_dir = "./paper"

    def get_all_dict(self):
        d = OrderedDict()
        d['plot'] = self.plot
        d['table'] = self.table
        d['val'] = self.val
        d['csv'] = self.csv
        return d
        
    def key2fn(self, key, kind):
        if kind=='plot':
            fn = 'fig_'+key+'.pdf'
        elif kind=='table':
            fn = 'tab_'+key+'.tex'
        elif kind=='csv':
            fn = 'tab_'+key+'.csv'
        elif kind=='val':
            fn = 'val_'+key+'.tex'
            
        fn = os.path.join(self.build_dir,fn)
        return fn

    def create_file(self, kind, name):
        i = 0
        all_dict = self.get_all_dict()
        os.system('mkdir -p {}'.format(self.build_dir))
        for key, func in all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif key.count(name)==1:
                    func()
                else:
                    continue
                    
                fn = self.key2fn(key, 'plot')
                plt.gcf().savefig(fn,dpi=300)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif key.count(name)==1:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'table')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='csv':
                if name=='all':
                    lines = func()
                elif key.count(name)==1:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'csv')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{",
                ]

                lines2 = [
                    "}[XX]",
                    "}"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("%\n".join(lines))

            i+=1

        if i==0:
            assert False, name + " not a valid key"

    def update_paper(self):
        all_dict = self.get_all_dict()
        os.system('mkdir -p {}'.format(self.paper_dir))
        for kind, d in all_dict.iteritems():
            for key, val in d.iteritems():
                src = self.key2fn(key, kind)
                cmd = 'cp {} {}'.format(src, self.paper_dir)
                print cmd
                os.system(cmd)
