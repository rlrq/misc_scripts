import os
import re

def get_effects(fdir,fout):
    fnames = filter(lambda x:x[-4:]==".bed",os.listdir(fdir))
    effects = []
    for fname in fnames:
        effects.extend(re.findall("([0-9a-zA-Z_]+)[&.]",fname))
    effects = sorted(list(set(effects)))
    f = open(fout,"w+")
    f.write('\n'.join(effects))
    f.close()


def match(s,c_type,phrases):
    if c_type == "intersect" or c_type == "i":
        for phrase in phrases:
            if not phrase in s:
                return False
        return True
    elif c_type == "union" or c_type == "u":
        for phrase in phrases:
            if phrase in s:
                return True
        return False
    elif c_type == "complement" or c_type == "c":
        for phrase in phrases:
            if phrase in s:
                return False
        return True
    return False

def combine_files(fdir,fout,combine_type="intersect",*patterns):
    fnames = os.listdir(fdir)
    sel_f = [fdir+'/'+x for x in fnames if match(x,combine_type,patterns)]
    os.system("cat '{}' > {}".format("' '".join(sel_f),fout))
