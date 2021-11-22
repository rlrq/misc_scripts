def gc_content(seq):
    seq = str(seq).upper().replace('-', '')
    gc_count = seq.count('G') + seq.count('C')
    return gc_count/len(seq)
    
