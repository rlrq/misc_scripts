def btop_to_cigar_generic(btop):
    output = []
    i = 0
    btop = str(btop)
    while i < len(btop):
        c = btop[i]
        if c.isdigit():
            count = ''
            while c.isdigit() and i < len(btop):
                count += c
                i += 1
                if i >= len(btop):
                    break
                c = btop[i]
            output.append((count, '='))
        else:
            count = 0
            c1, c2 = btop[i:i+2]
            c1_ref, c2_ref = c1, c2
            while i < len(btop) and not (c1.isdigit() or c2.isdigit()) and \
                  (c1_ref == c1 == '-' or c2_ref == c2 == '-' or \
                   (c1_ref != '-' != c2_ref and c1 != '-' != c2)):
                    count += 1
                    i += 2
                    if i + 2 > len(btop):
                        break
                    c1, c2 = btop[i:i+2]
            output.append((count, ('D' if c1_ref == '-' else 'I' if c2_ref == '-' else 'X')))
    return output

def btop_to_extended_cigar(btop):
    return ''.join([''.join(x) for x in btop_to_cigar_generic(btop)])

# def btop_to_cigar(btop):
#     output = ''
#     total_count = 0
#     for count, char in btop_to_cigar_generic(btop):
        
    

# def btop_to_extended_cigar_defunct(btop):
#     output = ''
#     i = 0
#     btop = str(btop)
#     while i < len(btop):
#         c = btop[i]
#         if c.isdigit():
#             while c.isdigit() and i < len(btop):
#                 output += c
#                 i += 1
#                 if i >= len(btop):
#                     break
#                 c = btop[i]
#             output += '='
#         else:
#             count = 0
#             c1, c2 = btop[i:i+2]
#             c1_ref, c2_ref = c1, c2
#             while i < len(btop) and not (c1.isdigit() or c2.isdigit()) and \
#                   (c1_ref == c1 == '-' or c2_ref == c2 == '-' or \
#                    (c1_ref != '-' != c2_ref and c1 != '-' != c2)):
#                     count += 1
#                     i += 2
#                     if i + 2 > len(btop):
#                         break
#                     c1, c2 = btop[i:i+2]
#             output += str(count) + ('D' if c1_ref == '-' else 'I' if c2_ref == '-' else 'X')
#     return output
