with open('./scf_log') as infile, open('mo_energies.txt', 'w') as outfile:
    copy = False
    for line in infile:
        d = line.split()
        if len(d)> 2 and d[1]=="mo_num":
            outfile.write(d[2]+"\n")

        if len(d)> 3 and d[1]=="SCF" and d[2]=="energy":
            outfile.write(d[3])
        if len(d)> 1 and d[0]=="========" and d[1]=="================" and copy==False:
            copy = True
            print(copy)
            continue
        elif len(d)> 1 and d[0]=="========" and d[1]=="================":
            copy = False
            print(copy)
            continue
        elif copy:
            outfile.write(line)
