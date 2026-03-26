import tskit
import argparse

def ts_analyse_file(filename, n=25, binwidth=5):
    ts = tskit.load(filename)
    ts_analyse(ts, filename, n, binwidth)

def ts_analyse(ts, name="name", n=24, binwidth = 5):
    treelengths = []
    for tree in ts.trees():
        distt = int((tree.interval[1] - tree.interval[0]) // 1) # simplify... but still float
        treelengths.append(distt)
    lset = sorted(list(set(treelengths)))
    lcounts = [(val, len([t for t in treelengths if t == val])) for val in lset]
    countprint(lcounts, name, n, binwidth)
    return(lcounts)

def countprint(lcounts, name="name", num=50, binwidth = 5):
    ncounts = list(reversed(sorted(lcounts, key=lambda x: x[1])))
    maxpair = ncounts[0]
    [print(maxpair)]
    norm = 100
    scalefactor = norm / (ncounts[0][1] * binwidth * 0.5)

    # make the major one... first 25
    ocounts = list(sorted(lcounts, key=lambda x: x[0]))
    print(f'longest 5 sequences:\t{ocounts[-5:]}')
    print(f"~~ Analysing tree sequence: {name}\n")
    print("----- ORDERED BY POSITION -----")
    bincount = 1
    low = 1
    count = 0
    scale=True
    subscale = 1
    for n, opair in enumerate(ocounts):
        if n >= num*binwidth:
            break
        lab = opair[0]
        if lab > num * binwidth:
            break
        if lab > bincount * binwidth:
            obase = f'{low}-{low + binwidth - 1}:'
            for _ in range(2 - (len(obase)) // 4):
                obase += '\t'
            ostring = f"{obase}\t{count}\t"
            if scale and count > 0:
                subscale = norm / (count * scalefactor)
                scale = False
            for n in range(int((count * scalefactor * subscale)//1)):
                ostring +="|"
            print(ostring)
            while lab > bincount * binwidth:
                low = bincount * binwidth + 1
                bincount += 1
                count = 0
                if lab > bincount * binwidth:
                    obase = f'{low}-{low + binwidth - 1}:'
                    for _ in range(2 - (len(obase)) // 4):
                        obase += '\t'
                    ostring = f"{obase}\t{count}\t"
                    for n in range(int((count * scalefactor)//1)):
                        ostring +="|"
                    print(ostring)
                    
        count += opair[1]
    # scale = True
    # print("----- MOST FREQUENT LENGTHS -----")
    # for n in range(10):
    #     cpair = ncounts[n]
    #     lab = cpair[0]
    #     count = cpair[1]
    #     ostring = f"{lab}:\t{count}\t"
    #     if scale:
    #         subscale = norm / (count * scalefactor)
    #         scale = False
    #     for n in range(int((count * scalefactor * subscale)//1)):
    #         ostring +="|"
    #     print(ostring)


def main():
    parser = argparse.ArgumentParser(description="Visualize basic attributes of .trees file")
    parser.add_argument("-tree", type=str, help=".trees file to analyse")
    parser.add_argument("-n", type=int, default=25, help="number of lines in frequency graph")
    parser.add_argument("-binwidth", type=int, default=5, help="binwidth of frequency graph")
    args = parser.parse_args()
    ts_analyse_file(args.tree, args.n, args.binwidth)



if __name__ == "__main__":
    main()