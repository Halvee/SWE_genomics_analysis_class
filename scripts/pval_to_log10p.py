

import sys
import math

## PARAM
P_0_SUBSTITUTE=4e-324

def main():

    try:
        ARGS = sys.argv[1:]
        infile = ARGS[0]
        pval_column = int(ARGS[1])
    except:
        print("pval_to_log10p.py <infile|stdin> <pval_column>")
        sys.exit(1)

    # open filehandle
    if infile == "stdin":
        in_fh = sys.stdin
    else:
        in_fh = open(infile, "r")

    # for each line in input file ..
    i = 0
    for line in in_fh:
        
        # increment i
        i += 1

        # split on whitespace 
        data = line.rstrip().split()

        # if this is row 1, then just print and go to next line
        if i == 1:
            print("\t".join(data))
            continue

        # get column with p-value in it
        pval = float(data[pval_column-1])

        # if pval too low (ie. 0) then rewrite to value that won't cause issues
        if pval < P_0_SUBSTITUTE:
            pval = P_0_SUBSTITUTE

        # transform p to -log10(p) and insert into row
        neg_log10_pval = -(math.log10(pval))
        data[pval_column-1] = neg_log10_pval

        # print row to stdout
        out = "\t".join([str(x) for x in data])
        print(out)

    return

if __name__ == "__main__":
    main()
