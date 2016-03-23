#!/usr/bin/env python

import sys
line = sys.stdin.readline()
ba = line.split()
with open(sys.argv[1]) as x:
    x_line = x.readline()
    xa = x_line.split()

    while line and x_line:

        if (int(ba[1]) >= int(ba[2])):
            line = sys.stdin.readline()
            ba = line.split()
            continue

        # Check for different chromosomes
        if (xa[0] != ba[0]):
            if (xa[0] < ba[0]):
                x_line = x.readline()
                xa = x_line.split()
            else:
                print("\t".join(ba))
                line = sys.stdin.readline()
                ba = line.split()
            continue

        # print and fast-forward bed
        if (int(ba[2]) <= int(xa[1])):
            print("\t".join(ba))
            line = sys.stdin.readline()
            ba = line.split()
            continue

        # fast-forward excludes
        if (int(xa[2]) <= int(ba[1])):
            x_line = x.readline()
            xa = x_line.split()
            continue

        # If exclude enitrely covers bed, skip to next bed
        if int(xa[2]) >= int(ba[2]) and int(xa[1]) <= int(ba[1]):
            line = sys.stdin.readline()
            ba = line.split()
            continue

        # If bed entirely covers exclude, print out first part and modify
        if int(ba[2]) >= int(xa[1]) and int(ba[1]) <= int(xa[1]):
            print("\t".join([ba[0], ba[1], xa[1], ba[3], ba[4]]))
            ba[1] = xa[2]
            continue

        # If exclude covers first part of bed, mask it and go to next exclude
        if (int(xa[2]) > int(ba[1])):
            ba[1] = xa[2]
            continue

        # If exclude covers last part of bed, print modified and read next
        if (int(xa[1]) <= int(ba[2])):
            ba[2] = xa[1]
            print("\t".join(ba))
            continue

        raise Exception("Coding error: Should never get here")

# Print out remaining lines
while line:
    sys.stdout.write(line)
    line = sys.stdin.readline()
