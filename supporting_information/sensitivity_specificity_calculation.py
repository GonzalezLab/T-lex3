#!/usr/bin/python


# Created 31st May 2018
# Sensitivity and Specificity for results in Tlex2 and Tlex2.3

import sys
import csv


TEs = ["FBti0019372","FBti0019386","FBti0019415","FBti0019065","FBti0019627","FBti0020091","FBti0019430","FBti0019056","FBti0020046","FBti0018880","FBti0020119"]

pcr = list(csv.reader(open(sys.argv[1], 'r'), delimiter='\t')) # List of PCR results: pcr_results.txt
tlex = list(csv.reader(open(sys.argv[2], 'r'), delimiter='\t')) # List of Tlex2.3 results: Tlex3_results.txt

pcr_abs = 0
pcr_pres = 0
absent = 0
present = 0

te_number=1

for te in TEs:
    #te_number=11
    for p in pcr:
        strain = p[0]
        te_name = p[te_number]
        print strain
        print te_name
        print te
        for t in tlex:
            if strain in t[0] and te in t[1]:
                if te_name == "+" and t[4] != "no_data":
                    pcr_pres = pcr_pres + 1
                    if t[4] == "present" or t[4] == "polymorphic":
                        present = present +1
                        print t[4]
                elif te_name == "-" and t[4] != "no_data":
                    pcr_abs = pcr_abs + 1
                    if t[4] == "absent" or t[4] == "polymorphic":
                        absent = absent + 1
                        print t[4]
                elif te_name == "p":
                    pcr_pres = pcr_pres + 1
                    pcr_abs = pcr_abs + 1
                    if t[4] != "no_data":
                        present = present + 1
                        absent = absent + 1
                        print t[4]
    te_number = te_number+1


print pcr_pres
print pcr_abs
print present
print absent

sensitivity = float(present)/float(pcr_pres)
specificity = float(absent)/float(pcr_abs)

print("Sensitivity is " + str(sensitivity*100))
print("Specificity is " + str(specificity*100))