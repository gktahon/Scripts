#Lijst met stopcodons (capital letters)
lijststopcodons = ["TAG", "TAA", "TGA"]

#Uitvoer
codonfoundfile = open("glnA_withstopcodons.fasta", "w")
codonNOTfoundfile = open("glnA_nostopcodons.fasta", "w")

#Invoer
#bestand met DNA-sequenties overlopen
invoerdna = open("P_gregormendelii_glnA_inFrame1.fasta", "r")
codonteller = 0
for lijn in invoerdna:
    c = False
    #Start scanning
    for x in range(0,len(lijn)-2,3):
        codon = lijn[x:x+3]
        if codon.upper() in lijststopcodons:
            c = True
            #write to codonfoundfile
            codonfoundfile.write(lijn)
            print("found: ", lijn)
            codonteller += 1
            break #stop scanning
    if c == False:
        #write to codonNOTfoundfile
        codonNOTfoundfile.write(lijn)
        print("not found: ", lijn)

print("{} sequences with at least one stop codon".format(codonteller))        
invoerdna.close()
codonfoundfile.close()
codonNOTfoundfile.close()
