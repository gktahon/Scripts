'''
Created on 24-jan.-2013
Plaatst elke sequentie in reading frame 1 op basis van de correcte aminozuursequentie.
@author: pstragie
Version: 0.3

#!/usr/bin/env python
#
# FastaFrameSetter.py
# written by
###################
# Pieter Stragier #
###################
#
# DEPENDENCIES
# ============
# None
#
#
# USAGE
# =====
#
# FastaFrameSetter.py [options]
#
# Options:
# -h, --help show this help message and exit
# -o OUTDIRNAME, --outdir=OUTDIRNAME
# Output directory
# -i INDIRNAME, --indir=INDIRNAME
# Input directory name
# -v, --verbose Give verbose output
# -f, --force Force file overwriting
# --noclobber Don't nuke existing files
#=============
'''

#### Veranderen indien van toepassing #####
#Dit bestand bevat de codon omzettingstabel
codoncode = 'standaard_code.txt'



######### Script (afblijven) ############


import os, shutil, sys, logging, logging.handlers, time
from optparse import OptionParser

def parse_cmdline(args):
    """ Parse command-line arguments for the script
    """
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", dest="outdirname",
                        action="store", default=None,
                        help="Output directory")
    parser.add_option("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Location of input file")
    
    parser.add_option("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_option("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Force file overwriting")
    
    parser.add_option("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    
    parser.add_option("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't nuke existing files")
    
    return parser.parse_args()


def make_outdir():
    """ Make the output directory, if required.
        This is a  little involved. If the output directory already exists,
        we take the safe option by default, and stop with an error. We can,
        however, choose to force  the program to go on, in which  case we can
        either clobber or not the existing directory. The options turn out
        as the following, if the directory exists:

        DEFAULT:    stop
        FORCE: continue,  and remove the existing    output directory
        NOCLOBBER+FORCE:    continue,  but do not remove the existing    output
    """
    if os.path.exists(options.outdirname):
        if not options.force:
            logger.error("Output directory {} would ".format(options.outdirname + "overwrite existing files (exiting)"))
            sys.exit(1)
        else:
            logger.info("Removing directory {} and everything below it".format(options.outdirname))
            if options.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(options.outdirname)
    logger.info("Creating directory {}".format(options.outdirname))
    try:
        os.makedirs(options.outdirname) #  We make    the directory  recursively
    except OSError:
        #  This    gets    thrown if the directory  exists. If we've  forced overwrite/
        #  delete and we're  not clobbering, we let things slide
        if options.noclobber  and options.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)
            

def ambigue(codon):
    """Ambigue codons omzetten indien mogelijk"""
    amb_lijst = ['W', 'S', 'R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    AA1 = ""
    AA2 = ""
    AA3 = ""
    AA4 = ""
    for i in codon:
        if i in ["A", "G", "C", "T"]:
            AA1 += i
            AA2 += i
            AA3 += i
            AA4 += i
        elif i == "W":
            AA1 += "A"
            AA2 += "T"
            AA3 += "A"
            AA4 += "T"
        elif i == "S":
            AA1 += "G"
            AA2 += "C"
            AA3 += "G"
            AA4 += "C"
        elif i == "R":
            AA1 += "A"
            AA2 += "G"
            AA3 += "A"
            AA4 += "G"
        elif i == "Y":
            AA1 += "C"
            AA2 += "T"
            AA3 += "C"
            AA4 += "T"
        elif i == "K":
            AA1 += "G"
            AA2 += "T"
            AA3 += "G"
            AA4 += "T"
        elif i == "M":
            AA1 += "A"
            AA2 += "C"
            AA3 += "A"
            AA4 += "C"
        elif i == "B":
            AA1 += "C"
            AA2 += "G"
            AA3 += "T"
            AA4 += "C"
        elif i == "D":
            AA1 += "A"
            AA2 += "G"
            AA3 += "T"
            AA4 += "A"
        elif i == "H":
            AA1 += "A"
            AA2 += "C"
            AA3 += "T"
            AA4 += "A"
        elif i == "V":
            AA1 += "A"
            AA2 += "C"
            AA3 += "G"
            AA4 += "A"
        elif i == "N":
            AA1 += "A"
            AA2 += "C"
            AA3 += "G"
            AA4 += "T"
    return [AA1, AA2, AA3, AA4]


class GenetischeCode:
    '''Genetische code'''
    
    def __init__(self, bestand):
        self.codon_file = bestand
        codons = open(self.codon_file, 'r')
        codon_dict = {}
        for regel in codons:
            code = regel.split()
            codon_dict[code[0].upper()] = code[1]
            RNA_streng = ''
            for i in code[0]:
                if i == "T" or i == "t":
                    RNA_streng += 'U'
                elif i == 'U' or i == 'u':
                    RNA_streng += 'T'
                else:
                    RNA_streng += i.upper()
            codon_dict[RNA_streng.upper()] = code[1]
                    
        self.code = codon_dict    


    def aminozuur(self, codon):
        assert(codon.upper() in self.code), "'{}' is geen geldig codon.".format(codon)
        return self.code[codon.upper()]


    def eiwit(self, streng):

        alfabet = "ACTUGWSRYKMBDHVN\*"
        result = ""
        
        for i in streng:
            assert(i.upper() in alfabet), "ongeldige DNA- of RNA-sequentie: {}".format(i)
        
        L = len(streng)
        if L % 3 == 0:
            streng = streng.upper()
        elif L % 3 == 1:
            streng = streng[:-1].upper()
        else:
            streng = streng[:-2].upper()
        for i in range(0, len(streng), 3):
            if not streng[i:i+3].upper() in self.code.keys():
                AA = ambigue(streng[i:i+3])
                if len(AA) == 2:
                    if self.code[AA[0]] == self.code[AA[1]]:
                        eiwit = self.code[AA[0]]
                        result += eiwit
                    else:
                        result += "?"
                elif len(AA) == 3:
                    if self.code[AA[0]] == self.code[AA[1]] == self.code[AA[2]]:
                        eiwit = self.code[AA[0]]
                        result += eiwit
                    else:
                        result += "?"
                elif len(AA) == 4:
                    if self.code[AA[0]] == self.code[AA[1]] == self.code[AA[2]] == self.code[AA[3]]:
                        eiwit = self.code[AA[0]]
                        result += eiwit
                    else:
                        result += "?"
                else:
                    result += "?"
            else:
                eiwit = self.code[streng[i:i+3].upper()]
                result += eiwit
        
        return result


def complement(nucleotide):
    '''
    Geef de complementaire nucleotide terug
    '''
    if nucleotide == "G":
        return "C"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "R":
        return "Y"
    elif nucleotide == "Y":
        return "R"
    elif nucleotide == "W":
        return "W"
    elif nucleotide == "S":
        return "S"
    elif nucleotide == "K":
        return "M"
    elif nucleotide == "M":
        return "K"
    elif nucleotide == "B":
        return "V"
    elif nucleotide == "D":
        return "H"
    elif nucleotide == "H":
        return "D"
    elif nucleotide == "V":
        return "B"
    elif nucleotide == "N":
        return "N"
    else:
        return "-"


def reverseComplement(sequence):
    '''
    Geef de reverse complementaire sequentie terug
    '''
    inverscomplement = ""
    for nucleotide in sequence[::-1]:
        inverscomplement += complement(nucleotide)
    return inverscomplement



def motieven(sequentie, subsequentie, window=11, begin=1, einde=0):
    if einde == 0:
        einde = len(subsequentie)
    a = begin-1
    b = einde
    r = False

    for j in range(0, len(subsequentie)-window):
        #allow mismatches and calculate similarity
        query = subsequentie[j:j+window]
        for k in range(0, len(sequentie)-window):
            sim = 0
            subject = sequentie[k:k+window]
            for i, l in enumerate(subject):
                if l == query[i]:
                    sim += 1
            percsim = sim/window
            if percsim >= 0.80:
                r = True
                return r
    return r


def verifyTemplate(AAtemplate, DNAtemplate):
    #Ask verification of AA sequence
    correct = False
    while correct == False:
        stopcodons = len([x for x in AAtemplate if x == "*"])
        if stopcodons == 1:
            cod = "codon"
        else:
            cod = "codons"
        cor = input("{} stop {} found. Correct frame? Continue? (y/n)".format(stopcodons, cod))

        if not (cor == "y" or cor == "yes"):
            nc = '0'
            while nc != '2' and nc != '3' and nc != '1':
                nc = input("Try other frame? (2 or 3): ")
            AAtemplate = code.eiwit(template_DNA[int(nc)-1:])
            print(AAtemplate)
            stopcodons = len([x for x in AAtemplate if x == "*"])
            if stopcodons == 1:
                cod = "codon"
            else:
                cod = "codons"

        else:
            nc = 1
            if cor == 'y' or cor == 'yes':
                break
                correct = True
            else:
                correct = False
    return (AAtemplate, DNAtemplate[nc-1:])


def tempFasta(directory, fastabestand):
    invoer = open(fastabestand, "r")
    #Remove \n from between sequences
    outfile = open(os.path.join(directory, 'fastaTmp.txt'), 'w')

    totalsequences = 0
    
    
    d = {}
    sequence = ""
    for regel in invoer:
        if regel.startswith(">"):
            if len(sequence) > 0:
                if not header in d.keys():
                    d[header] = sequence
                    sequence = ""
                else:
                    seq = d[strain]
                    seq += sequence
                    d[header] = seq
                    sequence = ""
            header = regel.rstrip("\n")
        else:
            sequence += regel.rstrip("\n").upper()
    d[header] = sequence
    invoer.close()
    """
    n = 80 #tekens per lijn
    for k, v in d.items():
        outfile.write("{}\n".format(k))
        L = len(v)
        for x in range(L//n):
            lijn = v[x*n:n*(x+1)]
            outfile.write(lijn+"\n")
        #laatste lijn toevoegen
        lijn = v[(x+1)*n:]
        outfile.write(lijn+"\n")
    """
    for k, v in d.items():
        outfile.write("{}\n{}\n".format(k,v))
        
    
    outfile.close()
    print("Total number of sequences found in file: {}".format(totalsequences))
    return os.path.join(directory, 'fastaTmp.txt')


def collectSequences(bestand):
    invoer = open(bestand, 'r')
    sequence = ""
    lijst = []
    for line in invoer:
        if line.startswith(">"):
            if not sequence == "":
                lijst.append((header, sequence))
                sequence = ""
            header = line.lstrip(">").rstrip("\n")
        else:
            if header:
                sequence = line.rstrip("\n")
    lijst.append((header, sequence))
    invoer.close()
    return lijst


if __name__ == '__main__':
    # Parse command-line
    # Options are all options - no arguments
    options, args = parse_cmdline(sys.argv)
    #  We set up logging,    and modify loglevel    according  to whether we need
    #  verbosity  or not
    #  err_handler points to sys.stderr
    #  err_handler_file    points to a  logfile,    if named
    logger = logging.getLogger('FastaFrameSetter.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.logfile is not None:
        try:
            logstream = open(options.logfile,  'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                            options.logfile)
            sys.exit(1)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# FastaFrameSetter.py logfile')
    logger.info('# Run: %s' % time.asctime())

    #  Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    #  Have    we got an input  and output directory? If not,    exit.
    if options.indirname  is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    
    logger.info("Input directory: {}".format(options.indirname))
    if options.outdirname is None:
        logger.error("No output directory name (Input directory will be used)")
        options.outdirname = os.path.dirname(options.indirname)
    if options.outdirname != os.path.dirname(options.indirname):
        make_outdir()
    
    logger.info("Output directory: {}".format(options.outdirname))

    directory = os.path.dirname(options.indirname)
    base = os.path.basename(options.indirname)
    shortname, extension = base.split(".") # split inputfilename into name and extension
    outputfile = shortname + "_inFrame1.fasta" # change name
    outputfilenoframe = shortname + "_noFrame.fasta" # change name
    
    #ask for DNA template in the correct frame
    template_DNA = str(input("Geef template DNA in reading frame +1:"))

    #translate template to amino acid sequence
    code = GenetischeCode(codoncode)
    AAtemplate = code.eiwit(template_DNA)
    logger.info(AAtemplate)

    #Ask verification of AA sequence
    AAtemplate, DNAtemplate = verifyTemplate(AAtemplate, template_DNA)

    print("fasta-bestand lezen...")


    #use AA sequence as template motif
    #Create temporarily fasta file with "-" removed
    invoer = tempFasta(options.outdirname, options.indirname)
    #Collect all sequences from the temp file, list
    sequences = collectSequences(invoer)


    uitvoerbestand = open(outputfile, 'w')
    uitvoerbestandnoframe = open(outputfilenoframe, 'w')
    counttotal, count1, count2, count3, countmin1, countmin2, countmin3, countnonefound = 0, 0, 0, 0, 0, 0, 0, 0


    nomatch_list = []
    for seq in sequences:
        header, seqline = seq
        seqline = seqline.rstrip("\n")

        totalmatched = count1+count2+count3+countmin1+countmin2+countmin3
        print("+1: {:6} +2: {:6} +3: {:6} -1: {:6} -2: {:6} -3: {:6}  Total matched: {:6}      No matching reading frame: {}     Total: {:6}   \n".format(count1, count2, count3, countmin1, countmin2, countmin3, totalmatched, countnonefound, counttotal), end = "\r")

        counttotal += 1
        if motieven(code.eiwit(seqline[:]), AAtemplate):
            count1 += 1
            uitvoerbestand.write(">"+header +"\n"+ seqline + "\n")
        elif motieven(code.eiwit(seqline[1:]), AAtemplate):
            count2 += 1
            uitvoerbestand.write(">"+header +"\n"+ seqline[1:] + "\n")
        elif motieven(code.eiwit(seqline[2:]), AAtemplate):
            count3 += 1
            uitvoerbestand.write(">"+header +"\n"+ seqline[2:] + "\n")
        elif motieven(code.eiwit(reverseComplement(seqline[:])), AAtemplate):
            countmin1 += 1
            uitvoerbestand.write(">"+header +"\n"+ reverseComplement(seqline[:]) + "\n")
        elif motieven(code.eiwit(reverseComplement(seqline[:-1])), AAtemplate):
            countmin2 += 2
            uitvoerbestand.write(">"+header +"\n"+ reverseComplement(seqline[:-1]) + "\n")
        elif motieven(code.eiwit(reverseComplement(seqline[:-2])), AAtemplate):
            countmin3 += 3
            uitvoerbestand.write(">"+header +"\n"+ reverseComplement(seqline[:-2]) + "\n")
        else:
            print("None found")
            nomatch_list.append(header)
            uitvoerbestandnoframe.write(">"+header +"\n"+ seqline + "\n")
            countnonefound += 1
    totalmatched = count1+count2+count3+countmin1+countmin2+countmin3
    print("+1: {:6} +2: {:6} +3: {:6} -1: {:6} -2: {:6} -3: {:6}  Total matched: {:6}      No matching reading frame: {}     Total: {:6}   \n".format(count1, count2, count3, countmin1, countmin2, countmin3, totalmatched, countnonefound, counttotal), end = "\r")

    uitvoerbestandnoframe.close()
    uitvoerbestand.close()

    print("Finished! Wrote output to {}".format(outputfile))
    print("No match found in: ", nomatch_list)

    bestandnframe = open(outputfile, 'r')
    nframe = 0
    for line in bestandnframe:
        if line.startswith(">"):
            nframe += 1
    bestandnframe.close()

    bestandoutframe = open(outputfilenoframe, 'r')
    outframe = 0
    for line in bestandoutframe:
        if line.startswith(">"):
            outframe += 1
    bestandoutframe.close()


    print("Number of sequences written to 'in frame' file: {}/{}".format(nframe, len(sequences)))
    print("Number of sequences written to 'no correct frame found' file: {}/{}".format(outframe, len(sequences)))
