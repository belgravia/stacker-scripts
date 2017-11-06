#!/usr/bin/env python3
"""
ALISON TANG
10/3/2016

USAGE INSTRUCTIONS
Example command:
$ python fastqPrep.py [options] < myfastq.fastq
At most one input and output can be specified. Input can be inferred and
output defaults to PHRED+33. Input is from stdin, the new FASTQ is in stdout
and error messages are in stderr.

options:
-P33in (--PHRED33input)
-P64in (--PHRED64input)
-P64Bin (--PHRED64Binput with B offset in quality values) 
-P64SOLin  (--PHRED64SOLinput with SOLEXA interpretation of Q score)
-P33out (--PHRED33output) [default]
-P64out (--PHRED64output)
--verbose
--fasta-out
--require-unique-ids

Explanation of the additional options:
A FASTA file can be outputted instead of a FASTQ and is mutually exclusive
with the FASTQ output options. The verbose option prints out certain anomalous
sequences for further inspection. Require unique IDs performs the conversion,
cataloguing number of unique and repeat entries, then outputting a FASTQ with
repeat ID entries discarded.

THE POINT OF THIS PROGRAM
This program can be used to convert between 4 different types of FASTQ files,
Sanger/Illumina 1.8+ (Phred+33), Solexa (Solexa+64), Illumina 1.3+ (Phred+64),
and Illumina 1.5+ (Phred+64 with scores below B unused).
It takes a FASTQ file as input and outputs another FASTQ file either in
Phred+33 or Phred+64 format.
It's useful considering all 4 filetypes exist in the bioinformatic space, so
this program can ease the exchange of fastq data across labs.
"""

import sys
import argparse
import math
import re
class fastqConverter:
    """ The class handles command line arguments in __init__ and the conversion
    of FASTQ files between various quality score formats with the readConvert
    method, which calls one of 6 conversion functions and prints the converted
    quality data, along with the rest of the FASTQ file, to stdout.
    The conversion functions are named such that .convert*_*() converts the
    quality scores from the format of the first wildcard to that of the second.
    """

    def __init__(self, fle):
        """ Parses the command line for the options provided by the user.
        If no input given, then the inferred filetype will be used as the
        conversion to use. After the input is determined/parsed, the general
        use self.convert is set to the relevant conversion method.
        Help with argparse (specifically, groupIN and groupOUT) from Akshar
        Lohith.
        stdin.seek(0) method used to reset stdin back to first byte/line for
        conversion, courtesy of Jon Akutagawa.
        The idea of etting a general method to a more specific method such as
        on lines 69, 116, 118, etc. is as shown in class 9/28/16.
        """
        self.fle = fle
        self.anomalies = self.numRecords = self.inferredType = self.verbose = \
                        self.sequenceLen = self.numNonUnique = 0
        self.processRead = self.processReadStandard
        self.setID = set()
        
        parser = argparse.ArgumentParser(description='process i/o format \
            specification', usage='0-1 input formats can be specified. \
            PHRED+33 output is default')
        groupIN = parser.add_mutually_exclusive_group()
        groupOUT = parser.add_mutually_exclusive_group()
        groupIN.add_argument('-P33in', '--PHRED33input', action='store_true', \
                            dest='PHRED+33', default=False)
        groupIN.add_argument('-P64in', '--PHRED64input', action='store_true', \
                            dest='PHRED+64', default=False)
        groupIN.add_argument('-P64Bin', '--PHRED64Binput', \
                            action='store_true', dest='PHRED+64B', \
                            default=False)
        groupIN.add_argument('-P64SOLin', '--PHRED64SOLinput', \
                            action='store_true', dest='SOLEXA', default=False)
        groupOUT.add_argument('-P33out', '--PHRED33output', \
                            action='store_true', dest='p33o', default=False)
        groupOUT.add_argument('-P64out', '--PHRED64output', \
                            action='store_true', dest='p64o', default=False)
        groupOUT.add_argument('--fasta-out', action='store_true', dest='fo', \
                            default=False)
        parser.add_argument('--verbose', action='store_true', dest='v', \
                            default=False)
        parser.add_argument('--require-unique-ids', action='store_true', dest='uID', \
                            default=False)
        argDict = vars(parser.parse_args())

        self.inferType()
        self.fle.seek(0)
	print('hi')
	print(self.interredType)

        self.verbose = argDict['v']
        self.uniqueID = argDict['uID']  # boolean: discard nonunique IDs?

        #if not (argDict['PHRED+33'] or argDict['PHRED+64'] or \
        #        argDict['PHRED+64B'] or argDict['SOLEXA']):
            #if self.inferredType:
            #    argDict[self.inferredType] = True
            #else:
            #    sys.stderr.write('ERROR: no input format specified and type cannot be detected.\n')
            #    sys.exit(2)

        if argDict['fo']:
            self.readConvert = self.makeFasta
        elif argDict['PHRED+33']:
            if argDict['p64o']:
                self.convert = self.convert33_64  # PHRED+33 to PHRED+64
            else:  # p33o is default
                self.convert = self.noConversion  # if the user so wishes
        elif argDict['PHRED+64']:
            if argDict['p64o']:
                self.convert = self.noConversion
            else:
                self.convert = self.convert64_33
        elif argDict['SOLEXA']:
            if argDict['p64o']:
                self.convert = self.convert64so_64
            else:
                self.convert = self.convert64so_33
        else:  # PHRED+64B
            self.processRead = self.processReadB
            if argDict['p64o']:
                self.convert = self.convert64b_64
            else:
                self.convert = self.convert64b_33

    def readConvert(self):
        """ Reads the FASTQ file and then performs the relevant conversion.
        Reads lines from stdin (i.e. self.fle) by calling self.readEntry,
        looking 4 lines at a time, then converts quality scores with call to 
        the processRead method, and finally prints out a new fastq read entry
        to stdout.
        """
        line = self.fle.readline()
        while line:
            entry = self.readEntry(line, True)
            if not entry:
                line = self.fle.readline()
                continue
            header, seq, qualHeader, qual = entry
            if self.uniqueID:
                if header in self.setID:
                    sys.stderr.write('WARNING: {} not a unique ID, skipped.\n'\
                                    .format(header))
                    self.numNonUnique += 1
                    line = self.fle.readline()
                    continue
                else:
                    self.setID.add(header)
            seq, qual = self.processRead(seq, qual)
            print(header)
            print(seq)
            print(qualHeader)
            print(qual)
            self.numRecords += 1
            line = self.fle.readline()

    def makeFasta(self):
        """ The alternative to readConvert, should the output file be
        specified to be a FASTA file.
        """
        line = self.fle.readline()
        while line:
            entry = self.readEntry(line)
            if not entry:
                line = self.fle.readline()
                continue
            header, seq, qualHeader, qual = entry
            print('>{}'.format(header))
            print(seq)
            self.numRecords += 1
            line = self.fle.readline()

    def readEntry(self, line, report=False):
        """ Reads 4 lines of FASTQ at a time (i.e. an entry). Returns the 4
        lines if no error is detected. Since this method is used both to infer
        type and to convert, printing to stderr can be restricted to during
        conversion by toggling the report argument.
        """
        while line.isspace():
            line = self.fle.readline()
        if line.startswith("@"):
            splitLine = line.rstrip().split(' ')
            header = splitLine[0]  # comment removal
            if len(header) < 15 and len(splitLine) > 1:  # arbitrary; .split on ' ' in illumina1.8 file isn't sufficient for unique IDs
                header = ' '.join(splitLine[:2])
            seq = self.replaceUnknown(self.fle.readline().rstrip()).upper()
            qualHeader = self.fle.readline().rstrip()
            qual = self.fle.readline().rstrip()  # returns None if EOF
            if len(seq) == 0:
                if report:
                    sys.stderr.write('WARNING: Sequence field empty, {} skipped.\n'\
                                    .format(header))
                    self.anomalies += 1  # potentially missing qualHeader and qual
                line = self.fle.readline()
                return 0
            if len(seq) != len(qual):
                if report:
                    sys.stderr.write('WARNING: Sequence and quality lengths differ for {}, skipped.\n'\
                                    .format(header))
                    if self.verbose:
                        sys.stderr.write('Sequence:\n{}\nQuality:\n{}\n'\
                                        .format(seq, qual))
                    self.anomalies += 1  
                line = self.fle.readline()
                return 0
            if not qualHeader.startswith('+'):
                if report:
                    sys.stderr.write('WARNING: Unexpected quality header starting character for {}, skipped.\n'\
                                    .format(header))
                    self.anomalies += 1
                line = self.fle.readline()
                return 0
            if report and self.verbose and len(seq) != self.sequenceLen:
                sys.stderr.write("WARNING: {} sequence length {} differs from the first entry's seq len {}\n"\
                                .format(header, len(seq), self.sequenceLen))
        else:  # improperly formatted ID or the lines are out of order
            sys.stderr.write("ERROR: '@' is not in the first column for {}\n" \
                            .format(line.rstrip()))
            return 0
            
        return header, seq, qualHeader, qual

    def inferType(self):
        """ Infers the FASTQ filetype using the first NUMSEQ number of entries.
        NUMSEQ is a variable determined from somewhat arbitrary assumptions, as
        follows: (1) all sequence entries in the FASTQ are of similar length,
        (2) there is a uniform chance of encountering each quality score, which
        has span of ~45, (3) the probability of an incorrect inferrence given
        (1) and (2) should be less than 10e-100. 
        The equation looks like: 10e-100 = (43/45)^(NUMSEQ*len(line))
        NUMSEQ is calculated to find the probability of not encountering the
        true min and max of the quality score range in the hopes of
        NUMSEQ < total entries in file (i.e. save the step of reading through
        the whole file to determine whether the type is +64B or +64).
        """
        line = self.fle.readline()
        #numSeq = round(math.log10(10e-100)/math.log10(43/45)/len(line))
        numSeq = 1000
	typeCounts = dict.fromkeys(['PHRED+33', 'SOLEXA', 'PHRED+64', \
                                    'PHRED+64B'], 0)
        while line and numSeq:
            entry = self.readEntry(line)
            if not entry:
                line = self.fle.readline()
                continue
            header, seq, qualHeader, qual = entry
            if not self.sequenceLen:  # stores the length of the first sequence
                self.sequenceLen = len(seq)
            qual = [ord(q) for q in list(qual)]
            min_qual, max_qual = min(qual), max(qual)
            if min_qual >= 59 and max_qual <= 73:
                continue  # overlap in ASCII characters used between file types
            elif min_qual >= 33 and max_qual <= 73:
                typeCounts['PHRED+33'] += 1
            elif min_qual >= 66 and max_qual <= 126:  # typically not > 104
                typeCounts['PHRED+64B'] += 1
            elif min_qual >= 64 and max_qual <= 126:
                typeCounts['PHRED+64'] += 1
            elif min_qual >= 59 and max_qual <= 126:
                typeCounts['SOLEXA'] += 1
            else:
                sys.stderr.write("WARNING: Quality scores in range {} to {} for ID {}, type could not be inferred.\n"\
                                .format(min_qual, max_qual, header))
                return 0
            numSeq -= 1
            line = self.fle.readline()

        entriesProcessed = sum(list(typeCounts.values()))
        non64Counts = sum([typeCounts['PHRED+33'], typeCounts['SOLEXA']])
        if typeCounts['SOLEXA'] > 0 and typeCounts['PHRED+33'] == 0:
            self.inferredType = 'SOLEXA'
        elif typeCounts['PHRED+33'] == entriesProcessed:
            self.inferredType = 'PHRED+33'
        elif typeCounts['PHRED+64B'] == entriesProcessed:
            self.inferredType = 'PHRED+64B'
        elif typeCounts['PHRED+64'] > 0 and non64Counts == 0:
            self.inferredType = 'PHRED+64'
        else:
            sys.stderr.write("WARNING: Ambiguous filetype, multiple types detected.\n")

    def processReadStandard(self, seq, qual):
        """Returns the converted quality scores."""
        qual = self.convert(qual)
        return seq, qual

    def processReadB(self, seq, qual):
        """ Returns the converted quality scores, used only for PHRED+64 input
        with a B offset. Differs from the processReadStandard method in that
        the sequence and quality score are coerced to 'N' and 'QO' for bases
        with corresponding quality scores of QO (i.e. 'B'). All quality scores
        are coerced but only the trailing 'B' coerce the sequence to 'N'.

        Regex information from https://docs.python.org/3/library/re.html,
        namely for the '$' help.
        """
        maskN = 0
        trailingB = re.findall(r'B+$', qual)  # finds trailing Bs
        if trailingB:
            maskN = len(trailingB[-1])
        convertedSeq = seq[:(len(seq) - maskN)] + 'N' * maskN  # coerces to 'N'
        qual = qual.replace('B', '@')  # coerces to QO
        convertedQual = self.convert(qual)
        return convertedSeq, convertedQual

    def replaceUnknown(self, seq):
        """ Replace various unknown base characters to N.
        Adapted from http://stackoverflow.com/questions/3411771/multiple- \
        character-replace-with-python
        """
        for ch in ['*', '.']:
            if ch in seq:
                seq = seq.replace(ch, 'N')
        return seq

    def convert64so_33(self, qualSeq):
        """ Performs character by character Solexa to PHRED conversion on as
        seen in Cock et al.:
        http://nar.oxfordjournals.org/content/38/6/1767.full
        """
        convertedQual = [chr(round(10 * math.log10(pow(10, (ord(q) - 64)/10) \
                        + 1)) + 33) for q in list(qualSeq)]
        return ''.join(convertedQual)

    def convert64so_64(self, qualSeq):
        """ Uses the same paper as the conver64so_33 method."""
        convertedQual = [chr(round(10 * math.log10(pow(10, (ord(q) - 64)/10) \
                        + 1)) + 64) for q in list(qualSeq)]
        return ''.join(convertedQual)

    def convert33_64(self, qualSeq):
        """ Both in PHRED, but with an offset difference of 64 - 33 = 31."""
        convertedQual = [chr(ord(q) + 31) for q in list(qualSeq)]
        return ''.join(convertedQual)

    def convert64_33(self, qualSeq):
        convertedQual = [chr(ord(q) - 31) for q in list(qualSeq)]
        return ''.join(convertedQual)

    def convert64b_33(self, qualSeq):
        convertedQual = [chr(ord(q) - 31) for q in list(qualSeq)]
        return ''.join(convertedQual)

    def convert64b_64(self, qualSeq):
        return qualSeq

    def noConversion(self, qualSeq):
        return qualSeq

    def aftermath(self):
        """ Prints warnings to stderr after the whole fastq is traversed
        successfully.
        """
        if self.inferredType:
            sys.stderr.write('Inferred filetype: {}.\n'. \
                            format(self.inferredType))
        if self.numRecords:
            sys.stderr.write('{} records successfully processed.\n' \
                            .format(self.numRecords))
            if self.anomalies > 0:
                sys.stderr.write('WARNING: input has {} anomalous sequences.\n'\
                                .format(self.anomalies))
        else: 
            sys.stderr.write('WARNING: input has no sequences.\n')
        if self.numNonUnique:
            sys.stderr.write('WARNING: Total {} nonunique entries discarded.\n'\
                            .format(self.numNonUnique))
def main():
    converter = fastqConverter(sys.stdin)
    converter.readConvert()
    converter.aftermath()

if __name__ == '__main__':
    main()
