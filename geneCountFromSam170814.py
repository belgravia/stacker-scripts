#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb

#############
# CONSTANTS #
#############

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("-i",
                          dest="samp2sam",
                          type="string",
                          help="Sample name to SAM file",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_file",
                          type="string",
                          help="Output file matrix of gene counts",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")


    samp2samfile_file = open(options.samp2sam)
    outfile = open(options.output_file , "w")

    
    samp2samfile = {} 
    samp_order = []
#    for line in samp2samfile_file:   # i just made edits starting here to bypass the manifest file for 1 sample

 #       line = formatLine(line)

  #      samp, samfile = line.split("\t")
   #     samp2samfile[samp] = samfile
    #    samp_order.append(samp)

    
    #samp2samfile_file.close()
    samp2samfile['hi'] = 'hi'
    samp_order = ['hi']
    gene2samp2count = {}
    gene2txtType = {}

    for samp in samp2samfile:
        samfile = samp2samfile_file #open(samp2samfile[samp])

        for line in samfile:
            # remove headers
            if line.startswith("@"):
                continue
    
            line = formatLine(line)
            lineList = line.split("\t")

            samflag = int(lineList[1])
            
            # Only use primary alignments
            if samflag & 0x900 != 0:
                continue

            # If unmapped
            if samflag & 0x4 != 0:
                continue

            geneName = getGeneName(lineList[2])
            # geneName=lineList
            # gene2txtType[geneName] = txt_type   
 
            if geneName in gene2samp2count:
                if samp in gene2samp2count[geneName]:
                    gene2samp2count[geneName][samp]+= 1
                else:
                    gene2samp2count[geneName][samp] = 1
            else:            
                gene2samp2count[geneName] = {samp:1}

        samfile.close()

    
    # Processed all files
    #headerline = "isoform\tgene\ttranscript_type\t%s\n" % ("\t".join(samp_order)) 
    headerline = "isoform\tcounts\n"
    outfile.write(headerline)

    for gene in gene2samp2count:
        #outline = "%s\t%s\t%s" % (gene)
        outline = gene
        for samp in samp_order:
            try:
                outline += "\t%d" % gene2samp2count[gene][samp]
            except: # No gene count
                outline += "\t0"

        outline += "\n"
        outfile.write(outline)

    outfile.close()
                
             
    
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getGeneName(rname):
    return rname
    elems = rname.split("|")
    return elems[0], elems[5], elems[7]
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
