'''
    Script to identify small ORFs in intergenic regions of the E. coli genome (v3)
    using ribosome profiling data at GSE123675 and GSE122129.    
    
    Copyright (C) 2019  Allen Buskirk, buskirk@jhmi.edu
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
'''


from BCBio import GFF
import csv


# Function to write dict out to csv format.
def writelisttoexcel(genelist,filestring):
	writer = csv.writer(open(filestring+".csv", "wb"),delimiter=',')	
	for gene in genelist:
		writer.writerow(gene)
		
		
def readwig(wigfile):
	resultlist = []
	f = open(wigfile + '.wig')
	for line in f:
	    if line[0] in ['t', 'f']:
	        continue
	    else:
	        resultlist.append(float(line))
	f.close()
	return resultlist


# this function makes a list of strings with an entry for every nt in the chromosome
# if there is a feature annotated on that strand, it gives the gene name (e.g. lacZ)
# otherwise it gives a zero
# note that it needs to be run for each sense (plus = 1, minus = -1)
def geneplot(chrom, sense):

	genepositions = []
	chrlength = len(chrom.seq)
	for i in range(0, chrlength):
		genepositions.append('0')

	for feature in chrom.features:
		if feature.strand != sense:
		      continue
		elif feature.type in ('tRNA', 'rRNA'):
		    gene_name = feature.type  
                    start = feature.location.start.position
                    end = feature.location.end.position
                    for i in range(start, end):
                        genepositions[i] = gene_name
		elif feature.sub_features == []:
		    continue		
		elif feature.sub_features[0].type in ['CDS', 'ncRNA', 'tmRNA']:
		    try:
		        gene_name = feature.qualifiers['Name'][0]
		    except:
			gene_name = feature.sub_features[0].type

                    start = feature.location.start.position
                    end = feature.location.end.position
                    for i in range(start, end):
                        genepositions[i] = gene_name

	return genepositions


# these two functions use Onc112 and retapamulin data to find new initiation sites
# it looks for new ORFs with start codons that are not in a canonical ORF (or within 18nt from one)
# that are at least 8 amino acids long
def INIT_scan(genome, genelist, counts_retapa, counts_Onc112, strand):

    outlist = []
    overlap = False
    length = len(genome)
         
    for pos in range(length-21):
        codon = str(genome[pos:pos+3])
        if codon in ('ATG', 'GTG', 'TTG'): 
            for iC in range(pos-18,pos+18):
                if genelist[iC] != '0':
                    overlap = True
            if overlap:
                overlap = False
                continue  
    
            stop = False
            iL = pos + 3
    
            while not stop:
                if iL >= len(genome):
                    break
                nextcodon = genome[iL:iL+3]
                if nextcodon in ('TAA', 'TAG', 'TGA'):
                    stop = True
                iL += 3
                
            if (iL - pos) > 24:  
                reads1 = sum(counts_retapa[pos:pos+18])
                reads2 = sum(counts_Onc112[pos:pos+18])
                ORF = genome[pos:iL].translate()
                ORF = str(ORF)  
                if strand == 'plus':
                    left = pos + 1                 # convert to 1 based list
                    right = iL
                elif strand == 'minus':            # always left to right on chromosome (not start and stop)
                    right = length - pos
                    left = length - iL + 1
                outlist.append([left, right, reads1, reads2, codon, nextcodon, ORF, strand])
    
    return outlist



def FindInitSites():
    GFFgen = GFF.parse('/users/buskirk/documents/profiling/GFF/MG1655/coli3.gff') 		
    chrom = GFFgen.next()
    f_genome = chrom.seq
    r_genome = chrom.seq.reverse_complement()
    
    ORFlist = []
    ORFlist.append(['start','stop','retapa','onc112','codon1st','codonlast','peptide','strand'])
    
    # for the plus strand
    f1 = 'retapa'       #retapamulin data from Shura Mankin
    pathi = "/users/buskirk/documents/profiling/projects/sORFs/wigfiles/"
    density_filestring1 = pathi + f1
    counts_f1 = readwig(density_filestring1 + "_plus")	
    
    f2 = 'onc'         #Onc112 data
    pathi = "/users/buskirk/documents/profiling/projects/sORFs/wigfiles/"
    density_filestring2 = pathi + f2
    counts_f2 = readwig(density_filestring2 + "_plus")	
    
    gp_plus = geneplot(chrom, 1)         # 1 is plus
    
    plusORFs = INIT_scan(f_genome, gp_plus, counts_f1, counts_f2, 'plus')
    ORFlist.extend(plusORFs)
    
    
    # for the minus strand
    counts_f1 = readwig(density_filestring1 + "_minus")	
    counts_f1.reverse()
    
    counts_f2 = readwig(density_filestring2 + "_minus")	
    counts_f2.reverse()
    
    gp_minus = geneplot(chrom, -1)         # -1 is minus
    gp_minus.reverse()
    
    minusORFs = INIT_scan(r_genome, gp_minus, counts_f1, counts_f2, 'minus')
    ORFlist.extend(minusORFs)
    
    writelisttoexcel(ORFlist, '/users/buskirk/documents/profiling/projects/sORFs/init_sites')    
    


# using newly found initiation sites, look for evidence of ribosome density in normal profiling data
def AreInitSitesUsed():
        
    hitlist = []
    
    with open('/users/buskirk/documents/profiling/projects/sORFs/init_sites.csv', 'rU') as csvfile:
            reader = csv.reader(csvfile, delimiter = ',')
            for row in reader:
                hitlist.append(row)
    
    GFFgen = GFF.parse('/users/buskirk/documents/profiling/GFF/MG1655/coli3.gff') 		
    chrom = GFFgen.next()
    
    genes_plus = geneplot(chrom, 1)
    genes_minus = geneplot(chrom, -1)
    
    pathi = "/users/buskirk/documents/profiling/projects/sORFs/wigfiles/"
    filelist = ['r_control', 'o_control']
      
    for fname in filelist:
        density_filestring = pathi + fname
        counts_plus = readwig(density_filestring + "_plus")
        counts_minus = readwig(density_filestring + "_minus")
    
        for iH in hitlist:
            
            if iH[0] == 'start':
                iH.append(fname)
                continue
                    
            start = int(iH[0]) - 1      # back to 0 based list
            stop = int(iH[1]) - 1
            rpm = 0.0
            counter = 0
            skip = False
            
            if iH[7] == 'plus':
                for iL in range(start - 15, stop + 1 + 15):    
                    if iL >= len(genes_plus):
                        skip = True
                        break
                    
                    if genes_plus[iL] != '0':            # no overlap for 15 nt on either side of new sORF allowed
                        rpkm = genes_plus[iL]
                        skip = True
                        break

                if not skip:                   
                    for iP in range (start, stop + 1):    # + 1 to include last nt in stop codon
                        rpm += counts_plus[iP + 15]
                        counter += 1
        
            elif iH[7] == 'minus':
                for iL in range(start - 15, stop + 1 + 15):           # start is left and stop right in chromosome  
                    
                    if genes_minus[iL] != '0':
                        rpkm = genes_minus[iL]
                        skip = True
                        break                     
                          
                if not skip:            
                    for iP in range(start, stop + 1):     
                        rpm += counts_minus[iP - 15]
                        counter += 1
            
            if counter != 0:
                rpkm = rpm * 1000 / counter
            iH.append(rpkm)
    
    writelisttoexcel(hitlist, '/users/buskirk/documents/profiling/projects/sORFs/init_sites_rpkm')    

FindInitSites()
AreInitSitesUsed()