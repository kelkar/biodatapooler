#!/usr/bin/python
import re, os
import csv
import sys, getopt, Bio
sys.path.append('/software/numpy/1.7.1/b1/lib/python2.7/site-packages')
sys.path.remove('/software/numpy/1.7.1/b1/lib/python3.4/site-packages')
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.PAML import codeml 
import collections
from numpy import *
import scipy as sp
from scipy import stats
from Bio.Phylo.PAML.chi2 import cdf_chi2
import os.path

'''
YOGESHWAR D. KELKAR (y.d.kelkar@gmail.com)

This is a data integration pipeline with a lot of subroutines that are designed to 
integrate data for sets of ortholgous proteins accross multiple species. 
Input data is presumed to be in form of files of various formats (fasta, gff3, cuffdiff, 
for now), each containing information for only one species. 
All data files are to be listed in the control file that is supplied to this code. 
User also supplies the address of the working directory and control files for running PAML 
using the control file. In future, paths of all external tools being used here (PAML, 
translatorx, maftools etc.) will be supplied using this control file as well.
Control file has a tab-separated format that needs to be followed strictly. 

Example:

#speciesTag     type    characteristic  file
ALL     numspecies      param   4
ALL     orthomcl        ortholist       /scratch/ykelkar2/muni/EVM/3May2016/tempundone.txt
ALL     pamlctl pairwise        /scratch/ykelkar2/muni/EVM/3May2016/lysozyme_pairwise_4sp.ctl
ALL     pamlctl null    /scratch/ykelkar2/muni/EVM/3May2016/lysozyme_M0_4sp.ctl
ALL     pamlctl freeRatio       /scratch/ykelkar2/muni/EVM/3May2016/lysozyme_M1_4sp.ctl
MUNI    pamlctl neutral /scratch/ykelkar2/muni/EVM/3May2016/lysozyme_MuniBranchNeutral_4sp.ctl
MUNI    pamlctl special /scratch/ykelkar2/muni/EVM/3May2016/lysozyme_MuniBranch_4sp.ctl
ALL     workdir NA      /scratch/ykelkar2/muni/EVM/3May2016/
MRAP    annotation      gff     /scratch/ykelkar2/muni/EVM/MRAP.evm.out.pepp.gff3
MUNI    annotation      gff     /scratch/ykelkar2/muni/EVM/MUNI.evm.out.pepp.gff3
MUNI    cuffdiff        male_female     /scratch/ykelkar2/muni/annotation/muni_cufflinks_male/MUNI.cuffdiff.2/gene_exp.diff
MRAP    cuffdiff        male_female     /scratch/ykelkar2/rap/annotation/MRAP.cuffdiff.2/gene_exp.diff
MUNI    annotation      protein /scratch/ykelkar2/muni/EVM/MUNI.evm.out.pepp.cds.fasta
MRAP    annotation      protein /scratch/ykelkar2/muni/EVM/MRAP.evm.out.pepp.cds.fasta
MUNI    annotation      nucleotide      /scratch/ykelkar2/muni/EVM/MUNI.evm.out.pepp.cds.fasta
MRAP    annotation      nucleotide      /scratch/ykelkar2/muni/EVM/MRAP.evm.out.pepp.cds.fasta
MUNI    annotation      genome  /scratch/ykelkar2/muni/EVM/muni_sspace_gapfilled_13Nov2014.fa
MRAP    annotation      genome  /scratch/ykelkar2/rap/EVM/rapvS-gapf.gapfilled.final.fa



Integration of these
species-specific datasets is driven by a user-provided list of proteins; most likely
this list is of orthologous proteins. The ortholgous protein list is of the format of
OrthoMCL output:

<name of orthologous group>: <species tag>|<protein 1> <species tag>|<protein 2> ...
 
For example

MN123: hg19|abc1 hg19|abc2 mm|abc1 mm|abc3 

The orthology file is listed in control file as well.

This orthology information is loaded into a special class, which has multiple
subroutines of itself, and which connects to other subroutines that together can perform
multiple tasks for the user:

1. align at nucleotide level
2. align at amino acid sequence level
3. obtain pair-wise dN/dS estimates using PAML, from multiple sequence alignment
4. obtain lineage-specific dN/dS estimates
5. perform log-likelihood tests of selection and neutral evolution using PAML



'''

filehash={}
gff3dict={}
intergenomefilehash={}
gff3line = re.compile('^[a-zA-Z0-9]+')
pairline = re.compile('\-')

def main(argv):
    global Controlfile
    
    usage = '~/bin/ParseOrthoMCL_v3_expr.py -c <PAML control file list>'
    try:
        opts, args = getopt.getopt(argv,"hc:",[ "Controlfile="])
    except getopt.GetoptError, err:
        print >> sys.stderr, usage
        print >> sys.stderr, str(err)
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            print >> sys.stderr,  usage
            sys.exit()
        elif o in ("-c", "--Controlfile"):
            Controlfile = a
        else:
            print >> sys.stderr, usage
            raise RuntimeError("unhandled command line option")
            sys.exit(2)
    print >> sys.stderr,'Controlfile file is ', Controlfile
    print >> sys.stderr,'====================================================================='
    print >> sys.stderr,' READING FILE LIST '
    
    # OPEN CONTROL FILE TO CHECK IF ALL INPUT FILES EXIST, AND ARE PROPERLY FORMATTED
    
    with open(Controlfile) as fl:
        for line in fl:
	    if re.match('^#', line) or not re.search('[0-9a-zA-Z]', line):
	        continue
#	    print line
            f=line.split()
	    if not os.path.exists(f[3]) and f[2] != 'param':
	        print >> sys.stderr, 'file does not exist', f[3]
	        sys.exit(2)
		g = re.split('\-',f[0])
        if f[0] not in filehash:
		    filehash[f[0]]={}
        if f[1] not in filehash[f[0]]:
		    filehash[f[0]][f[1]]={}
		    filehash[f[0]][f[1]][f[2]]=f[3]
		    print >> sys.stderr, "adding to filehash: ", f[0], f[1], f[2], ' = ', f[3]
    print >> sys.stderr,'====================================================================='
    print >> sys.stderr,'================ LOADING FILE DATA INTO DICTIONARY =================='
    
    # LOADING ALL LOCATION, SEQUENCE, GENE EXPRESSION, ETC. DATA FOR EACH GENOME

    for genome in filehash.keys():
	 if genome == 'ALL': continue
	 print >> sys.stderr,'-------------------  genome = ', genome, ' -----------------------'
	 print >> sys.stderr,'reading gff file'
	 gff3dict[genome]=read_gff_files(filehash[genome]['annotation']['gff'], genome)
	 print >> sys.stderr,'reading cuffdiff file for male_female'
	 

	 if 'exonerateRYOs' in filehash[genome]: 
		 for source,target in filehash[genome]['exonerateRYOs'].iteritems():
		     print >> sys.stderr,'reading exonerate file for', genome, source, target
		     gff3dict[genome]=read_exonerate_RYOs(filehash[genome]['exonerateRYOs'][source], gff3dict[genome], genome,source)
         if  'frameshiftsGFF' in filehash[genome]:
		 for source,target in filehash[genome]['frameshiftsGFF'].iteritems():
		     print >> sys.stderr,'reading exonerate file for', genome, source, target
		     gff3dict[genome]=read_frameshift_gff3(filehash[genome]['frameshiftsGFF'][source], gff3dict[genome], genome,source)

	 if 'cuffdiff' in filehash[genome]: gff3dict[genome]=read_cuffdiff_fpkm_values(filehash[genome]['cuffdiff']['male_female'], gff3dict[genome], genome)
	 if 'columndata' in filehash[genome]: gff3dict[genome]=read_columns_by_names(filehash[genome]['columndata']['NA'], gff3dict[genome], genome)
	 gff3dict[genome]=read_cufflinks_fpkm_values(filehash[genome]['cufflinks']['male'], gff3dict[genome], 'male', genome)
	 gff3dict[genome]=read_cufflinks_fpkm_values(filehash[genome]['cufflinks']['female'], gff3dict[genome], 'female', genome)
	 print >> sys.stderr,'reading protein sequence file'
	 gff3dict[genome]=add_sequence_objects(filehash[genome]['annotation']['protein'], gff3dict[genome], genome, 'protein')
	 print >> sys.stderr,'reading nucleotide sequence file'
	 gff3dict[genome]=add_sequence_objects(filehash[genome]['annotation']['nucleotide'], gff3dict[genome], genome, 'nucleotide')
	 


# READ EXONERATE OUTPUTS INTO gff3dict 
def read_exonerate_RYOs(efile, edict, genome, targetgenome):
    with open(efile) as f:
	
        linecount=0
        for line in f:
            linecount=linecount+1
            line=line.rstrip('\n')
            fields=re.split('\s',line)
            if fields[0] == 'RYO':continue
	    if not 'exonerate:'+targetgenome in edict[fields[1]]: edict[fields[1]]['exonerate:'+targetgenome]={}
	    (start,end)=beddize_coordinates(fields[9],fields[10])
	    targetregion=[fields[2],start,end,fields[5],fields[16],fields[17]]
            edict[fields[1]]['exonerate:'+targetgenome][fields[2]+" "+start+" "+end+" "+fields[5]+"/"+fields[3]+" "+fields[16]+" "+fields[17]]=[]
    return edict



# READ FRAMESHIFT DATA INTO THE gff3dict
def read_frameshift_gff3(efile, edict, genome, targetgenome):
    with open(efile) as f:
        linecount=0
        for line in f:
            linecount=linecount+1
            line=line.rstrip('\n')
            fields=re.split('\t',line)
            attrs=re.split(' ; ',fields[8])
	    print >> sys.stderr, "attrs = ", attrs[0], attrs[1], attrs
            attrs[3] = attrs[3].strip()
            targetregion=[fields[0],fields[3],fields[4],fields[8],fields[7]]
	    targetstring=" ".join(targetregion)
            for region, info in edict[attrs[3]]['exonerate:'+targetgenome].iteritems():
		rields = re.split(" ",region)
		if rields[0]==fields[0] and fields[3] >= rields[1] and rields[2] <=fields[4]:
		    print >> sys.stderr, 'printting ', " ".join(targetregion)
		    edict[attrs[3]]['exonerate:'+targetgenome][region].append(targetstring)
    return edict

def print_indel_frameshift_data(grpname, proteins_list, missing):
	frameshifts={}
        print >> sys.stderr, 'extracting exonerate info for ', grpname, proteins_list, missing
        for prot in proteins_list:
	    [genome, pname] = re.split('\|',prot)
            if not genome in frameshifts:     ryos[genome]={}


# SENDING EXONERATE DATA TO PRINTING (printhash command)

def print_exonerate_data(grpname, proteins_list, missing):
    ryos={}
    print >> sys.stderr, 'extracting exonerate info for ', grpname, proteins_list, missing
    for prot in proteins_list:
	[genome, pname] = re.split('\|',prot)
	if not genome in ryos:     ryos[genome]={}
	print >> sys.stderr, 'exonerate info is ', genome, pname, missing, gff3dict[genome][pname]
	if not 'exonerate:'+missing in gff3dict[genome][pname]: continue
	for linecount, lval in gff3dict[genome][pname]['exonerate:'+missing].iteritems():
	    ryos[genome][linecount]=  lval
    printhash(ryos, grpname, 'exonerate:'+missing)
    
def beddize_coordinates(startc, endc):
    if endc < startc: 
            temp=endc
	    endc=startc
	    startc=temp
    return (startc,endc)


# RETREIVE NUCLEOTIDE OBJECTS FROM gff3dict, FOR A LIST OF PROTEIN SEQUENCES
# PRESUMABLY PICKING OUT SEQUENCE DATA FOR ORTHOLOGOUS PROTEINS

def retreive_nucl_objs_from_gff3dict(protein_list):
    prots=[]
    print 'prots'
    for prot in protein_list:
        pparts = re.split('\|',prot)
	print >> sys.stderr, 'getting nucleotide entry for', pparts[0],pparts[1]
        prots.append(gff3dict[pparts[0]][pparts[1]]['mRNA']['nucleotide'])
	print >> sys.stderr, 'Done ', pparts[0], pparts[1], 'mRNA','nucleotide'
    return prots

# PERFORM LOG LIKELIHOOD TEST FOR SELECTION VS. NEUTRAL CONDITIONS USING PAML 

def test_special_neutral(dr, marker1, marker2):
    print >> sys.stderr, 'doing test with', marker1 ,' and ' , marker2, 
    D1 = 2*(dr[marker1]['start:lnL'] - dr[marker2]['start:lnL'])
 #   df1 = dr[marker1]['start:numParams'] - dr[marker2]['start:numParams']
    df1 = dr[marker1]['start:numParams'] - (dr[marker1]['start:numParams'] - 1)
    chip = 1 - stats.chi2.cdf(D1, df1)
    print >> sys.stderr,  ' with D1 and df1 being ', D1, df1 , 'chip is ', chip    
    return chip

def read_columns_by_names(efile, edict, genome):
    with open(efile) as f:
        linecount=0
        gene_id_pos=0
        nameshash={}
        for line in f:
            linecount=linecount+1
            line=line.rstrip('\n')
            fields=re.split('\t',line)
            if linecount == 1:
		for idx, h in enumerate(fields):
		    nameshash[idx]=h
		    if h == 'gene_id' : gene_id_pos=idx
		continue
	    for idx, h in enumerate(fields):
	        if idx == gene_id_pos: 
			gene_id = h
			continue
		if not fields[gene_id_pos] in edict: continue
		if not 'columndata' in edict[fields[gene_id_pos]]:edict[fields[gene_id_pos]]['columndata']={}
		
  	        edict[fields[gene_id_pos]]['columndata'][nameshash[idx]]=h

    print >> sys.stderr, 'done adding ',efile,' to edict'
    return edict

# SENDING CASUAL INPUT COLUMN DATA TO PRINTING

def print_columndata(grpname, proteins_list):
    column = {}
    print >>sys.stderr, 'printing columns for', grpname, proteins_list
    for prot in proteins_list:
        [genome, pname] = re.split('\|',prot)
        if not 'columndata' in gff3dict[genome][pname] : continue
        if not genome in column:column[genome]={}
        for ckey, cvalue in gff3dict[genome][pname]['columndata'].iteritems():
	    column[genome][ckey]=  cvalue
    printhash(column, grpname, 'columndata')


# SEND ALIGNED SEQUENCES TO PAML FOR OBTAINING LINEAGE-SPECIFIC AND PAIRWISE dN, dS AND
# dN/dS VALUES 

def get_dn_ds(aligned_phylip, control_file, marker, dr):
    print >> sys.stderr, 'marker = in getdnds', marker
    if not os.path.isfile(aligned_phylip + "." +marker + "." +".phylip"):
        cml = codeml.Codeml(alignment = aligned_phylip, out_file = aligned_phylip + "." +marker + "." +".phylip", \
        working_dir = filehash['ALL']['workdir']['NA'])	
        print >>sys.stderr, 'cml=',cml
        cml.read_ctl_file(control_file)
        cml.get_option("NSsites")
        results = cml.run()
    else:
        results = codeml.read(aligned_phylip + "." +marker + "." +".phylip")	
	print >> sys.stderr, results
    dr2=rprint(results, 'start', marker, dr)
    return dr2    

# THE MAIN PRINTING FUNCTION FOR PAML OUTPUTS

def rprint(d, head, marker, dr):
    for k, v in d.iteritems():
        if isinstance(v, dict):
	    if k == 'branches':
		message=''
	    elif k == 'start':
		message=''
	    elif k == 'parameters':
		message=''
	    elif k == 'tree length':
		message='tree length'
	    elif k == 'codon model':
		message='codon model'
	    elif head == 'codon model':
		message='codon model'
	    elif k == 'NSsites0':
		head=''
	    else:
		type(head)
	        message = str(head) + str(k)
	    rprint(v, message, marker, dr)
	else:
	    if k == 'parameter list':
		words = v.split()
		v=len(words)
		k = 'numParams'
		print >> sys.stderr, 'words = ',words, 'length = ',v
	    if head == '':
	        head='start'
	    elif head == 'startNSsites0':
	        head='start'
            message = str(head) + ":" + str(k)
	    print >> sys.stderr, marker, message
	    if marker not in dr:
		dr[marker]={}
		v=re.sub('\s','__',v)
		message=re.sub('\s','__',message)
	    dr[marker][message]=v
    return dr

# ALIGN PROTEIN SEQUENCES USING MUSCLE

def align_prot_objs(prots, grpname):
    alignmentInputFile= filehash['ALL']['workdir']['NA'] +'/' + grpname + '.prot.fasta'
    alignmentOutputFile= filehash['ALL']['workdir']['NA'] +'/' + grpname + '.prot.align.fasta'
    print >> sys.stderr, alignmentInputFile, "\t", alignmentOutputFile
    output_handle = open(alignmentInputFile, "w")
    SeqIO.write(prots, output_handle, "fasta")
    output_handle.close()   
    cline = MuscleCommandline(input=alignmentInputFile, out=alignmentOutputFile)
    print >> sys.stderr, cline
    cline()
    handleM = open(alignmentOutputFile, "rU")
    alignedProts = []
    for record in SeqIO.parse(handleM, "fasta") :
        alignedProts.append(record.seq)
    handleM.close()
    return alignedProts


# ALIGN NUCLEOTIDE SEQUENCES USING AMINO ACID SEQUENCE MATCHES, WITH THE HELP OF
# TOOL translatorx (USES MUSCLE)

def align_nucl_objs_phylip(prots, grpname):
    print >> sys.stderr, prots
    path=filehash['ALL']['workdir']['NA'] +'/' + grpname
    alignmentInputFile=path + '.nucl.fasta'
    alignmentOutputFile=path + '.nucl.align'
    print >> sys.stderr, "output phy file = ", alignmentOutputFile + ".nt_ali.fasta"+".phy";
    if os.path.isfile(alignmentOutputFile + ".nt_ali.fasta"+".phy"): 
        print >> sys.stderr, "file exists";
	return alignmentOutputFile + ".nt_ali.fasta"+".phy"
    
    if not os.path.isfile(alignmentInputFile):
        output_handle = open(alignmentInputFile, "w")
        SeqIO.write(prots, output_handle, "fasta")
        output_handle.close()   
    translatorx = "perl ~/bin/translatorx -i " + alignmentInputFile + " -o " + alignmentOutputFile
    alignmentOutputFile+='.nt_ali.fasta'
    print >> sys.stderr, translatorx,"\t",  alignmentOutputFile
    if not os.path.isfile(alignmentOutputFile): os.system(translatorx) 
    if not os.path.isfile(alignmentOutputFile + ".1"):
        fin = open(alignmentOutputFile)
        fout = open(alignmentOutputFile + ".1", "wt")
        for line in fin:
            fout.write( re.sub('\|[a-zA-Z0-9\.\-_]+','',line) )
        fin.close()
        fout.close()

    comm = 'perl ReplaceStopCodonsWithGaps_sort.pl -nuc ' + alignmentOutputFile + ".1"
    print >> sys.stderr, comm
    os.system(comm)
    if not os.path.isfile(alignmentOutputFile + ".phy"):
        input_handle = open(alignmentOutputFile + ".1_nostop.fasta", "rU")
	output_handle = open(alignmentOutputFile + ".phy" , "w")
	alignments = AlignIO.parse(input_handle, "fasta")
	print >> sys.stderr, alignments
	AlignIO.write(alignments, output_handle, "phylip-sequential")
	output_handle.close()
	input_handle.close()

    lfiles=[path + '.nucl.align.aaseqs', path + '.nucl.align.aaseqs.fasta', path + '.nucl.align.nt12_ali.fasta', \
    path + '.nucl.align.nt1_ali.fasta', path + '.nucl.align.nt2_ali.fasta', path + '.nucl.align.nt3_ali.fasta', \
    path + '.nucl.align.aa_based_codon_coloured.html', path + '.nucl.align.html']

    return alignmentOutputFile + ".phy"


def retreive_prot_objs(protein_list):
    print >> sys.stderr, protein_list
    prots=[]
    for prot in protein_list:
	prots.append(protsDict[prot].protobj)
    return prots

def read_gff_files(gfffile, genome):
    print >> sys.stderr,'reading gff file ', gfffile
    gff3dict={}
    with open(gfffile) as f:
	counter=0
        for line in f:
#	    print line
	    line=line.rstrip('\n')
	    if not gff3line.match(line):
                continue
            fields=re.split('\t',line)
            if fields[2] == 'mRNA':
               	attrs=re.split('[;=]', fields[8])
		if not attrs[1] in gff3dict: gff3dict[attrs[1]]={}
		if not 'mRNA' in gff3dict[attrs[1]]: gff3dict[attrs[1]]['mRNA']={}
		if not 'CDS' in gff3dict[attrs[1]]: gff3dict[attrs[1]]['CDS']={}
		
	      	gff3dict[attrs[1]]['mRNA']['chrom']=fields[0]
	      	gff3dict[attrs[1]]['mRNA']['start']=fields[3]
	      	gff3dict[attrs[1]]['mRNA']['end']=fields[4]
	      	gff3dict[attrs[1]]['mRNA']['strand']=fields[6]
		
	      	gff3dict[attrs[1]]['mRNA']['male.FPKM']=0
	      	gff3dict[attrs[1]]['mRNA']['female.FPKM']=0
#		print >> sys.stderr, 'adding to gff3dict for ',genome,'::',attrs[1], 'mRNA : maleFPKM '
		counter=0
            elif fields[2] == 'CDS':
		counter+=1
               	attrs=re.split('[;=]', fields[8])
#		print attrs
		maletupp = (fields[3], fields[4], 'malecov')
		femaletupp = (fields[3], fields[4], 'femalecov')
#		print maletupp
                gff3dict[attrs[3]]['CDS'][maletupp]=0
                gff3dict[attrs[3]]['CDS'][femaletupp]=0
#		print 'adding to gff3dict for ',genome,'::',attrs[1], 'CDS', maletupp

    return gff3dict


# ADD SEQUENCE OBJECTS TO gff3dict FROM FASTA FILE

def add_sequence_objects(sfile, sdict, genome, stype):
    handle = open(sfile, "rU")
    pdict=sdict
    print >> sys.stderr,  'reading file', sfile
    for record in SeqIO.parse(handle, "fasta") :
        idsplit=re.split('\|',record.id)
        if len(idsplit) == 2:
            if idsplit[0] != genome:
	        print >> sys.stderr, "in function add_sequence_objects, the sequence file does not match genome"
            else:
		if idsplit[1] in pdict:
	            if not record.id in pdict: pdict[record.id]={}
	            if not record.id in pdict[idsplit[1]]: pdict[idsplit[1]]['mRNA']={}
                    pdict[idsplit[1]]['mRNA'][stype]=record
	else:
	    if not record.id in pdict: pdict[record.id]={}
	    if not 'mRNA' in pdict[record.id]: pdict[record.id]['mRNA']={}
	    idcopy=record.id
	    record.id=genome+"|"+idcopy
	    pdict[idcopy]['mRNA'][stype]=record
            if idcopy == "evm_model_ellu5Jan2014_scaffold717_size71769_2" :print >> sys.stderr, idcopy, record.id
    handle.close()
    print  >> sys.stderr, "DONE adding"
    return pdict

# READ CUFFLINKS VALUES

def read_cufflinks_fpkm_values(efile, edict, source, genome):
    print  >> sys.stderr,'reading cufflinks file',efile,'into dict, which is of size', len(edict)
    with open(efile) as f:
        for line in f:
	    if not gff3line.match(line):
                continue
            fields=re.split('\t',line)
	    fields[8]=re.sub('\s*;\s*',';',fields[8])
	    fields[8]=re.sub(';$','',fields[8])
            attrs =re.split(';',fields[8])
            names={}
	    for attr in attrs:
	        n = re.split('\s+',attr)
		while 1:
		    n = [re.sub('[\s"]+','',s) for s in n]
                    if not re.match('\"', n[0]):
			break
	        names[n[0]]=n[1]
            if fields[2] == 'transcript':
		fpkmattr=source+'.FPKM'
	        edict[names['reference_id']]['mRNA'][fpkmattr]=names['FPKM']
            if fields[2] == 'exon':
		tupp = (fields[3], fields[4], source+'cov')
	        if tupp in edict[names['reference_id']]['CDS']:
		    edict[names['reference_id']]['CDS'][tupp]=names['cov']
    print >> sys.stderr, 'done adding ',efile,' to edict'
    return edict

# READ CUFFDIFF VALUES

def read_cuffdiff_fpkm_values(efile, edict, genome):
    print  >> sys.stderr,'reading cufflinks file',efile,'into dict, which is of size', len(edict)
    with open(efile) as f:
        for line in f:
            fields=re.split('\t',line)
            if fields[0] == 'test_id': continue
	    edict[fields[0]]['mRNA']['male.FPKM']=fields[7]
	    edict[fields[0]]['mRNA']['female.FPKM']=fields[8]
    print >> sys.stderr, 'done adding ',efile,' to edict'
    return edict

def isfloat(value):
  try:
    float(value)
    return True
  except:
    return False

# PRINTING LOT OF EXPRESSION RELATED STATISTICS

def print_expression(grpname, proteins_list):
    expr = {}
    malelist=[]
    femalelist=[]
    foldchangelist=[]
    for prot in proteins_list:
	[genome, pname] = re.split('\|',prot)
	print >> sys.stderr, 'trying to extract maleFPKM using ',genome, pname, 'mRNA: maleFPKM'
	mval=float(gff3dict[genome][pname]['mRNA']['male.FPKM'])
	fval=float(gff3dict[genome][pname]['mRNA']['female.FPKM'])
	print >> sys.stderr, 'is float mval and fval', isfloat(mval), isfloat(fval)
	foldchange=0
        if not genome in expr: expr[genome]={}	
        if not 'all' in expr: expr['all'] = {}
	expr[genome]['male.FPKM']=mval
	expr[genome]['female.FPKM']=fval
	expr[genome]['proteinName']=pname
	expr[genome]['foldchange']=0
	expr[genome]['logFC']='NA'
	if mval + fval != 0:            
	    if mval > 0 and fval > 0: foldchange = mval/fval
	    if mval > 0 and fval == 0: foldchange = mval/1 	    
	    if mval == 0 and fval > 0: foldchange = 1/fval
            print >> sys.stderr,'mval, fval and foldchange =' , mval, fval, foldchange	
	    expr[genome]['foldchange']=foldchange
	    expr[genome]['logFC']=log(foldchange)
	malelist.append(mval)
	femalelist.append(fval)
	foldchangelist.append(foldchange)
    print >> sys.stderr, 'male and female and foldchange lists are ', malelist, femalelist, foldchangelist
	
    expr['all']['foldchangemean']=  mean(foldchangelist)
    expr['all']['foldchangemin']=  min(foldchangelist)
    expr['all']['foldchangemax']=  max(foldchangelist)
    expr['all']['malemean']=  mean(malelist)
    expr['all']['malemin']=  min(malelist)
    expr['all']['malemax']=  max(malelist)
    expr['all']['femalemean']= mean(femalelist)
    expr['all']['femalemin']=  min(femalelist)
    expr['all']['femalemax']=  max(femalelist)
    printhash(expr, grpname, 'expression')

def print_missing_genome_section(missing_species, proteins_list):
    for prot in proteins_list:
	[genome, pname] = re.split('\|',prot)
	gene_start = gff3dict[genome][pname]['mRNA']['start']
	gene_end = gff3dict[genome][pname]['mRNA']['end']
	gene_chrom = gff3dict[genome][pname]['mRNA']['chrom']
	align_results=gff3dict[genome]['idx'].search()


# MAIN PRING SUBROUTINE

def printhash(d, filen, type):
    #print len(d)
    for k, v in d.iteritems():
        for k2, v2 in v.iteritems():
	    print >> sys.stderr, 'OUTdata', filen, k, k2, v2
	    outlist = ['data',filen, type, k, k2, v2]
            print >> sys.stdout, '\t'.join(map(str, outlist))

if __name__ == "__main__":
    main(sys.argv[1:])

class orthogroup(object):

    def __init__(self, name, size, species, num_species, proteins_list):
        self.name = name
        self.size = size
        self.species = species
        self.num_species = num_species
        self.proteins_list = proteins_list

    def prot_align(self):
        prots=retreive_prot_objs(self.proteins_list)
        aligned_prots=align_prot_objs(prots, self.name)
        return aligned_prots    

    def nucl_align_file(self):
        nucs=retreive_nucl_objs_from_gff3dict(self.proteins_list)
        aligned_nucs_file=align_nucl_objs_phylip(nucs, self.name)
        return aligned_nucs_file

    def get_pairwise_dn_ds(self):
        aligned_phylip=self.nucl_align_file()
        dr = collections.OrderedDict()
        dr['pairwise']={}	
        nullctl = 	filehash['ALL']['pamlctl']['pairwise']
	outfile = aligned_phylip + "." +'pairwise' + "." +".phylip" 
        if not os.path.isfile(outfile):
            print >> sys.stderr, 'working_dir = ', filehash['ALL']['workdir']['NA'],'aligned_phylip=', \
            aligned_phylip
            cml = codeml.Codeml(alignment = aligned_phylip, out_file = aligned_phylip + "." +'pairwise' + "." +".phylip", \
            working_dir = filehash['ALL']['workdir']['NA'])
	    cml.read_ctl_file(nullctl)
	    cml.get_option("NSsites")
	    print 'cml=',cml
            results = cml.run()
	else:
	    results = codeml.read(outfile)
	    print >> sys.stderr, results

        paircapture = "cat " + outfile + " | perl -p -e \'s/\n/\t/g\' | grep -oP \"(?<=\t)([0-9]+)\s+\([a-zA-Z0-9]+\)\s+\.\.\.\s+([0-9]+)\s+\([a-zA-Z0-9]+\)\tlnL\s*=\s*[0-9\-\.]+\t\s+[0-9\.\-]+\s+[0-9\.\-]+\t\t([0-9a-zA-Z=\S \.\-]+)\t\" | perl -p -e \'s/\s*=\s*/\t=\t/g\' | perl -p -e \'s/[ \t]+/\t/g\' | cut -f 1,4,22,25,28 > " + outfile+".table"
        os.system(paircapture) 
        with open(outfile+".table") as f:
            for line in f:
		print >> sys.stderr,  line
                line=line.rstrip('\n')
                fields=line.split()
        	dr['pairwise'][fields[0]+".."+fields[1]+":dN/dS"]=fields[2]
        	dr['pairwise'][fields[0]+".."+fields[1]+":dN"]=fields[3]
        	dr['pairwise'][fields[0]+".."+fields[1]+":dS"]=fields[4]
        printhash(dr, self.name, 'paml')

    def test_selection(self):
        aligned_phylip=self.nucl_align_file()
        dr = collections.OrderedDict()  
        nullctl = 	filehash['ALL']['pamlctl']['null']
	dr = get_dn_ds(aligned_phylip, nullctl , 'null', dr)
    	for genome in filehash.keys():
            if genome == 'ALL':
               dr = get_dn_ds(aligned_phylip, filehash[genome]['pamlctl']['null'] , 'null', dr)
               continue
            if not 'pamlctl' in filehash[genome]: continue
            specialctl = filehash[genome]['pamlctl']['special']
            neutralctl = filehash[genome]['pamlctl']['neutral']
            dr = get_dn_ds(aligned_phylip, specialctl , genome+'.special', dr)
            dr = get_dn_ds(aligned_phylip, neutralctl , genome+'.neutral', dr)
            special_p=test_special_neutral(dr, genome+'.special' , 'null')
            neutral_p=test_special_neutral(dr, genome+'.special' , genome+'.neutral')
	    if not genome+'.special' in dr: dr[genome+'.special']={}
	    if not genome+'.neutral' in dr: dr[genome+'.neutral']={}
	    
            dr[genome+'.special']['LRT']=special_p
            dr[genome+'.neutral']['LRT']=neutral_p  	
	
	printhash(dr, self.name, 'paml')

    def collect_expression(self):
	print_expression(self.name, self.proteins_list)
    def collect_columndata(self):
	print_columndata(self.name, self.proteins_list)
    def collect_exonerates(self):
        species_in_set =set(self.species)
        all_species_set =set(gff3dict.keys())
	missing_species = all_species_set - species_in_set
	print_exonerate_data(self.name, self.proteins_list,  list(missing_species)[0] )
    def genome_region_correspondence_from_orthomcl_clusters(self):
        species_in_set =set(self.species)
        all_species_set =set(gff3dict.keys())
	missing_species = all_species_set - species_in_set
	print >> sys.stderr,'speciesinset, allspeciesst and missing species are ', species_in_set, all_species_set, missing_species
	print >> sys.stderr, 'missing species = ', missing_species
	for prot in self.proteins_list:
            [genome,pname]=re.split('\|',prot)
	    print >> sys.stderr, 'looking at protein', genome, pname, gff3dict[genome][pname]['mRNA'].keys()
	    
  	    gene_start = gff3dict[genome][pname]['mRNA']['start']
	    gene_end = gff3dict[genome][pname]['mRNA']['end']
	    gene_chrom = gff3dict[genome][pname]['mRNA']['chrom']
	    genomelc = genome.lower()
            for missing in missing_species:
		if genome in intergenomefilehash[missing]:
	            alignfile = intergenomefilehash[missing][genome]['bialignment']
                    missing_command = '/home/ykelkar2/privatemodules/mafTools/bin/mafExtractor --maf '+ \
		    alignfile+' --seq '+ genomelc+'.'+gene_chrom + ' --start ' + gene_start + ' --stop ' + gene_end + \
		    ' > ' + filehash['ALL']['workdir']['NA'] +'/' + self.name + '_missing_' + missing + '_from_' + genome + '-' + pname + '.maf'
		
		    print  >> sys.stderr, missing_command
		else:
		    continue



groupsDict={}


orthomclfile=filehash['ALL']['orthomcl']['ortholist']
print >> sys.stderr, 'LOADING groupsDict from file', orthomclfile

ungroupsDict={}

with open(orthomclfile) as f:
    for line in f:
#        print >> sys.stderr, line
        line=line.rstrip('\n')
        line = re.sub(':','',line)
        allspecies=re.findall('[A-Z]+\|', line)
        allspecies= [ re.sub("\|$", "", elem) for elem in allspecies ]
        species=list(set(allspecies))
        fields=line.split()
        grname=fields.pop(0)
        size=len(fields)
        groupsDict[grname]=orthogroup(grname, len(fields), list(set(allspecies)), len(list(set(allspecies))), fields)


print >> sys.stderr, 'DONE'

for key in groupsDict:
 #   print >> sys.stderr, key, filehash['ALL']['numspecies']['param'] , groupsDict[key].size,  groupsDict[key].num_species
    print >> sys.stderr, "consdering", key
    print groupsDict[key].nucl_align_file()
    print >> sys.stderr, "examining", key
    print groupsDict[key].get_pairwise_dn_ds()	
    print groupsDict[key].test_selection()
    print groupsDict[key].collect_expression()
    print groupsDict[key].collect_columndata()


print >> sys.stderr, "--------------------Finished run-------------------"

