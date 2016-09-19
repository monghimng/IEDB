#!/usr/bin/python
print 'Content-type: text/html'
print '''
<html>
<head>
    <title>EETT | Results</title>
</head>
<body>
<h1>Organism: Ebolavirus<br>
'''
import cgi, cgitb
import cgi
import sys
from cgi import FieldStorage
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import Bio.Align.Applications
from Bio.Align.Applications import ClustalwCommandline
from sys import stdout, stderr
import os
import subprocess
#import mysql.connector
Entrez.email = "lynn.bi1225@gmail.com"
cgitb.enable()

def isValidInput(form):
    if "fdat" not in form:
        return False
    elif "mhci" not in form and "mhcii" not in form:
        return False
    else:
        return True

class Input():     
    def __init__(self, HTMLform):
        self.orgn = "ebolavirus"
        self.fdat = str(HTMLform["fdat"])
        if "tdat" not in HTMLform:
            self.tdat = 3000
        else:
            self.tdat = str(HTMLform["tdat"])
        self.usemhci =  "mchi" in HTMLform
        self.mhciLen = HTMLform["mhciLen"]
        self.usemhcii = "mchii" in HTMLform
        self.outputFormat = HTMLform["outputFormat"]

def UniqueSeqFromSeqL1(seqL1, seqL2):
    '''compare old Ebola sequences before 2015 to Ebola sequences published after 01/01/2015,
    obtain all 'new' sequences
    for our purposes, seqL1 = new list, seqL2 = old list of ebola sequences.'''
    returnL = list()        
    for seqRecord1 in seqL1:
        for seqRecord2 in seqL2:
            seq1 = str(seqRecord1.seq)
            seq2 = str(seqRecord2.seq)
            if seq1 == seq2:
                break
        else:
            #if it never breaks, it is unique and we append it
            returnL.append(seqRecord1)
    return returnL

def epitope_query():
    '''return a list of sequenceRecord of ebola epitope'''
    config = {
    'user': 'lynnb',
    'password': 'Nqg+x3xh2r+UJHfV',
    'host': 'iedb-mysql.liai.org',
    'database': 'iedb_query',
    'port': '33306'
    }
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()
    query = '''select distinct linear_sequence, source_organism_name from epitope_list where source_organism_name like "%ebola%";'''
    cursor.execute(query)  
    epitope_list = list()
    for epitope in cursor:
        name = epitope[1]
        if epitope[0]==None:
            continue
        seq = Seq(epitope[0])
        seqR = SeqRecord(seq)
        seqR.description = name
        epitope_list.append(seqR)
    cursor.close()
    cnx.close()
    return epitope_list

def getConservedEpitope(seqL, epitopeL):
    '''return a list of conserved epitope'''
    returnL = list()
    for epitope in epitopeL:
        eSeq = str(epitope.seq)
        for seqRecord in seqL:
            seq = str(seqRecord.seq)
            if eSeq in seq:
                returnL.addEpitope(epitope)
    return returnL


##########################################################################
#taking user's input specification from html form
form = cgi.FieldStorage()
if not isValidInput(form):
    print "<br>Please enter both the collection date and the prediction type.<br>"
    sys.exit()
input = Input(form)
######################################################################################
#collect all Ebola sequences published before the "from date" by quering GenBank
term = input.orgn+"[ORGN] AND (0[PDAT]:"+input.tdat+"[PDAT])"
term = "ebolavirus[ORGN] NOT 2015[PDAT]"
handlesearch_old = Entrez.esearch(db="protein", term=term, usehistory="y")
recordsearch_old = Entrez.read(handlesearch_old)
handlesearch_old.close()

query_key = recordsearch_old["QueryKey"]
web_env = recordsearch_old["WebEnv"]

handlefetch_old = Entrez.efetch(db="protein", rettype="fasta", retmode="text", query_key=query_key, webenv=web_env)
older_seqlist = SeqIO.parse(handlefetch_old, "fasta")
older_seqlist = list(older_seqlist)

SeqIO.write(older_seqlist, "oldSeqList.fasta", "fasta")

#################################################################################################################################
#collect all Ebola sequences published in the user's specified date ranage by quering GenBank
term = input.orgn+"[ORGN] AND ("+input.fdat+"[PDAT]:"+input.tdat+"[PDAT])"
handlesearch_new = Entrez.esearch(db="protein", term=term, usehistory="y")
recordsearch_new = Entrez.read(handlesearch_new)
handlesearch_new.close()

query_key = recordsearch_new["QueryKey"]
web_env = recordsearch_new["WebEnv"]

handlefetch_new = Entrez.efetch(db="protein", rettype="fasta", retmode="text", query_key=query_key, webenv=web_env)
newer_seqlist = SeqIO.parse(handlefetch_new, "fasta")
newer_seqlist = list(newer_seqlist)

newer_seqlist = UniqueSeqFromSeqL1(newer_seqlist, older_seqlist)
SeqIO.write(newer_seqlist, "newSeqList.fasta", "fasta")
#################################################################################################################################

epitopeL = epitope_query()
conservedEpitopeList = getConservedEpitope(newer_seqlist, epitopeL)
print '''
Date range: from "+input.fdat+" to "+input.tdat+"</h2>"
</body>
</html>
'''