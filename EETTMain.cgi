#!/usr/bin/python -u

'''
Authors: Mong H. Ng, Lynn Bi, Alec Walls
Mentor: Dr. Sinu Paul, Dr. Sergei Pond, Dr. Perter Bjeorn, Dr. Ward Fleri
Institute: La Jolla Institute of Allergy and Immunology
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
import os
import subprocess
import mysql.connector
from libxml2 import predefinedEntity
from datetime import date, datetime
import random
from lxml.html.builder import AREA
Entrez.email = "lynn.bi1225@gmail.com"
cgitb.enable()

print 'Content-type: text/html'
print '''
<html>
<head>
    <title>EETT | Results</title>
    <link rel="stylesheet" type="text/css" href="/eett.css" />
</head>
<body>
<h2>Ebolavirus Epitope Tracking Tool - Results:</h2><br>
'''
random.seed()
def isValidDate(inputDate):
    '''This is function is needed to check if user's input date is a correctly
    formatted and logical date. No future date.'''
    values = inputDate.split("/")
    #if input is not date at all
    if not len(values)==3:
        return False
    try:
        year = int(values[0])
        month = int(values[1])
        day = int(values[2])
        inputDate = date(year, month, day)
    except:
        return False
    today = date.today()
    if inputDate > today:
        return False
    else:
        return True

def validateInput(form):
    '''This function validate the user's input and return the appropriate message
    back. Currently, organism name cannot be ommited, but included just for future
     change'''
    
    if "organismname" not in form:
        return "Please specify an organism name"
    elif "cutoffdate" not in form or not isValidDate(form["cutoffdate"].value):
        return "Please specify a valid cutoff date"
    elif form["mhcclass"].value=="mhci" and "classialleles" not in form:
        return "Please specify the MHC I alleles you want to predict with"
    elif form["mhcclass"].value=="mhcii" and "classiialleles" not in form:
        return "Please specify the MHC II alleles you want to predict with"
    elif form["mhcclass"].value=="mhci" and "classilength9" not in form and "classilength10" not in form:
        return "Please specify peptide length for MHC I"
    #i do not know why Dr. Sinu had length 15 as checkbox, 
    #since it is the only choice, but i check it anyway
    elif form["mhcclass"].value=="mhcii" and "classiilength15" not in form:
        return "Please specify peptide length for MHC II"
    else:
        return "Valid"

class Input():
    '''This is an encapsulation of the form the user submitted. The str()
    is overriden for debugging purposes'''   
    def __init__(self, form):
        '''The value are string type'''
        self.orgn = form["organismname"].value
        self.pdat = form["cutoffdate"].value
        mhcClass = form["mhcclass"].value
        if(mhcClass=="mhci"):
            self.useMHCI = True
            allelesList = form.getlist("classialleles")
            self.alleles = ",".join(allelesList)
            #if the checkbox is not checked, it is not even present in form
            self.length = list()
            if "classilength9" in form:
                self.length.append("9")
            if "classilength10" in form:
                self.length.append("10")
            self.length = ",".join(self.length)
                
        else:
            self.useMHCI = False
            allelesList = form.getlist("classiialleles")
            self.alleles = ",".join(allelesList)
            self.length = "15"
        self.predictionMethod = form["predictionMethod"].value
    #This function is for debugging purposes
    def __str__(self):
        return '''
        Organism name: %s\n<br/>
        Cutoff date: %s\n<br/>
        Using MHC I: %s\n<br/>
        Alleles: %s\n<br/>
        Peptide length: %s\n<br/>
        Prediction method: %s\n<br/>''' % (self.orgn, self.pdat, self.useMHCI, 
                                    self.alleles, self.length, self.predictionMethod)


def limitResult(idList, num):
    '''This function is for debugging purposes. It limits the records retrieved 
    from NCBI by shortening the UID list input to len 'num'. '''
    l = list()
    for i in range(0,num):
        if i >= len(idList):
            break
        l.append(idList[i])
    return l

def generateFileName():
    '''This function generate a unique file name so that each user can access his
    own unique result files.'''
    return str(datetime.today())+'.'+str(random.randint(0, 999))
    

def queryNCBI(database, term, debugging=False, debuggingSeqLimit=0):
    '''This function queries NCBI and returns a list of SeqRecord. If debugging is 
    set to true, then idlist will be used for limiting result by debuggingSeqLimit,
    thus not query thousand of records. '''
    debuggingSeqLimit = random.randint(3, 8)
    if debugging:
        handleSearch = Entrez.esearch(db=database, term=term)
        recordListSearch = Entrez.read(handleSearch)
        handleSearch.close()
        idList = recordListSearch["IdList"]
        idList = limitResult(idList, debuggingSeqLimit)
        handleFetch = Entrez.efetch(db="protein", rettype="fasta", retmode="text",id=idList)
    else:
        handleSearch = Entrez.esearch(db=database, term=term, usehistory='y')
        recordListSearch = Entrez.read(handleSearch)
        handleSearch.close()
        query_key = recordListSearch["QueryKey"]
        web_env = recordListSearch["WebEnv"]
        handleFetch = Entrez.efetch(db="protein", rettype="fasta", retmode="text", query_key=query_key, webenv=web_env)
    
    seqRecordList = SeqIO.parse(handleFetch, "fasta")
    seqRecordList = list(seqRecordList)
    return seqRecordList

def UniqueSeqFromSeqL1(seqL1, seqL2):
    '''This function takes in two lists of SeqRecord and returns those SeqRecords
    that are only present in seqL1. It is needed to eliminate repeated old sequences
    from the new sequences list'''
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

def queryIEDB():
    '''This function queries iedb and retrieves the epitope sequences(for 
    elobavirus for now). It currently uses intern Lynn Bi's password, which might
    requires another password when the internship ends. This function return a list
     of sequenceRecord of epitope'''
    config = {
    'user': 'lynnb',
    'password': 'Nqg+x3xh2r+UJHfV',
    'host': 'iedb-mysql.liai.org',
    'database': 'iedb_query',
    'port': '33306'
    }
    cnx = mysql.connector.connect(**config)
    cursor = cnx.cursor()
    query = '''select epitope_list_id, source_organism_name, linear_sequence 
    from epitope_list 
    where source_organism_name like "%ebola%" 
    group by linear_sequence limit 20;'''
    cursor.execute(query)  
    epitope_list = list()
    for epitope in cursor:
        id = epitope[0]
        name = epitope[1]
        seq = epitope[2]
        if seq==None:
            continue
        seq = Seq(seq)
        seqR = SeqRecord(seq)
        seqR.description = name
        seqR.id = id
        epitope_list.append(seqR)
    cursor.close()
    cnx.close()
    return epitope_list

def getConservedEpitope(seqL, epitopeL):
    '''return a list of conserved epitope by finding if the string of the epitope 
    sequences is in the string of input sequences'''
    returnL = list()
    for epitope in epitopeL:
        eSeq = str(epitope.seq)
        for seqRecord in seqL:
            seq = str(seqRecord.seq)
            if eSeq in seq:
                returnL.append(epitope)
                break
    return returnL


def mhci_predict(method, inSeq, inAllele, length):
    '''returns a list of mhc I predicted epitope on inSeq. Put into seqRecord so that 
    other functions such as uniqueSeqFromSeqL1 can use it too. The threshold to 
    predict epitope is 1% for mhc I. percentile[7] is the percentile rank; 
    percentile[5] is the sequence of epitope'''
    eList = list()
    command = 'curl --data "method='+method+'&sequence_text='+inSeq+'&allele='+inAllele+'&length='+length+'" http://tools-api.iedb.org/tools_api/mhci/'
    result = subprocess.check_output(command, shell=True)
    rank = result.splitlines( )
    for line in rank:
        percentile= line.split("\t")
        #this line is needed because someimtes "<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">" is printed, idk why
        if len(percentile) <= 7 or percentile[7]== "percentile_rank":
            continue
        elif float(percentile[7]) < 1.0:
            seq = Seq(percentile[5])
            seqRecord = SeqRecord(seq)
            eList.append(seqRecord)
        else:
            #over 1% has been reached and the rest will not be epitope
            break
    return eList

def mhcii_predict(method, inSeq, inAllele):
    '''returns a list of mhc II predicted epitope on inSeq. Put into seqRecord so 
    that other functions such as uniqueSeqFromSeqL1 can use it too. The threshold to 
    predict epitope is 1% for mhc I. percentile[6] is the percentile rank; 
    percentile[4] is the sequence of epitope'''
    eList = list()
    command = 'curl --data "method='+method+'&sequence_text='+inSeq+'&allele='+inAllele+'" http://tools-api.iedb.org/tools_api/mhcii/'
    result = subprocess.check_output(command, shell=True)
    rank = result.splitlines( )
    for line in rank:
        percentile= line.split("\t")
        #this line is needed because someimtes "<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">" is printed, idk why
        if len(percentile) <= 6 or percentile[6]== "percentile_rank":
            continue
        elif float(percentile[6]) < 10.0:
            seq = Seq(percentile[4])
            seqRecord = SeqRecord(seq)
            eList.append(seqRecord)
        else:
            #over 10% has been reached and the rest will not be epitope
            break
    return eList

def deleteRedundantSeq(seqList):
    '''returns a list of seqRecord that does not repeat''' 
    '''seqRecord object are the sequences in the list'''
    returnL = list()
    for a in seqList:
        aSeqString = str(a.seq)
        for b in returnL:
            bSeqString = str(b.seq)
            if aSeqString == bSeqString:
                break
        else:
            returnL.append(a)
    return returnL
#############################################################################################
# taking user's input specification from html form
form = cgi.FieldStorage()
message = validateInput(form)
if not message=="Valid":
    print message
    sys.exit()
input = Input(form)
#print back the input so the user know(at least I myself sometimes forgot what i input)
print '<table class="entry">'
print "<tr><th colspan='2'>User's Input</th></tr>"
print '<tr><td>Organism</td><td>%s</td></tr>'%( input.orgn)
print '<tr><td>Cutoff Date</td><td>%s</td></tr>'%(input.pdat)
if input.useMHCI:
    print '<tr><td>MHC Class</td><td>%s</td></tr>'%('MHC I')
else:
    print '<tr><td>MHC Class</td><td>%s</td></tr>'%('MHC II')
print '<tr><td>Alleles</td><td>%s</td></tr>'%(input.alleles)
print '<tr><td>Peptide length</td><td>%s</td></tr>'%(input.length)
print '<tr><td>Prediction method</td><td>%s</td></tr>'%(input.predictionMethod)
######################################################################################
term = input.orgn+"[ORGN] AND (0[PDAT]:"+input.pdat+"[PDAT])"
fileNameOldSeq = 'new%s.fasta' % generateFileName()
oldSeqRecordList = queryNCBI("protein", term, debugging=True)
oldSeqRecordList = deleteRedundantSeq(oldSeqRecordList)

SeqIO.write(oldSeqRecordList, "/var/www/html/Sequences File/"+fileNameOldSeq, "fasta")
# term = "ebolavirus[ORGN] NOT 2015[PDAT]"
print "<tr><th colspan='2'>Old sequences<br/>(sequences in Genbank before the user's requested date)</th></tr>"
print "<tr><td>GenBank search query</td><td>", term, "</td></tr>"
print "<tr><td>GenBank query</td><td> OK</td></tr>"
print "<tr><td>No. of old sequences (sequences in GenBank before the user's requested date)</td><td>", len(oldSeqRecordList)
print ' <a href="../Sequences File/'+fileNameOldSeq+'">Download</a></td></tr>'
#################################################################################################################################
#collect all Ebola sequences published in the user's specified date range by quering GenBank
term = input.orgn+"[ORGN] AND ("+input.pdat+"[PDAT]:3000[PDAT])"
fileNameNewSeq = 'old%s.fasta' % generateFileName()

newSeqRecordList = queryNCBI("protein", term, debugging=True)
newSeqRecordList = deleteRedundantSeq(newSeqRecordList)
repeatingNewSeqRecordListLen = len(newSeqRecordList)
newSeqRecordList = UniqueSeqFromSeqL1(newSeqRecordList, oldSeqRecordList)

SeqIO.write(newSeqRecordList, "/var/www/html/Sequences File/"+fileNameNewSeq, "fasta")

print "<tr><th colspan='2'>New sequences<br/>sequences in Genbank after the user's requested date)</th></tr>"
print "<tr><td>GenBank search query</td><td>", term, "</td></tr>"
print "<tr><td>GenBank query</td><td> OK</td></tr>"
print "<tr><td>No. of new sequences (Might repeat old sequences)</td><td>", repeatingNewSeqRecordListLen, '</td></tr>'
print "<tr><td>No. of unqie new sequences (Does not repeat old sequences)</td><td>", len(newSeqRecordList)
print ' <a href="../Sequences File/'+fileNameNewSeq+'">Download</a></td></tr>'
#################################################################################################################################
epitopeList = queryIEDB()
epitopeList = deleteRedundantSeq(epitopeList)
conservedEpitopeList = getConservedEpitope(newSeqRecordList, epitopeList)
conservedEpitopeList = deleteRedundantSeq(conservedEpitopeList)
predictedEpitopeList = list()
for seqRecord in newSeqRecordList:
    if input.useMHCI:
        l = mhci_predict(input.predictionMethod, str(seqRecord.seq), input.alleles, input.length)
    else:
        l = mhcii_predict(input.predictionMethod, str(seqRecord.seq), input.alleles)
    predictedEpitopeList = predictedEpitopeList + l
predictedEpitopeList = deleteRedundantSeq(predictedEpitopeList)
predictedEpitopeList = UniqueSeqFromSeqL1(predictedEpitopeList, conservedEpitopeList)
print '<tr><th colspan="2">Epitope</th><tr>'
print '<tr><td>No. of known epitopes from IEDB</td><td>', len(epitopeList),'</td><tr>'
print '<tr><td>No. of conserved epitopes on new sequences</td><td>', len(conservedEpitopeList),'</td></tr>'
print '<tr><td>No. of predicted epitopes on new sequences</td><td>', len(predictedEpitopeList), '</td></tr>'
print '</table>'
print '<h3>The list of conserved epitopes</h3>'
for e in conservedEpitopeList:
    print str(e.seq)
    print "<br/>"
print "<br/>"
print '<h3>The list of predicted epitopes</h3>'
for e in predictedEpitopeList:
    print str(e.seq)
    print "<br/>"
print "<br/>"
print "<br/>"
print '''
</body>
</html>
'''




'''
Note:
the limit result num limit the max of record to fetch, 
but ncbi already limit the return idList to be 20,
unless retmax argument is set.

For debugging purpose, new seq ncbi query and old seq ncbi query and iedb query are
limited.'''
