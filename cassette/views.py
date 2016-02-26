# Django imports 
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from .models import Results
from .forms import *

import os
from django.conf import settings
PROJECT_ROOT = settings.PROJECT_ROOT

# Imports for cassette generation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp
from Bio.Alphabet import SingleLetterAlphabet
import copy
from intermine.webservice import Service
from pandas import *
from pandas import DataFrame, read_csv
import pandas as pd #this is how I usually import pandas
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# define global variables
HomologyLength = 1000
PrimerMaxTm = 55
PrimerMaxLen = 60
OverhangMaxFrac = 1

def index(request):
	# if this is a POST request we need to process the form data
	if request.method == 'POST':
		# create a form instance and populate it with data from the request
		form = TrueIndexForm(request.POST)
		# validate
		choice0 = request.POST['choices0']
		if choice0 == "1":
			cutsite = request.POST['cutsite']
			choice1 = request.POST['choices1']

			if choice1 == "1":
				name = request.POST['name1']
				sequence = request.POST['sequence1']
				ster = request.POST['sequence1']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editEmpty(name, sequence, cutsite)
				

			elif choice1 == "2":
				orfName = request.POST['orfName1']
				orfSeq = request.POST['orfSeq1']
				promoterName = request.POST['promoterName1']
				terminatorName = request.POST['terminatorName1']

				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editEmpty(orfName, orfSeq, cutsite, promoterName, terminatorName)
				
			request.session['Lup'] = Lup
			request.session['Rup'] = Rup
			request.session['Ldown'] = Ldown
			request.session['Rdown'] = Rdown
			request.session['L'] = L
			request.session['R'] = R
			request.session['seqLen'] = seqLen
			request.session['donorSeq'] = donorSeq
			return HttpResponseRedirect('/cassetteBuilder/results')

		elif choice0 == "2":
			locus = request.POST['locus']
			choice2 = request.POST['choices2']

			if choice2 == "1":
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editExisting(locus, 1)
			elif choice2 == "2":
				name2 = request.POST['name2']
				sequence2 = request.POST['sequence2']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editExisting(locus, 2, NewGeneName=name2, NewGeneSeq=sequence2)
			elif choice2 == "3":
				promoterName = request.POST['promoterName2']
				terminatorName = request.POST['terminatorName2']
				orfName = request.POST['orfName2']
				orfSeq = request.POST['orfSeq2']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq = editExisting(locus, 3, promoter=promoterName,terminator=terminatorName,NewGeneName=orfName,NewGeneSeq=orfSeq)

			elif choice2 == "4":
				pass
			elif choice2 == "5":
				pass

			request.session['Lup'] = Lup
			request.session['Rup'] = Rup
			request.session['Ldown'] = Ldown
			request.session['Rdown'] = Rdown
			request.session['L'] = L
			request.session['R'] = R
			request.session['seqLen'] = seqLen
			request.session['donorSeq'] = donorSeq
			return HttpResponseRedirect('/cassetteBuilder/results')

		elif choice0 == "3":
			choice3 = request.POST['choices3']
			if choice3 == "1":
				geneList = request.POST.getlist('geneList1[]')
				seqList = request.POST.getlist('seqList1[]')
				answer = variableCassette(geneList, seqList)

				request.session['customAnswer'] = answer

			elif choice3 == "2":
				geneList = request.POST.getlist('geneList2[]')
				seqList = request.POST.getlist('seqList2[]')
				toVary = request.POST['toVary']
				variants = request.POST.getlist('geneList3[]')
				variantSeq = request.POST.getlist('seqList3[]')
				answer = variableCassette(geneList, seqList, toVary, variants, variantSeq)

				request.session['customAnswer'] = answer

			return HttpResponseRedirect('/cassetteBuilder/customResults')

	# If a GET (or any other method)we'll create a blank form) 
	else:
		form = TrueIndexForm()
	df = pd.read_excel(os.path.join(PROJECT_ROOT, "cutsites.xlsx"))
	table = df.to_html(index=False)
	return render(request, 'cassette/index.html', {'form': form, 'table': table})

def results(request):
	Lup = request.session.get('Lup')
	Rup = request.session.get('Rup')
	Ldown = request.session.get('Ldown')
	Rdown = request.session.get('Rdown')
	L = request.session.get('L')
	R = request.session.get('R')
	seqLen = request.session.get('seqLen')
	donorSeq = request.session.get('donorSeq')

	return render(request, 'cassette/results.html', {'Lup': Lup, 'Rup': Rup, 'Ldown': Ldown, 'Rdown': Rdown, 'L': L, 'R': R, 'seqLen': seqLen, 'donorSeq': donorSeq})

def customResults(request):
	answer = request.session.get('customAnswer')
	original = answer[0][0]
	variants = answer[1:]
	return render(request, 'cassette/customResults.html', {'original' : original, 'variants' : variants})


#---------------  Helper Functions  ---------------------#
def editEmpty(name, sequence, cutname, promoter = None, terminator = None):
	df = pd.read_excel(os.path.join(PROJECT_ROOT, "cutsites.xlsx"))

	labels=df['name'].values
	ChrLetters=df['chrom. loc.'].values
	ExpValues=df['exp. lev.'].values
	cutSeqs=df['sequence'].values

	cutArray={'name' : Series(labels, index=labels),
				'exp. lev.' : Series(ExpValues, index=labels),
				'chrom. loc.' : Series(ChrLetters, index=labels),
				'sequence' : Series(cutSeqs, index=labels)
			 }

	cutFrame= DataFrame(cutArray)
	location=cutFrame.loc[cutname,'chrom. loc.']+".fasta"
	cutSequence=cutFrame.loc[cutname,'sequence']

	ChromosomeSeq=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\\" + location), "fasta").seq
	
	if ChromosomeSeq.find(cutSequence)==-1:
		ChromosomeSeq=ChromosomeSeq.reverse_complement()

	StartIndex=ChromosomeSeq.find(cutSequence)
	EndIndex=StartIndex+34
	
	UpSeq=ChromosomeSeq[StartIndex-HomologyLength:StartIndex]
	DownSeq=ChromosomeSeq[EndIndex:EndIndex+HomologyLength]
		
	UpHomRec = SeqRecord(UpSeq, id=cutname)
	DownHomRec = SeqRecord(DownSeq, id=cutname)

	orfRecord = SeqRecord(Seq(sequence, SingleLetterAlphabet()), id=name)
	if promoter is None: 
		fragments = [UpHomRec, orfRecord, DownHomRec]
	else:
		PromoterGeneRec = fetchGene(promoter)
		PromoterRec = fetchNeighbor(PromoterGeneRec,"upstream",600)
		PromoterRec.id = PromoterRec.id + "ps"

		TerminatorGeneRec = fetchGene(promoter)
		TerminatorRec = fetchNeighbor(TerminatorGeneRec,"upstream",600)
		TerminatorRec.id = TerminatorRec.id + "ts"

		fragments = [UpHomRec, PromoterRec, orfRecord, TerminatorRec, DownHomRec]
	return stitch(fragments)

def editExisting(name, option, promoter = None, terminator = None, NewGeneName = "", NewGeneSeq = ""):
	OrigGeneRecord = fetchGene(name)
	UpHomRec = fetchNeighbor(OrigGeneRecord, "upstream", HomologyLength)
	DownHomRec = fetchNeighbor(OrigGeneRecord, "downstream", HomologyLength)

	if option == 1:
		fragments = [UpHomRec, DownHomRec]
		
	elif option == 2:
		InsertRec = SeqRecord(Seq(NewGeneSeq, SingleLetterAlphabet()), id=NewGeneName)
		fragments = [UpHomRec, InsertRec, DownHomRec]
	elif option == 3:
		PromoterRec, orfRecord, TerminatorRec = standardCassette(promoter, terminator, NewGeneName, NewGeneSeq)
		fragments = [UpHomRec, PromoterRec, orfRecord, TerminatorRec, DownHomRec]
	elif option == 4:
		pass
	elif option == 5:
		pass

	return stitch(fragments)

def fetchGene(GeneName):
	
	service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
	template = service.get_template('Gene_GenomicDNA')

	rows = template.rows(
		E = {"op": "LOOKUP", "value": GeneName, "extra_value": "S. cerevisiae"}
	)
	
	# this service seems to return multiple similar genes but we want the first one only, so count
	# and it returns information about the gene you want
	count=0
	for row in rows:
		
		count=count+1
		if count==1:
			descr= row["description"]
			GeneSeq=Seq(row["sequence.residues"])
			GeneSysName=row["secondaryIdentifier"]
			#print(" ")
			#print("I think you want...... "+row["secondaryIdentifier"])
			#print(row["description"])
			#print(" ")
			#print(row["sequence.residues"])
			#print(" ")
			#print("Good choice! I have a feeling you're going to get lucky with this one.")
			#print(" ")
			#print("Give me a second to put some of my ducks in a circle...")
	   

			
	#let's create a record for the oldGene
	GeneRecord = SeqRecord(GeneSeq, id=GeneSysName)
	
	#now let's add some more information to make it useful
	GeneRecord.name=GeneName
	GeneRecord.features=GeneSysName
	GeneRecord.description=descr

	return GeneRecord 

def fetchNeighbor(NeighborRecord, direction, distance):


	# let's load the appropriate chromosome file. The record of the gene we looked up
	# contains in the "features" the systematic name, wherein the second letter
	# corresponds to chromosome number, e.g., 1=A etc
	if NeighborRecord.features[1]=="A":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer01.fasta"), "fasta")
	if NeighborRecord.features[1]=="B":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer02.fasta"), "fasta")
	if NeighborRecord.features[1]=="C":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer03.fasta"), "fasta")
	if NeighborRecord.features[1]=="D":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer04.fasta"), "fasta")
	if NeighborRecord.features[1]=="E":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer05.fasta"), "fasta")
	if NeighborRecord.features[1]=="F":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer06.fasta"), "fasta")
	if NeighborRecord.features[1]=="G":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer07.fasta"), "fasta")
	if NeighborRecord.features[1]=="H":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer08.fasta"), "fasta")
	if NeighborRecord.features[1]=="I":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer09.fasta"), "fasta")
	if NeighborRecord.features[1]=="J":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer10.fasta"), "fasta")
	if NeighborRecord.features[1]=="K":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer11.fasta"), "fasta")
	if NeighborRecord.features[1]=="L":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer12.fasta"), "fasta")
	if NeighborRecord.features[1]=="M":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer13.fasta"), "fasta")
	if NeighborRecord.features[1]=="N":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer14.fasta"), "fasta")
	if NeighborRecord.features[1]=="O":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer15.fasta"), "fasta")
	if NeighborRecord.features[1]=="P":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes\Scer16.fasta"), "fasta") 

	
	
	# let's explicitely name the sequences from the seq record
	NeighborSeq=NeighborRecord.seq
	ChromosomeSeq=ChromosomeRec.seq
	
	# flip the sequence to orient with respect to the old gene
	if ChromosomeSeq.find(NeighborSeq)==-1:
		ChromosomeSeq=ChromosomeSeq.reverse_complement()

	StartIndex=ChromosomeSeq.find(NeighborSeq)
	EndIndex=StartIndex+len(NeighborSeq)
	
	if direction=="upstream":
		DesiredSeq=ChromosomeSeq[StartIndex-distance:StartIndex]
	if direction=="downstream":
		DesiredSeq=ChromosomeSeq[EndIndex:EndIndex+distance]

	
	
	
	NeighborRec = SeqRecord(DesiredSeq, id=NeighborRecord.name)
	
	return NeighborRec

	#print(NeighborRec)

def getPrimer(currRecord):
	

	mp = 0
	length = 0
	primer = Seq("")

	seq=currRecord.seq
	
	while mp <= PrimerMaxTm and length <= PrimerMaxLen:
		primer = primer + seq[length]
		mp = MeltingTemp.Tm_staluc(primer)
		length += 1

	return primer           
		
def overhangPrimer(currRecord,prevSeq):
	#let's get the template-binding primer first
	primer=getPrimer(currRecord)
	
	
	#OK let's work on the overhang
	maxOhLen=PrimerMaxLen-len(primer)    
	maxFrac=1
	
	#let's decide on a max overhang length
	if round(len(primer)*(OverhangMaxFrac+1)) < 60:
			 maxOhLen=round(len(primer)*OverhangMaxFrac)
	
	#the index must be an integer!!!
	maxOhLen=int(maxOhLen)
	ohprimer=prevSeq.seq[-maxOhLen:]+primer #we add the .seq so that it returns a string
	
	return ohprimer      

def stitch(fragments):
	#this function takes seq records and prints primers

	#let's make an empty sequence file
	Nfrags=len(fragments)
	donor=Seq("")
	index=[]
	print("")
	for i in range (0, Nfrags):
		donor=donor+fragments[i]
	# Dummy assignment setup to allow for compilation
	Lup = ""
	Rup = ""
	Ldown = ""
	Rdown = ""
	L = ""
	R = ""

	for i in range (0, Nfrags):
		if i==0:
			Lup = "Lup"+ fragments[i].id + " " + getPrimer(donor)
			Rup = "Rup"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement())
		elif i==Nfrags-1:
			Ldown = "Ldown"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1])
			Rdown = "Rdown"+ fragments[i].id + " " + getPrimer(donor.reverse_complement())
		else:
			L = "L"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1])
			R = "R"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement())

	sequenceLength = len(donor.seq)
	donorSequence = donor.seq

	return str(Lup), str(Rup), str(Ldown), str(Rdown), str(L), str(R), "Sequence Length: " + str(sequenceLength), "Sequence: " + str(donorSequence)

# Modified functions for input
def standardCassette(PromoterName,TerminatorName, orfName, orfSeq):
	
	#first, the promoter
	print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
	print("First, which PROMOTER do you want to use, e.g., TDH3")
	
	PromoterGeneRec=fetchGene(PromoterName)
	PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
	PromoterRec.id=PromoterRec.id+"ps"
	
	
	#second, the terminator
	print("Which TERMINATOR do you want to use, e.g., ADH1")
	TerminatorGeneRec=fetchGene(TerminatorName)
	TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
	TerminatorRec.id=TerminatorRec.id+"ts"
	
	
	#and last, the gene
	print("What is the name of your gene, e.g., KlGapDH")
	
	print("What's the sequence")
	
	orfRecord=SeqRecord(Seq(orfSeq, SingleLetterAlphabet()), id=orfName)
	
	insertRec=[PromoterRec,orfRecord,TerminatorRec]
	return PromoterRec, orfRecord, TerminatorRec

def buildCassette(PromoterName, TerminatorName, orfName, orfSeq):
	
	#first, the promoter
	print("I'm going to build a standard cassette in which promoter is 600nt, terminator 250nt.") 
	print("First, which PROMOTER do you want to use, e.g., TDH3")
	
	PromoterGeneRec=fetchGene(PromoterName)
	PromoterRec=fetchNeighbor(PromoterGeneRec,"upstream",600)
	PromoterRec.id=PromoterRec.id+"ps"
	
	
	#second, the terminator
	print("Which TERMINATOR do you want to use, e.g., ADH1")
	TerminatorGeneRec=fetchGene(TerminatorName)
	TerminatorRec=fetchNeighbor(TerminatorGeneRec,"downstream",250)
	TerminatorRec.id=TerminatorRec.id+"ts"
	
	#and last, the gene
	print("What is the name of your gene, e.g., KlGapDH")
	
	print("What's the sequence")
	
	orfRecord=SeqRecord(Seq(orfSeq), id=orfName)
	
	insertRec=[PromoterRec,orfRecord,TerminatorRec]
	return PromoterRec, orfRecord, TerminatorR

def variableCassette(geneList, seqList, toVary="", variants=[], variantSeq=[]):
	# Store both name and sequence in a SeqRecord
	# Append them to a list
	# Return list as fragments to be stitched
	if toVary != "":
		toVary = int(toVary)
	records = []

	counter = 0
	for gene in geneList:
		name = gene;
		sequence = seqList[counter]
		Rec = SeqRecord(Seq(sequence, SingleLetterAlphabet()), id = str(counter+1))
		Rec.name = name
		records.append(Rec)
		counter += 1

	variantRecords = []
	variantRecords.append(records)
	
	# Executes if variants is not empty
	counter = 0
	if variants != []:
		for variant in variants:
			name = variant
			sequence = variantSeq[counter]
			Rec = SeqRecord(Seq(sequence, SingleLetterAlphabet()), id = str(counter+1))
			Rec.name = name
			# Make a copy of the original, switch the fragments and add it to the list. 
			# Deep-copy ensures there are no pointer issues
			tempVariant = copy.deepcopy(records)
			tempVariant[toVary - 1] = Rec
			variantRecords.append(copy.deepcopy(tempVariant))
			counter += 1

	# Returns a list of lists of the answers.
	answer = [[stitch(variantRecords[0])]]
	variants = []
	for n in range(len(variantRecords)-1):
		frags = variantRecords[n + 1][toVary - 2: toVary]
		variantStitch = [stitch(frags)]
		answer.append(variantStitch)

	return answer
