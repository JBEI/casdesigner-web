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
import numpy as np
from textwrap import fill
from pygenome import sg


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
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq, rendered = editEmpty(name, sequence, cutsite)


			elif choice1 == "2":
				orfName = request.POST['orfName1']
				orfSeq = request.POST['orfSeq1']
				promoterName = request.POST['promoterName1']
				terminatorName = request.POST['terminatorName1']
				NtagName = request.POST['NtagName1']
				CtagName = request.POST['CtagName1']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq, rendered = editEmpty(orfName, orfSeq, cutsite, promoterName, terminatorName, NtagName, CtagName)

			request.session['Lup'] = Lup
			request.session['Rup'] = Rup
			request.session['Ldown'] = Ldown
			request.session['Rdown'] = Rdown
			request.session['L'] = L
			request.session['R'] = R
			request.session['seqLen'] = seqLen
			request.session['donorSeq'] = rendered
			request.session['rendered'] = rendered
			return HttpResponseRedirect('/cassette/results')

		elif choice0 == "2":
			locus = request.POST['locus']
			choice2 = request.POST['choices2']

			if choice2 == "1":
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq, rendered = editExisting(locus, 1)
			elif choice2 == "2":
				name2 = request.POST['name2']
				sequence2 = request.POST['sequence2']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq, rendered = editExisting(locus, 2, NewGeneName=name2, NewGeneSeq=sequence2)
			elif choice2 == "3":
				promoterName = request.POST['promoterName2']
				terminatorName = request.POST['terminatorName2']
				orfName = request.POST['orfName2']
				orfSeq = request.POST['orfSeq2']
				NtagName = request.POST['NtagName2']
				CtagName = request.POST['CtagName2']
				Lup, Rup, Ldown, Rdown, L, R, seqLen, donorSeq, rendered = editExisting(locus, 3, promoter=promoterName,terminator=terminatorName,NewGeneName=orfName,NewGeneSeq=orfSeq, Ntag=NtagName, Ctag=CtagName)

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
			request.session['donorSeq'] = rendered
			request.session['rendered'] = rendered
			return HttpResponseRedirect('/cassette/results')

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

			return HttpResponseRedirect('/cassette/customResults')

	# If a GET (or any other method)we'll create a blank form)
	else:
		form = TrueIndexForm()
	return render(request, 'cassette/index.html', {'form': form})

def results(request):
	Lup = request.session.get('Lup')
	Rup = request.session.get('Rup')
	Ldown = request.session.get('Ldown')
	Rdown = request.session.get('Rdown')
	L = request.session.get('L')
	R = request.session.get('R')
	seqLen = request.session.get('seqLen')
	donorSeq = request.session.get('donorSeq')
	rendered = request.session.get('rendered')

	return render(request, 'cassette/results.html', {'Lup': Lup, 'Rup': Rup, 'Ldown': Ldown, 'Rdown': Rdown, 'L': L, 'R': R, 'seqLen': seqLen, 'donorSeq': donorSeq, 'rendered': rendered })

def customResults(request):
	answer = request.session.get('customAnswer')
	original = answer[0][0]
	variants = answer[1:]
	return render(request, 'cassette/customResults.html', {'original' : original, 'variants' : variants})


#---------------  Helper Functions  ---------------------#
def processTags (nTerminalTags, cTerminalTags):
	nTerminalTags = nTerminalTags.split(',')
	nTerminalTags = list(filter(None, map(lambda x: x.strip(), nTerminalTags)))
	cTerminalTags = cTerminalTags.split(',')
	cTerminalTags = list(filter(None, map(lambda x: x.strip(), cTerminalTags)))
	return nTerminalTags, cTerminalTags

def makeTags(text):
	def validTags():
		tagFrame = pd.read_excel("/usr/local/casdesigner/ProteinTags.xlsx")
		return tagFrame["tagName"].tolist()

	def tagFinder(tagName):
		tagFrame = pd.read_excel("/usr/local/casdesigner/ProteinTags.xlsx", index_col="tagName")

		sequence = tagFrame.loc[tagName, 'sequence']

		return sequence

	if text in validTags():
		return SeqRecord(Seq(tagFinder(text)), id=text + "tag")
	else:
		n, s = text.split()
		return SeqRecord(Seq(s), id=n + "tag")


def moveCodons (orfSeq, orfName, nTerminalTags, cTerminalTags):
	def removeStartCodon(CDS):
		# this alternate method would be useful for error handling
		#         i = CDS.find('ATG')
		#         return (CDS[:i], CDS[:i])
		return CDS[:3], CDS[3:]

	def removeStopCodon(CDS):
		#this alternate method would be useful for error handling
		#         stops = ['TAG', 'TAA', 'TAG']
		#         for codon in stops:
		#             i = CDS.rfind(codon)
		#             print(i)
		#             if i != -1:
		#                 return (CDS[:i], CDS[:i])
		l = len(CDS)
		return CDS[:l - 3], CDS[l - 3:]
	nTerminals = len(nTerminalTags)
	cTerminals = len(cTerminalTags)

	if nTerminals > 0:
		nTerminalTags = list(map(makeTags, nTerminalTags))
		start, orfSeq = removeStartCodon(orfSeq)
		nTerminalTags[0] = start + nTerminalTags[0]
	if cTerminals > 0:
		cTerminalTags = list(map(makeTags, cTerminalTags))
		orfSeq, stop = removeStopCodon(orfSeq)
		cTerminalTags[cTerminals - 1] += stop
	areNTags = [i + 1 for i in range(nTerminals)]
	areCTags = [i + 2 + nTerminals for i in range(cTerminals)]

	orfRecord = SeqRecord(Seq(orfSeq), id=orfName)

	return nTerminalTags + [orfRecord] + cTerminalTags, areNTags, areCTags

def editEmpty(name: object, sequence: object, cutname: object, promoter: object = None, terminator: object = None, Ntag = "", Ctag = "") -> object:
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

	ChromosomeSeq=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", location), "fasta").seq

	if ChromosomeSeq.find(cutSequence)==-1:
		ChromosomeSeq=ChromosomeSeq.reverse_complement()

	StartIndex=ChromosomeSeq.find(cutSequence)
	EndIndex=StartIndex+34

	UpSeq=ChromosomeSeq[StartIndex-HomologyLength:StartIndex]
	DownSeq=ChromosomeSeq[EndIndex:EndIndex+HomologyLength]

	UpHomRec = SeqRecord(UpSeq, id=cutname)
	DownHomRec = SeqRecord(DownSeq, id=cutname)
	Ntag, Ctag = processTags(Ntag, Ctag)
	if len(Ntag)>0 or len(Ctag)>0:
		orfRecords, areNTags, areCTags  = moveCodons(sequence, name, Ntag, Ctag)
	else:
		orfRecords = [SeqRecord(Seq(sequence, SingleLetterAlphabet()), id=name)]
		areNTags = []
		areCTags = []

	if promoter is None:
		fragments = [UpHomRec] + orfRecords + [DownHomRec]
	else:
		PromoterGeneRec = fetchGene(promoter)
		PromoterRec = fetchNeighbor(PromoterGeneRec,"upstream",600)
		PromoterRec.id = PromoterRec.id + "ps"

		TerminatorGeneRec = fetchGene(terminator)
		TerminatorRec = fetchNeighbor(TerminatorGeneRec,"downstream",250)
		TerminatorRec.id = TerminatorRec.id + "ts"

		fragments = [UpHomRec, PromoterRec] + orfRecords + [TerminatorRec, DownHomRec]
		areNTags = list(map(lambda x: x + 1, areNTags))
		areCTags = list(map(lambda x: x + 1, areCTags))
		# JPNTODO the build cassette function should be used for when building a cassette...
	return stitch(fragments, areNTags, areCTags)

def editExisting(name, option, promoter = None, terminator = None, NewGeneName = "", NewGeneSeq = "", Ntag="", Ctag=""):
	OrigGeneRecord = fetchGene(name)
	UpHomRec = fetchNeighbor(OrigGeneRecord, "upstream", HomologyLength)
	DownHomRec = fetchNeighbor(OrigGeneRecord, "downstream", HomologyLength)
	cleanDeletion = False
	areNTags = []
	areCTags = []
	if option == 1:
		fragments = [UpHomRec, DownHomRec]
		cleanDeletion = True
	elif option == 2:
		InsertRec = SeqRecord(Seq(NewGeneSeq, SingleLetterAlphabet()), id=NewGeneName)
		fragments = [UpHomRec, InsertRec, DownHomRec]
	elif option == 3:
		PromoterRec, orfRecord, TerminatorRec = standardCassette(promoter, terminator, NewGeneName, NewGeneSeq)
		Ntag, Ctag = processTags(Ntag, Ctag)
		if len(Ntag) > 0 or len(Ctag) > 0:
			orfRecords, areNTags, areCTags = moveCodons(NewGeneSeq, NewGeneName, Ntag, Ctag)
		else:
			orfRecords = [orfRecord]
		fragments = [UpHomRec, PromoterRec] + orfRecords + [TerminatorRec, DownHomRec]
		areNTags = list(map(lambda x: x + 1, areNTags))
		areCTags = list(map(lambda x: x + 1, areCTags))
		#JPNTODO the build cassete function should be actually used... and take tags into account
	elif option == 4:
		pass
	elif option == 5:
		pass

	return stitch(fragments, areNTags, areCTags, deletion=cleanDeletion)

def fetchGene(GeneName):

#let's create a record for the oldGene
    DesiredSeq = sg.gene[GeneName].cds.seq
    
    GeneRecord = SeqRecord(DesiredSeq, name=GeneName)
    
    #now let's add some more information to make it useful
    GeneRecord.features=str(sg.gene[GeneName])[-7:]
    GeneRecord.description=sg.gene[GeneName].short_description
    
    return GeneRecord

def fetchNeighbor(NeighborRecord, direction, distance):


	# let's load the appropriate chromosome file. The record of the gene we looked up
	# contains in the "features" the systematic name, wherein the second letter
	# corresponds to chromosome number, e.g., 1=A etc
	if NeighborRecord.features[1]=="A":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer01.fasta"), "fasta")
	if NeighborRecord.features[1]=="B":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer02.fasta"), "fasta")
	if NeighborRecord.features[1]=="C":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer03.fasta"), "fasta")
	if NeighborRecord.features[1]=="D":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer04.fasta"), "fasta")
	if NeighborRecord.features[1]=="E":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer05.fasta"), "fasta")
	if NeighborRecord.features[1]=="F":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer06.fasta"), "fasta")
	if NeighborRecord.features[1]=="G":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer07.fasta"), "fasta")
	if NeighborRecord.features[1]=="H":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer08.fasta"), "fasta")
	if NeighborRecord.features[1]=="I":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer09.fasta"), "fasta")
	if NeighborRecord.features[1]=="J":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer10.fasta"), "fasta")
	if NeighborRecord.features[1]=="K":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer11.fasta"), "fasta")
	if NeighborRecord.features[1]=="L":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer12.fasta"), "fasta")
	if NeighborRecord.features[1]=="M":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer13.fasta"), "fasta")
	if NeighborRecord.features[1]=="N":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer14.fasta"), "fasta")
	if NeighborRecord.features[1]=="O":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer15.fasta"), "fasta")
	if NeighborRecord.features[1]=="P":
		ChromosomeRec=SeqIO.read(os.path.join(PROJECT_ROOT, "chromosomes", "Scer16.fasta"), "fasta")



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

	while mp <= PrimerMaxTm and length < PrimerMaxLen:
		primer = primer + seq[length]
		mp = MeltingTemp.Tm_staluc(primer)
		length += 1

	return primer

def overhangPrimer(currRecord,prevSeq, linker = ""):
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
	ohprimer=prevSeq.seq[-maxOhLen:] + linker + primer #we add the .seq so that it returns a string
	if len(ohprimer):
		ohprimer = ohprimer[-60:]
	return ohprimer

def stitch(fragments, areNTags = [], areCTags = [], deletion=False):
	nLinker = Seq("GGAGGTGGTGGAGGTGGA")
	cLinker = Seq("GGTAGCGGTAGCGGCAGC")
	#this function takes seq records and prints primers

	#let's make an empty sequence file
	Nfrags=len(fragments)
	donor=Seq("")
	index=[]
	print("")
	for i in range(0, Nfrags):
		if i in areCTags:
			# add the linker before bin
			#JPNTODO check that this annotates the linkers correctly in GB output
			donor += SeqRecord(cLinker, id="cLinker")
		donor = donor + fragments[i]
		if i in areNTags:
			# add linker after bin
			# JPNTODO check that this annotates the linkers correctly in GB output
			donor += SeqRecord(nLinker, id="nLinker")
	# Dummy assignment setup to allow for compilation
	Lup = ""
	Rup = ""
	Ldown = ""
	Rdown = ""
	L = ""
	R = ""

	# JPNTODO oh my god this needs some work, I put my stuff into both loops since the first is returned (to where???)
	# and the second one seems to make the page output

	# I'm not sure what this loop does (don't confuse it with the other in this method. If you can
	# factor it out, be my guest.
	for i in range (0, Nfrags):
		FLinker = ""
		RLinker = ""
		# this section introduces the protein linkers for tags
		if i in areNTags:
			# reverse primer will need a liniker
			RLinker = nLinker.reverse_complement()
		if i - 1 in areNTags:
			# forward primer will need a linker
			FLinker = nLinker
		if i + 1 in areCTags:
			# reverse primer will need a linker
			RLinker = cLinker.reverse_complement()
		if i in areCTags:
			# forward primer will need a linker
			FLinker = cLinker
		if i==0:
			Lup = "Lup"+ fragments[i].id + " " + getPrimer(donor)
			Rup = "Rup"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement(), RLinker)
		elif i==Nfrags-1:
			Ldown = "Ldown"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1], FLinker)
			Rdown = "Rdown"+ fragments[i].id + " " + getPrimer(donor.reverse_complement())
		else:
			L = "L"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1], FLinker)
			R = "R"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement(), RLinker)

	sequenceLength = len(donor.seq)
	donorSequence = donor.seq

	# Separate rendering stage for custom cassettes
	rendered = "<pre>"

	rendered = rendered +"Here are the primers to amplify your fragments and construct your donor DNA cassette:\n\n"
	# The names include information on the homology provided by the overhang
	# Note that some primers don't have overhangs
	for i in range (0, Nfrags):
		FLinker = ""
		RLinker = ""
		# this section introduces the protein linkers for tags
		if i in areNTags:
			# reverse primer will need a liniker
			RLinker = nLinker.reverse_complement()
		if i - 1 in areNTags:
			# forward primer will need a linker
			FLinker = nLinker
		if i + 1 in areCTags:
			# reverse primer will need a linker
			RLinker = cLinker.reverse_complement()
		if i in areCTags:
			# forward primer will need a linker
			FLinker = cLinker

		if i==0:
			if deletion:
					rightFragmentId = "del"
			else:
				rightFragmentId = fragments[i+1].id

			rendered = rendered +"F-up"+ fragments[i].id + " " + getPrimer(donor) + "\n"
			rendered = rendered +"R-up"+ fragments[i].id + "(" + rightFragmentId + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement(), RLinker) + "\n"
		elif i==Nfrags-1:
			if deletion:
					leftFragmentId = "del"
			else:
				leftFragmentId = fragments[i-1].id

			rendered = rendered +"F-dn"+ fragments[i].id + "(" + leftFragmentId + ") " + overhangPrimer(fragments[i],fragments[i-1], FLinker) + "\n"
			rendered = rendered +"R-dn"+ fragments[i].id + " " + getPrimer(donor.reverse_complement()) + "\n"
		else:
			rendered = rendered +"F-"+ fragments[i].id + "(" + fragments[i-1].id + ") " + overhangPrimer(fragments[i],fragments[i-1], FLinker) + "\n"
			rendered = rendered +"R-"+ fragments[i].id + "(" + fragments[i+1].id + ") " + overhangPrimer(fragments[i].reverse_complement(),fragments[i+1].reverse_complement(), RLinker) + "\n"

	rendered = rendered +"\n\nThe size and sequence of your donor DNA is below.\n\n"

	rendered = rendered + "> " + str(len(donor.seq)) + "\n"

	rendered = rendered + fill(str(donor.seq), 80)

	rendered = rendered + "</pre>"

	return str(Lup), str(Rup), str(Ldown), str(Rdown), str(L), str(R), "Sequence Length: " + str(sequenceLength), "Sequence: " + str(donorSequence), str(rendered)


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
