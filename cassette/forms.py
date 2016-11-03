from django import forms
from django.utils.translation import ugettext as _
from django.utils.safestring import mark_safe
from django.forms.widgets import RadioSelect

class IndexForm(forms.Form):
	CHOICES = [['1','Insert DNA into a characterized locus'],
		       ['2','Edit an existing gene, e.g., delete or replace'],
		  	   ['3','Build a custom cassette']]

	choices = forms.ChoiceField(required=True, 
								widget = forms.RadioSelect, 
								choices = CHOICES,
								label = False
								)

class CharLocusForm(forms.Form):
	cutsite = forms.CharField(max_length = 10)

	CHOICES = [['1', 'Insert a pre-built custom cassette'],
			   ['2', 'Construct a cassette using standard promoters']]

	choices = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES, 
								label = 'What would you like to do to this cutsite?')

class CharLocusPrebuiltForm(forms.Form):
	name = forms.CharField(max_length = 100,
						   label = 'What is the name of the ORF?:')
	sequence = forms.CharField(max_length = 1000000,
							   label = 'What is the sequence of the ORF?:')

class CharLocusStandardForm(forms.Form):
	promoterName = forms.CharField(max_length = 100,
								   label = 'Which PROMOTER would you like to use (e.g. TDH3)?: ')
	terminatorName = forms.CharField(max_length = 100,
									 label = 'Which TERMINATOR would you like to use (e.g. ADH1)?: ')
	orfName = forms.CharField(max_length = 100, 
						   label = 'What is the name of your custom gene (e.g. K1GapDH)?: ')
	orfSeq = forms.CharField(max_length = 1000000, 
							 label = 'What is the sequence of the ORF?: ')

class ExistingLocusForm(forms.Form):
	locus = forms.CharField(max_length = 100,
							label = 'Which locus do you want to edit (e.g. OAF1): ')
	CHOICES = [['1', 'Delete the protein coding DNA sequence (CDS)'],
			   ['2', 'Replace the CDS with another CDS, or pre-made fragment'],
			   ['3', 'Replace the CDS with a standard cassette I will build for you'],
			   ['4', 'Replace the CDS with a custom cassette (NOT YET IMPLEMENTED)'],
			   ['5', 'Replace a specified region near your target gene (NOT YET IMPLEMENTED)']]

	choices = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES,
								label = 'What do you want to do to this gene?: '
								)

	insertGene = forms.CharField(max_length = 100, 
							     label = 'Name of gene to insert: ')
	insertSeq = forms.CharField(max_length = 100,
								label = 'Sequence of gene to insert: ')

	promoter = forms.CharField(max_length = 100,
							   label = 'Name of PROMOTER (e.g. TDH3): ')
	terminator = forms.CharField(max_length = 100,
								 label = 'Name of TERMINATOR (e.g. ADH1): ')
	Ntag = forms.CharField(max_length = 1000000,
							label = 'What N-terminal tag do you want to use: ')
	Ctag = forms.CharField(max_length = 1000000,
							label = 'What C-terminal tag do you want to use: ')

class TrueIndexForm(forms.Form):
	# Initial Choice
	CHOICES0 = [['1','Insert DNA into a characterized locus'],
		       ['2','Edit an existing gene, e.g., delete or replace'],
		  	   ['3','Build a custom cassette']]

	choices0 = forms.ChoiceField(required=True, 
								widget = forms.RadioSelect, 
								choices = CHOICES0,
								label = False
								)

	# Choosing 0.1
	cutsite = forms.CharField(max_length = 10)

	CHOICES1 = [['1', 'Insert a pre-built custom cassette'],
			   ['2', 'Construct a cassette using standard promoters']]

	choices1 = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES1, 
								label = 'What would you like to do to this cutsite?')

	# Choosing 1.1
	name1 = forms.CharField(max_length = 100,
						   label = 'What is the name of the ORF?:')
	sequence1 = forms.CharField(max_length = 1000000,
							   label = 'What is the sequence of the ORF?:')

	# Choosing 1.2
	promoterName1 = forms.CharField(max_length = 100,
										label = 'Which PROMOTER would you like to use (e.g. TDH3)?: ')
	terminatorName1 = forms.CharField(max_length = 100,
										label = 'Which TERMINATOR would you like to use (e.g. ADH1)?: ')
	orfName1 = forms.CharField(max_length = 100, 
										label = 'What is the name of your custom gene (e.g. K1GapDH)?: ')
	orfSeq1 = forms.CharField(max_length = 1000000, 
										label = 'What is the sequence of the ORF?: ')
	NtagName1 = forms.CharField(max_length = 1000000,
										label = 'What N-terminal tag do you want to use: ')
	CtagName1 = forms.CharField(max_length = 1000000,
										label = 'What C-terminal tag do you want to use: ')

	# Choosing 0.2
	locus = forms.CharField(max_length = 100,
							label = 'Which locus do you want to edit (e.g. OAF1): ')

	CHOICES2 = [['1', 'Delete the protein coding DNA sequence (CDS)'],
			   ['2', 'Replace the CDS with another CDS, or pre-made fragment'],
			   ['3', 'Replace the CDS with a standard cassette I will build for you'],
			   ['4', 'Replace the CDS with a custom cassette (NOT YET IMPLEMENTED)'],
			   ['5', 'Replace a specified region near your target gene (NOT YET IMPLEMENTED)']]

	choices2 = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES2,
								label = 'What do you want to do to this gene?: '
								)

	# Choosing 2.1 => Goes to results screen.

	# Choosing 2.2
	name2 = forms.CharField(max_length = 100,
						   label = 'What is the name of the ORF?:')
	sequence2 = forms.CharField(max_length = 1000000,
							   label = 'What is the sequence of the ORF?:')

	# Choosing 2.3
	promoterName2 = forms.CharField(max_length = 100,
										label = 'Which PROMOTER would you like to use (e.g. TDH3)?: ')
	terminatorName2 = forms.CharField(max_length = 100,
										label = 'Which TERMINATOR would you like to use (e.g. ADH1)?: ')
	orfName2 = forms.CharField(max_length = 100, 
										label = 'What is the name of your custom gene (e.g. K1GapDH)?: ')
	orfSeq2 = forms.CharField(max_length = 1000000, 
										label = 'What is the sequence of the ORF?: ')
	NtagName2 = forms.CharField(max_length = 1000000,
										label = 'What N-terminal tag do you want to use: ')
	CtagName2 = forms.CharField(max_length = 1000000,
										label = 'What C-terminal tag do you want to use: ')

	# Choosing 2.4 (NOT IMPLEMENTED)

	# Choosing 2.5 (NOT IMPLEMENTED)

	# Choosing 0.3
	CHOICES3 = [['1', 'Insert a pre-built custom cassette'],
			   ['2', 'Construct a cassette using standard promoters']]
	choices3 = forms.ChoiceField(required = True,
								widget = forms.RadioSelect,
								choices = CHOICES3,
								label = 'How do you want your cassettes?: '
								)

	# Choosing 3.1
	N1 = forms.CharField(max_length = 100,
				  label = 'How many pieces are in your custom cassette?: ')
	# Choosing 3.2
	N2 = forms.CharField(max_length = 100,
				   label = 'How many pieces are in your custom cassette?: ')

	toVary = forms.CharField(max_length = 100,
					   label = 'Which piece do you want to vary?: ')

	variants = forms.CharField(max_length = 100, 
						 label = 'How many variants of that piece would you like?: ')
