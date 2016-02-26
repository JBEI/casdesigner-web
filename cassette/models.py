from django.db import models

"""
This is where we create the models. Note that the first question will
be a choice based question.

If the question is a question that acquires a 
"""

# Create your models here.
class Results(models.Model):
	sequenceLength = models.CharField(max_length=1000000)
	donorSequence = models.CharField(max_length=1000000)
	Lup = models.CharField(max_length=1000000)
	Rup = models.CharField(max_length=1000000)
	Ldown = models.CharField(max_length=1000000)
	Rdown = models.CharField(max_length=1000000)
	L = models.CharField(max_length=1000000)
	R = models.CharField(max_length=1000000)