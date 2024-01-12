from django.db import models

# Create your models here.

class SmilesData(models.Model):
    smiles = models.CharField(max_length=256, unique=True)
    pic50 = models.FloatField()
    bio_class = models.BooleanField(default=False)

class CurrSmilesData(models.Model):
    smiles = models.CharField(max_length=256, unique=True)
    pic50 = models.FloatField()
    bio_class = models.BooleanField(default=False)