from django.db import models

# Create your models here.

class SmilesData(models.Model):
    smiles = models.CharField(max_length=256)
    pic50 = models.IntegerField()