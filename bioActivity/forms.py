from django import forms
from django.core.exceptions import ValidationError
import os

class InsertFile(forms.Form):
    file = forms.FileField(label='File')

    def clean_file(self):
        file = self.cleaned_data['file']
        
        ext = os.path.splitext(file.name)[1]
        valid_extensions = ['.csv', '.xlsx', '.txt']
        if not ext.lower() in valid_extensions:
            raise ValidationError('Unsupported file extension. Only CSV, XLSX, TXT files are allowed.')

        return file