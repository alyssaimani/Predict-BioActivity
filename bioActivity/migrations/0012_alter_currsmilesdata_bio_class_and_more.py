# Generated by Django 5.0.1 on 2024-01-11 16:16

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioActivity', '0011_alter_currsmilesdata_bio_class_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='currsmilesdata',
            name='bio_class',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='smilesdata',
            name='bio_class',
            field=models.BooleanField(default=False),
        ),
    ]
