# Generated by Django 5.0.1 on 2024-01-08 08:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioActivity', '0006_smilesdata_bio_class_alter_smilesdata_pic50'),
    ]

    operations = [
        migrations.AlterField(
            model_name='smilesdata',
            name='pic50',
            field=models.FloatField(),
        ),
    ]
