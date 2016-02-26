# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('cassette', '0002_auto_20150717_1112'),
    ]

    operations = [
        migrations.CreateModel(
            name='Results',
            fields=[
                ('id', models.AutoField(serialize=False, auto_created=True, primary_key=True, verbose_name='ID')),
                ('sequenceLength', models.CharField(max_length=1000000)),
                ('donorSequence', models.CharField(max_length=1000000)),
                ('Lup', models.CharField(max_length=1000000)),
                ('Rup', models.CharField(max_length=1000000)),
                ('Ldown', models.CharField(max_length=1000000)),
                ('Rdown', models.CharField(max_length=1000000)),
                ('L', models.CharField(max_length=1000000)),
                ('R', models.CharField(max_length=1000000)),
            ],
        ),
        migrations.RemoveField(
            model_name='answer',
            name='question',
        ),
        migrations.RemoveField(
            model_name='choice',
            name='question',
        ),
        migrations.DeleteModel(
            name='Answer',
        ),
        migrations.DeleteModel(
            name='Choice',
        ),
        migrations.DeleteModel(
            name='Question',
        ),
    ]
