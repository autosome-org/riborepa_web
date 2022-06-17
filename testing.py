from flask import Flask, render_template, url_for, flash, request, redirect
from flask_bootstrap import Bootstrap
import pandas as pd

classes_table = pd.read_table('./static/tables/classes_table.tsv', sep='\t')
classes_dict = classes_table.set_index('class').to_dict('index')
web_classes_dict = dict()
for key in classes_dict:
    web_class = classes_dict[key]['web_class']
    if web_class not in web_classes_dict.keys():
        web_classes_dict[web_class] = []
    web_classes_dict[web_class].append(key)

