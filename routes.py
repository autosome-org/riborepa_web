from flask import Flask, render_template, url_for, flash, request, redirect
from flask_bootstrap import Bootstrap
import pandas as pd
from wtforms import Form, StringField, SelectField
from Bio.SeqUtils.ProtParam import ProteinAnalysis

app = Flask(__name__)  # , template_folder='/home/artyom/Desktop/web/templates')


@app.route('/', methods=['GET', 'POST'])
def home():
    q = request.args.get('q')
    classes_table = pd.read_table('./static/tables/classes_table.tsv', sep='\t')
    classes_dict = classes_table.set_index('class').to_dict('index')
    filtered_classes, filtered_families, filtered_elements = [], [], []
    for key in classes_dict:
        filtered_classes.append(key)
        filtered_families += classes_dict[key]['families'].split(',')
        filtered_elements += classes_dict[key]['elements'].split(',')
    found_classes, found_families, found_elements = [], [], []
    if q:
        found_classes = [i for i in filtered_classes if q.lower() in i.lower()]
        found_families = [i for i in filtered_families if q.lower() in i.lower()]
        found_elements = [i for i in filtered_elements if q.lower() in i.lower()]
    else:
        q = 'NA'
    if found_classes + found_families + found_elements == [] and q != 'NA':
        q = 'empty'
    return render_template('index.html', q=q, found_classes=found_classes, found_families=found_families,
                           found_elements=found_elements)


@app.route('/class/<class_s>')
def class_rb(class_s):
    classes_table = pd.read_table('./static/tables/classes_table.tsv', sep='\t')
    classes_dict = classes_table.set_index('class').to_dict('index')
    web_classes_dict = dict()
    for key in classes_dict:
        web_class = classes_dict[key]['web_class']
        classes_dict[key]['families'] = classes_dict[key]['families'].split(',')
        classes_dict[key]['elements'] = classes_dict[key]['elements'].split(',')
        if web_class not in web_classes_dict.keys():
            web_classes_dict[web_class] = []
        web_classes_dict[web_class].append(key)
    selected_classes = web_classes_dict[class_s]
    return render_template('class.html', class_s=class_s, selected_classes=selected_classes, classes_dict=classes_dict)


@app.route('/family/<family>')
def family(family):
    families_table = pd.read_table('./static/tables/families_table.tsv', sep='\t')
    families_dict = families_table.set_index('family').to_dict('index')
    elements = families_dict[family]['elements'].split(',')
    paths = []
    for element in elements:
        paths.append(element + '.png')
    return render_template('family.html', family=family, elements=elements)

@app.route('/element/<element>')
def element(element):
    merged_statistics_table = pd.read_table(
        './static/tables/Merged_statistics_filtered_with_classes_for_visualization.tsv', sep='\t')
    merged_statistics_dict = merged_statistics_table.set_index('ORF_id').to_dict('index')
    element_tables = pd.read_table('./static/tables/element_table.tsv', sep='\t')
    element_dict = element_tables.set_index('element').to_dict('index')
    element_frames = element_dict[element]['frame'].split(',')
    for frame in element_frames:
        extended_form = merged_statistics_dict[frame]['the_longest_ORF']
        if extended_form not in element_frames:
            element_frames.append(extended_form)
    element_class = element_dict[element]['el_class']
    element_family = element_dict[element]['el_family']
    frames_dict = dict()
    for frame in element_frames:
        frames_dict[frame] = dict()
        try:
            table = pd.read_csv("./static/tables/repeat_masker_blast_results/" + frame + '.tsv', sep='\t')
            table = table.to_html(header="true", index_names='false', justify='center', index=False)
            repmask_table = table[table.find('\n') + 1:table.rfind('\n')]
        except:
            repmask_table = 'NA'
        try:
            table = pd.read_csv("./static/tables/nblast_results/" + frame + '.tsv', sep='\t')
            table = table.to_html(header="true", index_names='false', justify='center', index=False)
            nblast_table = table[table.find('\n') + 1:table.rfind('\n')]
        except:
            nblast_table = 'NA'
        try:
            table = pd.read_csv("./static/tables/blastp_results/" + frame + '.tsv', sep='\t')
            table = table.to_html(header="true", index_names='false', justify='center', index=False)
            pblast_table = table[table.find('\n') + 1:table.rfind('\n')]
        except:
            pblast_table = 'NA'
        datasets = merged_statistics_dict[frame]['Datasets_sorted_by_priority'].split(',')
        frames_dict[frame]['repmask_table'] = repmask_table
        frames_dict[frame]['nblast_table'] = nblast_table
        frames_dict[frame]['pblast_table'] = pblast_table
        frames_dict[frame]['datasets'] = datasets
    return render_template('element.html', element=element, element_family=element_family, element_class=element_class,
                           frames_dict=frames_dict)

@app.route('/ORF/<orf>')
def orf(orf):
    merged_statistics_table = pd.read_table(
        './static/tables/Merged_statistics_filtered_with_classes_for_visualization.tsv', sep='\t')
    merged_statistics_dict = merged_statistics_table.set_index('ORF_id').to_dict('index')[orf]
    frame_dict = dict()
    try:
        table = pd.read_csv("./static/tables/repeat_masker_blast_results/" + orf + '.tsv', sep='\t')
        table = table.to_html(header="true", index_names='false', justify='center', index=False)
        repmask_table = table[table.find('\n') + 1:table.rfind('\n')]
    except:
        repmask_table = 'NA'
    try:
        table = pd.read_csv("./static/tables/nblast_results/" + orf + '.tsv', sep='\t')
        table = table.to_html(header="true", index_names='false', justify='center', index=False)
        nblast_table = table[table.find('\n') + 1:table.rfind('\n')]
    except:
        nblast_table = 'NA'
    try:
        table = pd.read_csv("./static/tables/blastp_results/" + orf + '.tsv', sep='\t')
        table = table.to_html(header="true", index_names='false', justify='center', index=False)
        pblast_table = table[table.find('\n') + 1:table.rfind('\n')]
    except:
        pblast_table = 'NA'

    fasta_nucleotide = '>' + orf + '_nt | length:' + str(
        merged_statistics_dict['length']) + ' | Start codon sequence [-3:3], context:  ' + merged_statistics_dict[
                           'start_codon_nucl'] + ', ' + merged_statistics_dict['start_codon_context'] + '\n' + \
                       merged_statistics_dict['nucleotide_seq']
    frame_dict['fasta_nucleotide'] = fasta_nucleotide

    aa_seq = merged_statistics_dict['aa_seq']
    fasta_aa = '>' + orf + '_aa | length:' + str(int(merged_statistics_dict['length'] / 3)) + ' | Molecular weight (Da): ' + \
               str("%0.2f" % ProteinAnalysis(aa_seq.replace('X', '').replace('J','')).molecular_weight()) + '\n' + \
               aa_seq
    frame_dict['fasta_aa'] = fasta_aa

    filters_stat = pd.read_csv("./static/tables/filtering_stat/" + orf + '.tsv', sep='\t')
    filters_stat = filters_stat.to_html(header="true", index_names='false',  index=False)
    filters_stat = filters_stat[filters_stat.find('\n') + 1:filters_stat.rfind('\n')]

    datasets = merged_statistics_dict['Datasets_sorted_by_priority'].split(',')
    frame_dict['repmask_table'] = repmask_table
    frame_dict['nblast_table'] = nblast_table
    frame_dict['pblast_table'] = pblast_table
    frame_dict['filtering_stat_table'] = filters_stat
    frame_dict['datasets'] = datasets
    element = merged_statistics_dict['element_name']
    element_family = merged_statistics_dict['element_family']
    element_class = merged_statistics_dict['element_class']
    return render_template('orf.html', orf=orf, element=element, element_family=element_family,
                           element_class=element_class,
                           frame_dict=frame_dict)
    return render_template('orf.html')

@app.route('/ORFs')
def orfs():
    labels_dict = {'1': 'Only FLOSS or FOOPS filter is passed at least in one dataset',
                   '2': 'Only MS filter is passed',
                   '3': 'Both FLOSS and FOOPS filters are passed at least in one dataset',
                   '4': 'MS filter & Both FLOSS and FOOPS filters are passed at least in one dataset'}
    classes_tables_dict = dict()
    for orf_class in labels_dict.keys():
        table = pd.read_csv("./static/tables/ORFs_tables/class_" + orf_class + ".tsv", sep='\t')
        table = table.to_html(header="true", index_names='false', justify='center', index=False)
        table = table[table.find('\n') + 1:table.rfind('\n')]
        classes_tables_dict[orf_class] = table

    return render_template('orfs.html', labels_dict = labels_dict, classes_tables_dict = classes_tables_dict)

@app.route('/MS')
def ms():
    MS_table = pd.read_csv("./static/tables/MS_filtered_ORFs.tsv", sep='\t')
    MS_table = MS_table.to_html(header="true", index_names='false', justify='center', index=False)
    MS_table = MS_table[MS_table.find('\n') + 1:MS_table.rfind('\n')]
    MS_table = MS_table.replace('<tr style="text-align: right;">', '<tr style="text-align: left;">')

    return render_template('MS_data.html', MS_table = MS_table)

@app.route('/help')
def help():
    return render_template('help.html')


@app.route('/about')
def about():
    return render_template('about.html')






if __name__ == '__main__':
    app.run(debug=True)
