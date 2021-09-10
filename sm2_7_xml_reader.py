
import xml.etree.ElementTree as ET
import pandas as pd
import glob

#############################################################################
# Important variables to set for each run
# The DIRECTORY variable sets the working directory or folder to look in for STRmix decon runs
# Paste the path into
DIRECTORY = r'Test2/Decon'
#############################################################################


# Hard coded information about the mixtures used in the experiments
mix_dictionary = {
    'M1-2': ['1A', '1B', '', '', ''],
    'M1-3': ['2A', '2B', '3B', '', ''],
    'M1-4': ['3A', '4B', '5B', '6B', ''],
    'M1-5': ['4A', '7B', '8B', '9B', '10B'],
    'M2-2': ['1A', '22F', '', '', ''],
    'M2-3': ['2A', '2B', '3B', '', ''],
    'M2-4': ['3A', '4B', '5B', '6B', ''],
    'M2-5': ['4A', '7B', '8B', '9B', '10B']
}


def contributors_list(mix_num, contributors):
    mix_code = mix_num + '-' + contributors
    return mix_dictionary[mix_code]


def parse_results_xml_file(file):
    """Takes an XML file and parses it for values needed
    to construct a table of the various LR values. It
    returns a list of the values of interest in the following
    order: Case Number, Sample ID, Case Notes, Seed, Gelman Rubin,
    AfAm Total LR, Af Am HPD LR, Af Am Unified LR, ...
    for the three remaining NIST populations. If the XML file
    is an LR from previous results file it will place a '0.0' in
    the Gelman Rubin column as that data is not carried into the XML file
    from the previous run."""

    tree = ET.parse(file)

    lr_data = []

    if tree.find('./analysisResults') is not None:
        lr_data.append('Interpretation')
    else:
        lr_data.append('LR From Previous')

    case_number = tree.find('.//caseNumber').text
    lr_data.append(case_number)

    sample_id = tree.find('.//sampleId').text
    lr_data.append(sample_id)

    if tree.find('.//caseNotes') is not None:
        case_notes = tree.find('.//caseNotes').text
        lr_data.append(case_notes)
    else:
        lr_data.append('')

    sample_id_parts = sample_id.split('_')
    if sample_id_parts[0] == 'SS':
        lr_data.extend([sample_id_parts[-1], "", "", "", "", sample_id_parts[1]])
    elif sample_id_parts[0] in ['M1', 'M2']:
        mix_number = sample_id.split('_')[0]
        true_contributor_count = str(sample_id.split('_')[1].count('-') + 1)
        contrib_list = contributors_list(mix_number, true_contributor_count)
        lr_data.extend(contrib_list)
        lr_data.append(sample_id_parts[2])
    else:
        lr_data.extend(["", "", "", "", "", ""])

    seed = tree.find('.//seed').text
    lr_data.append(seed)

    contributors = tree.find('.//contributors').text
    lr_data.append(contributors)

    if tree.find('.//totalIterations') is not None:
        total_iterations = tree.find('.//totalIterations').text
        lr_data.append(total_iterations)
    else:
        lr_data.append('0.0')

    if tree.find('.//effectiveSampleSize') is not None:
        effective_sample_size = tree.find('.//effectiveSampleSize').text
        lr_data.append(effective_sample_size)
    else:
        lr_data.append('0.0')

    if tree.find('.//averageLogLikelihood') is not None:
        average_log_likelihood = tree.find('.//averageLogLikelihood').text
        lr_data.append(average_log_likelihood)
    else:
        lr_data.append('0.0')

    if tree.find('.//gelmanRubin') is not None:
        gelman_rubin = tree.find('.//gelmanRubin').text
        lr_data.append(gelman_rubin)
    else:
        lr_data.append('0.0')

    if tree.find('.//lsaeVariance') is not None:
        lsae_variance = tree.find('.//lsaeVariance').text
        lr_data.append(lsae_variance)
    else:
        lr_data.append('0.0')

    if tree.find('.//variance') is not None:
        allele_variance = tree.find('.//variance').text
        lr_data.append(allele_variance)
    else:
        lr_data.append('0.0')

    if tree.find('.//stutterVariances') is not None:
        stutters = {x.attrib['stutter']: x.text for x in tree.findall('.//stutterVariances/variance')}
        lr_data.append(stutters['Back Stutter'])
        lr_data.append(stutters['Forward Stutter'])
        lr_data.append(stutters['Half Back (-2bp) Stutter'])
        lr_data.append(stutters['Double Back (-8bp) Stutter'])
    else:
        lr_data.extend(['0','0','0','0'])

    if tree.find('.//dnaAmount') is not None:
        dna_amounts = [x.text for x in tree.findall('.//dnaAmount')]
        while len(dna_amounts) < 5:
            dna_amounts.append('0')
        lr_data.extend(dna_amounts)
    else:
        lr_data.extend(['0', '0', '0', '0', '0'])

    if tree.find('.//mixtureProportion') is not None:
        proportions = [x.text for x in tree.findall('.//mixtureProportion')]
        while len(proportions) < 5:
            proportions.append('0')
        lr_data.extend(proportions)
    else:
        lr_data.extend(['0', '0', '0', '0', '0'])

    lr_values = {
        'point': [],
        'hpd': [],
        'unified': []
    }
    pops = tree.findall('.//population')

    pop_order = ['afam', 'asian', 'cauc', 'hisp']
    read_order = []
    pop_contributor_orders = {}

    for pop in pops:
        # only taking the population stripping the nist and suffixes from the name
        read_order.append(pop.attrib['name'].split('_')[1])

        lr_total = pop.find('./locusLrTotal').text
        lr_values["point"].append(lr_total)

        hpd_lr = pop.find('.//relation[@type="TOTAL"]/hpdValue').text
        lr_values["hpd"].append(hpd_lr)

        unified_lr = pop.find('.//relation[@type="HD_UNIFIED"]/hpdValue').text
        lr_values["unified"].append(unified_lr)

        if tree.find('.//bestLRContributorPositions') is not None:
            pop_contributor_orders[pop.attrib['name']] = [x.text for x in pop.findall('.//bestLRContributorPositions')]


    # This report an issue with the files not having the populations in the correct order
    if pop_order != read_order and read_order:
        raise Exception(f"The file: {file} has the populations in an incorrect order.")
    elif not read_order:
        # For files with no population data
        lr_data.extend(["", "", "", "", "", "", "", "", "", "", "", ""])
    else:
        for lr_type in lr_values:
            lr_data.extend(lr_values[lr_type])

    if not read_order:
        lr_data.extend(["", "", "", ""])
    else:
        for pop in pop_contributor_orders:
            lr_data.append(pop_contributor_orders[pop])

    # Find errors within the run
    if tree.find('.//missingStutterIssue') is not None:
        missing_stutter_issue = [{'Missing Stutter': x.attrib} for x in tree.findall('.//missingStutterIssue')]
        while len(missing_stutter_issue) < 10:
            missing_stutter_issue.append('')
        lr_data.extend(missing_stutter_issue)
    else:
        lr_data.extend(['', '', '', '', '', '', '', '', '', ''])
    return lr_data


def make_data_frame_and_export(array):
    """Takes an array of parsed XML data and turns it into
    a Pandas Data Frame for export as a CSV file."""

    df = pd.DataFrame(array,
                      columns=['Run Type', 'Case Number', 'Sample ID', 'Case Notes',
                               'True Contr. 1', 'True Contr. 2', 'True Contr. 3',
                               'True Contr. 4', 'True Contr. 5', 'Input DNA', 'Seed',
                               'Contributors', 'Total Iterations', 'Effective Sample Size',
                               'Average Log Likelihood', 'Gelman Rubin',
                               'LSAE variance', 'Allele Variance', 'BS Variance',
                               'FS Variance', 'HBS Variance', 'DBS Variance',
                               'DNA Amount 1', 'DNA Amount 2', 'DNA Amount 3',
                               'DNA Amount 4', 'DNA Amount 5',
                               'Proportion 1', 'Proportion 2', 'Proportion 3',
                               'Proportion 4', 'Proportion 5',
                               'P.LR AfAm', 'P.LR Asian', 'P.LR Cauc', 'P.LR Hisp',
                               'HPD AfAm', 'HPD Asian', 'HPD Cauc', 'HPD Hisp',
                               'Uni AfAm', 'Uni Asian', 'Uni Cauc', 'Uni Hisp',
                               'Best LR Contributor AfAm',
                               'Best LR Contributor Asian',
                               'Best LR Contributor Cauc',
                               'Best LR Contributor Hisp',
                               'Issue 1', 'Issue 2', 'Issue 3', 'Issue 4', 'Issue 5',
                               'Issue 6', 'Issue 7', 'Issue 8', 'Issue 9', 'Issue 10'
                               ])

    df['Seed'] = df['Seed'].astype(int)
    df['Contributors'] = df['Contributors'].astype(int)
    df['Total Iterations'] = df['Total Iterations'].astype(float)
    df['Effective Sample Size'] = df['Effective Sample Size'].astype(float)
    df['Average Log Likelihood'] = df['Average Log Likelihood'].astype(float)
    df['Gelman Rubin'] = df['Gelman Rubin'].astype(float)
    df['DNA Amount 1'] = df['DNA Amount 1'].astype(int)
    df['DNA Amount 2'] = df['DNA Amount 2'].astype(int)
    df['DNA Amount 3'] = df['DNA Amount 3'].astype(int)
    df['DNA Amount 4'] = df['DNA Amount 4'].astype(int)
    df['DNA Amount 5'] = df['DNA Amount 5'].astype(int)
    df['P.LR AfAm'] = df['P.LR AfAm'].astype(float, errors = 'ignore')
    df['P.LR Asian'] = df['P.LR Asian'].astype(float, errors = 'ignore')
    df['P.LR Cauc'] = df['P.LR Cauc'].astype(float, errors = 'ignore')
    df['P.LR Hisp'] = df['P.LR Hisp'].astype(float, errors = 'ignore')
    df['HPD AfAm'] = df['HPD AfAm'].astype(float, errors = 'ignore')
    df['HPD Asian'] = df['HPD Asian'].astype(float, errors = 'ignore')
    df['HPD Cauc'] = df['HPD Cauc'].astype(float, errors = 'ignore')
    df['HPD Hisp'] = df['HPD Hisp'].astype(float, errors = 'ignore')
    df['Uni AfAm'] = df['Uni AfAm'].astype(float, errors = 'ignore')
    df['Uni Asian'] = df['Uni Asian'].astype(float, errors = 'ignore')
    df['Uni Cauc'] = df['Uni Cauc'].astype(float, errors = 'ignore')
    df['Uni Hisp'] = df['Uni Hisp'].astype(float, errors = 'ignore')

    df.to_csv('LR_data.csv', index = False)
    # df.to_csv('LR_data_rounded.csv', float_format='%.5g', index = False)


def main():
    # Using Test path with subfolder names as * wildcard
    # Could specifically refer to file or use * wildcard
    lr_array = []
    for x in glob.glob(fr'{DIRECTORY}\*\results.xml', recursive=True):
        lr_array.append(parse_results_xml_file(x))

    make_data_frame_and_export(lr_array)


if __name__ == '__main__':
    main()
