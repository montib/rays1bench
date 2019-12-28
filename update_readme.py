import os
import glob
import shutil
import fileinput


def replace_text_in_file(file_to_search, text_to_search, results):

    with fileinput.FileInput(file_to_search, inplace=True) as file:
        for line in file:
            print(line.replace(text_to_search, results), end='')


def parse(dirs, size, add_table_header=False):
    baseline_mrays = -1.0

    if add_table_header:
        results = '\n|scene|version|time|total rays|performance|speedup|\n|--|--|--|--|--|--|\n'
    else:
        results = ''

    for i in dirs:
        os.chdir(i)
        # print('parsing: ' + i + '/' + str(size))

        f = open('out_' + size + '.txt','rt')
        results += '|' + size + '|' + i[4:] + ': '
        line = f.readline()
        results += line
        tokens = line.split('|')  # 123.45 mrays/s
        mrays = float(tokens[3].split()[0])   # 123.45

        if baseline_mrays < 0:
            baseline_mrays = mrays

        mrays_str = '{:.2f}'.format(mrays / baseline_mrays)
        if i == dirs[-1]:
            mrays_str = '**' + mrays_str + '**'   # last value is in bold
        results += mrays_str + '|\n'
        f.close

        os.chdir('../..')

    # print('RESULTS = \n' + results + '\n')

    return results


dirs = ['src/step1', 'src/step2', 'src/step3', 'src/step4', 'src/step5', 'src/step6', 'src/step7', 'src/step8', 'src/step9', 'src/step10', 'src/step11', 'src/step12', 'src/step13']


# all results
shutil.copyfile('README_template.md', 'README.md')
replace_text_in_file('README.md', '__RESULTS_MSVC2019_X64__', parse(dirs, 'large', True) + parse(dirs, 'medium', True) + parse(dirs, 'small', True))

# step1
steps_to_cmp = ['src/step1']
replace_text_in_file('README.md', '__RESULTS_MSVC2019_X64_step1__', parse(steps_to_cmp, 'large', True) + parse(steps_to_cmp, 'medium') + parse(steps_to_cmp, 'small'))

# stepPREV -> stepCURRENT

prev = dirs[0]
for i in dirs[1:]:
    step_name = i[4:]
    to_replace = '__RESULTS_MSVC2019_X64_' + step_name + '__'
    steps_to_cmp = [prev, i]
    replace_text_in_file('README.md', to_replace, parse(steps_to_cmp, 'large', True) + parse(steps_to_cmp, 'medium') + parse(steps_to_cmp, 'small'))

    prev = i
