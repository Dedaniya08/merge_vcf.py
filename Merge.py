# '*** The script is developed by : Akshay Dedaniya ***'
# 'The script can be used to merge two or more vcf '
# ' The script is design in such a way that it is scalable and computable'
# 'contact me at +919494558500 or email: akshaydedaniya075@gmail.com

import csv
import getopt
import sys
import dask.dataframe as dd
import pandas as pd


def output_header_modification(file_1, info_tag, format_tag):
    file_header_combined_df = pd.DataFrame([])
    for file_in_list in file_1:
        count = 0
        file_name = file_in_list.split('_')[0]
        all_tag = info_tag + format_tag
        info_tag_1 = list(map(lambda x: r'\b' + x + r'\b', all_tag))
        new_info_tag = list(map(lambda x: file_name + '_' + x, all_tag))
        changing_tag_dict = dict(zip(info_tag_1, new_info_tag))
        for line in open(file_in_list):
            if line.startswith('##'):
                count += 1
            if line.startswith('#C'):
                break

        with open(file_in_list, 'r') as f:
            file_header = [next(f) for x in range(count)]
        file_header = map(lambda s: s.strip(), file_header)
        dataframe_pd = pd.DataFrame(file_header)
        dataframe_pd = dataframe_pd.replace(changing_tag_dict, regex=True)
        file_header_combined_df = file_header_combined_df.append(dataframe_pd)
    file_header_combined_df = file_header_combined_df.drop_duplicates(keep='first').astype(str)
    return file_header_combined_df


def merging_the_table(output_file_final_main_with_count):
    final_table = pd.DataFrame([])
    groups = output_file_final_main_with_count.groupby('CHROM')
    for group in output_file_final_main_with_count['CHROM'].unique():
        chrwise_merging = groups.get_group(group).compute()
        unique_merge = chrwise_merging[chrwise_merging['chr_pos'] == 1]
        unique_merge = unique_merge.drop(['chr_pos'], axis=1)
        unique_merge = unique_merge.set_index('POS')
        common_merge = chrwise_merging[chrwise_merging['chr_pos'] > 1]
        common_merge = common_merge.drop(['chr_pos'], axis=1)
        col_common_merge = list(common_merge.columns.values)
        common_merge_1 = common_merge[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']]
        common_merge_1 = common_merge_1.groupby(['POS']).max().reset_index()
        common_merge_1 = common_merge_1.set_index('POS')
        common_merge_2 = common_merge.drop(['CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)
        common_merge_2 = common_merge_2.set_index('POS').fillna('')
        common_merge_1['INFO'] = common_merge.groupby(['POS'])['INFO'].apply(lambda x: ';'.join(x))
        common_merge_1['FORMAT'] = common_merge.groupby(['POS'])['FORMAT'].apply(lambda x: ':'.join(x))
        for col in common_merge_2:
            common_merge_1[col] = common_merge_2.groupby(['POS'])[col].apply(lambda x: ''.join(x))
        vcf_table_output_file = unique_merge.append(common_merge_1)
        final_table = final_table.append(vcf_table_output_file)
    final_table.reset_index(inplace=True)
    col_name = col_common_merge
    final_table = final_table[col_name]
    return final_table


def count_info_and_variant(file):
    count = 0
    header_columns = ''
    file_header = []
    for line in open(file):
        if line.startswith('##'):
            count += 1
        if line.startswith('#C'):
            header_columns = line[1:].split()
            break

    with open(file, 'r') as f:
        for i in range(count):
            file_header = next(f).strip()

    return count, header_columns


def finding_common_info_format(file_1):
    file_name_list_count = 0
    output_file_final = ''
    for vcf in file_1:
        info_line_count, header_columns = count_info_and_variant(vcf)
        file_name = vcf.split('_')[0]
        variant = dd.read_csv(vcf, skiprows=info_line_count, sep='\t', header=None, comment='#', names=header_columns, dtype='object')
        variant['INFO'] = variant['INFO'] + ';calledBy=' + file_name

        if file_name_list_count < 1:
            output_file_final = variant
        else:
            output_file_final = output_file_final.append(variant)
        file_name_list_count += 1
    output_file_final['chr_pos'] = output_file_final.CHROM + ':' + output_file_final.POS
    variant_count = output_file_final.chr_pos.value_counts()
    output_file_final = output_file_final.set_index('chr_pos')
    variant_file_with_count = output_file_final.join(variant_count, how='inner').compute()
    list_for_common_count = 1
    for i1, g1 in variant_file_with_count.groupby('CHROM'):
        chr_file = g1
        common = chr_file[chr_file['chr_pos'] > 1]
        common = common.replace({'INFO': r'calledBy=.+'}, {'INFO': 'calledBy=Freebayes+VarScan'}, regex=True)
        if list_for_common_count == 1:
            common = common.groupby('CHROM')
            common_tag_info_format = common.get_group(i1).head(file_name_list_count)
            pd_common_tag_info_format = pd.DataFrame(common_tag_info_format)
            for i3 in range(1, file_name_list_count):
                list_of_common_info = pd_common_tag_info_format.loc[:, 'INFO'].str.split(';').tolist()
                list_of_common_info = [j for i in list_of_common_info for j in i]
                common_info_tag_in_common_variant = list(set([x for x in list_of_common_info if list_of_common_info.count(x) > 1]))
                common_info_tag_in_common_variant.remove('calledBy=Freebayes+VarScan')

                list_of_common_format = pd_common_tag_info_format.loc[:, 'FORMAT'].str.split(':').tolist()
                list_of_common_format = [j for i in list_of_common_format for j in i]
                common_format_tag_in_common_variant = list(set([x for x in list_of_common_format if list_of_common_format.count(x) > 1]))
        list_for_common_count += 1
    return common_info_tag_in_common_variant, common_format_tag_in_common_variant


def vcf_merging(file_1, output_file_name):
    info_tag_label_comm, format_tag_label_comm = finding_common_info_format(file_1)
    output_file_final_main = ''
    file_name_list_count_1 = 0
    for vcf in file_1:
        info_line_count, header_columns = count_info_and_variant(vcf)
        file_name = vcf.split('_')[0]
        variant = dd.read_csv(vcf, skiprows=info_line_count, sep='\t', header=None, comment='#', names=header_columns, dtype='object')
        new_format_tag_list = list(map(lambda x: file_name + '_' + x, format_tag_label_comm))
        format_tag_label_comm_1 = list(map(lambda x: r'\b' + x + r'\b', format_tag_label_comm))
        format_dict = dict(zip(format_tag_label_comm_1, new_format_tag_list))
        variant['FORMAT'] = variant['FORMAT'].replace(format_dict, regex=True)
        info_dict = dict(zip(format_tag_label_comm_1, new_format_tag_list))
        variant['INFO'] = variant['INFO'].replace(info_dict, regex=True)
        variant['INFO'] = variant['INFO'] + ';calledBy=' + file_name
        if file_name_list_count_1 < 1:
            output_file_final_main = variant
        else:
            # variant['file_no'] = file_name_list_count_1
            output_file_final_main = output_file_final_main.append(variant)
        file_name_list_count_1 += 1

    file_header = output_header_modification(file_1, info_tag_label_comm, format_tag_label_comm)
    output_file_final_main['chr_pos'] = output_file_final_main.CHROM + ':' + output_file_final_main.POS
    variant_count = output_file_final_main.chr_pos.value_counts()
    output_file_final_main = output_file_final_main.set_index('chr_pos')
    output_file_final_main_with_count = output_file_final_main.join(variant_count, how='inner')
    final_table = merging_the_table(output_file_final_main_with_count)
    final_table.rename(columns={'CHROM': '#CHROM'}, inplace=True)
    file_header.to_csv(output_file_name, header=False, index=False, mode='w', quoting=csv.QUOTE_NONE, escapechar=' ')
    final_table.to_csv(output_file_name, header=True, index=None, sep='\t', mode='a', chunksize=1000)


def main(argv):
    if argv != []:
        try:
            opts, args = getopt.getopt(argv, 'hi:o:', ['help', 'i1file=', 'ofile='])
        except getopt.GetoptError as err:
            print(err)
            print('Use the script as :Merge.py -i <vcf_file_name>,<vcf_file2_name>,... -o <vcf_file_name>')
            sys.exit(2)
        input_vcf_file1 = ''
        output_vcf_file = ''
        for opt, arg in opts:
            if opt == '-h':
                print('Try : Merge.py -i <vcf_file_name>,<vcf_file2_name>,...  -o <vcf_file_name>')
                sys.exit()
            if opt in ('-i', '--ifile'):
                input_vcf_file1 = arg.split(',')
            elif opt in ('-o', '--ofile'):
                output_vcf_file = arg
            else:
                assert False, 'Unhandled Option'
        vcf_merging(input_vcf_file1, output_vcf_file)
    else:
        print('Use the script as: Merge.py -i <vcf_file>,<vcf_file2>,.... -o <vcf_name> or use python.py -h')


if __name__ == '__main__':
    print('*** The script is developed by : Akshay Dedaniya ***')
    print('The script can be used to merge two or more vcf ')
    main(sys.argv[1:])
