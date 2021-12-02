# '*** The script is developed by : Akshay Dedaniya ***'
# 'The script can be used to merge two or more vcf '
# ' The script is design in such a way that it is scalable and computable'
# 'contact me at +919494558500 or email: akshaydedaniya075@gmail.com


import getopt
import sys
import numpy as np


def merging_all_the_common_variant(
    info_tag, format_tag, common_dict_1, common_dict_2, file_name_1, file_name_2
):
    """merging the common variant dict"""

    def changing_tag_in_dict(info_tag, format_tag, dict_file):
        """changing the common tag in the dict file provided, using the tag"""
        for key, value in dict_file.items():
            value_info = value[2].split(";")
            value_format = value[3].split(":")
            value_info = [info_tag.get(item, item) for item in value_info]
            value_format = [format_tag.get(item, item) for item in value_format]
            value[2] = ";".join(value_info)
            value[3] = ":".join(value_format)
        return dict_file

    # for file 1
    new_info_tag = list(map(lambda x: file_name_1 + "_" + x, info_tag))
    changing_info_tag = dict(zip(info_tag, new_info_tag))
    new_format_tag = list(map(lambda x: file_name_1 + "_" + x, format_tag))
    changing_format_tag = dict(zip(format_tag, new_format_tag))
    common_dict_1 = changing_tag_in_dict(
        changing_info_tag, changing_format_tag, common_dict_1
    )

    # for file 2
    new_info_tag = list(map(lambda x: file_name_2 + "_" + x, info_tag))
    changing_info_tag = dict(zip(info_tag, new_info_tag))
    new_format_tag = list(map(lambda x: file_name_2 + "_" + x, format_tag))
    changing_format_tag = dict(zip(format_tag, new_format_tag))
    common_dict_2 = changing_tag_in_dict(
        changing_format_tag, changing_format_tag, common_dict_2
    )

    common_dict = {}
    for key, value in common_dict_1.items():
        new_value = []
        temp_format_dict_1 = {}
        temp_format_dict_2 = {}
        temp_format_detail_list_1 = []
        temp_format_detail_list_2 = []
        temp_format_list_1 = []
        temp_format_list_2 = []

        new_value.insert(0, max(common_dict_1[key][0], common_dict_2[key][0]))
        new_value.insert(1, max(common_dict_1[key][1], common_dict_2[key][1]))
        new_value.insert(
            2,
            (
                ";".join(
                    set(common_dict_1[key][2].split(";"))
                    | set(common_dict_2[key][2].split(";"))
                )
            ),
        )
        temp_format_list = set(common_dict_1[key][3].split(":")) | set(
            common_dict_2[key][3].split(":")
        )
        temp_format_list = [x for x in temp_format_list]

        temp_format_list_1 = common_dict_1[key][3].split(":")
        temp_format_list_2 = common_dict_2[key][3].split(":")
        temp_format_detail_list_1 = common_dict_1[key][4].split(":")
        temp_format_detail_list_2 = common_dict_2[key][4].split(":")

        temp_format_dict_1 = dict(zip(temp_format_list_1, temp_format_detail_list_1))
        temp_format_dict_2 = dict(zip(temp_format_list_2, temp_format_detail_list_2))

        for key1 in temp_format_list:
            if key1 not in temp_format_dict_1:
                temp_format_dict_1[key1] = ""
            if key1 not in temp_format_list_2:
                temp_format_dict_2[key1] = ""

        temp_format_dict_1 = dict(
            sorted(
                temp_format_dict_1.items(),
                key=lambda pair: temp_format_list.index(pair[0]),
            )
        )
        temp_format_dict_2 = dict(
            sorted(
                temp_format_dict_2.items(),
                key=lambda pair: temp_format_list.index(pair[0]),
            )
        )

        new_value.insert(3, (":".join(temp_format_list)))
        new_value.insert(4, (":".join(temp_format_dict_1.values())))
        new_value.insert(5, (":".join(temp_format_dict_2.values())))
        common_dict[key] = new_value

    return common_dict


def identifying_common_info_format_tag(list_vcf_dict_1, list_vcf_dict_2):
    """The function help in indentifing the common tag in INFO and FORMAT between two files"""

    print("Identifying the common INFO and FORMAT tag")
    # Get INFO tag from the list and identifying the common tag
    info_dict_1 = list_vcf_dict_1[2]
    info_dict_2 = list_vcf_dict_2[2]
    info_dict_1 = info_dict_1.split(";")
    info_dict_1 = [w.split("=", 1)[0] for w in info_dict_1]
    info_dict_2 = info_dict_2.split(";")
    info_dict_2 = [w.split("=", 1)[0] for w in info_dict_2]

    common_info_tag = set(info_dict_1) & set(info_dict_2)
    print("The common INFO tag between two files :", list(common_info_tag))

    # Get FORMAT tag from the list and identifying the common tag
    format_dict_1 = list_vcf_dict_1[3]
    format_dict_2 = list_vcf_dict_2[3]
    format_dict_1 = format_dict_1.split(":")
    format_dict_2 = format_dict_2.split(":")

    common_format_tag = set(format_dict_1) & set(format_dict_2)
    print("The common FORMAT tag between two files:", list(common_format_tag))

    return common_info_tag, common_format_tag


def variant_file_to_dict(file_variant):
    """Modifing the variant tabular into dict where key value is "chr_pos_ref_alt"
    and the rest of the information such as info, format are store as value"""

    vcf_dict = {}
    key_modified = []
    file_variant = np.array(file_variant)
    file_variant_chr_pos = file_variant[:, :5]
    file_variant_chr_pos = file_variant_chr_pos.tolist()
    for r in file_variant_chr_pos:
        file_1 = "_".join(r)
        key_modified.append(file_1)

    file_variant_info = file_variant[:, 5:]
    file_variant_info = file_variant_info.tolist()
    vcf_dict = dict(zip(key_modified, file_variant_info))
    return vcf_dict


def header_seperation(vcf_file):
    """This function mainly seperate into two divisions and return the sample_info_header,
    header as array and varaint table as array"""

    header_array = []
    variant_array = []
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("##"):
                header_array.append(line.strip())
            if line.startswith("#C"):
                header = list(line.split("\t"))
                header[-1] = header[-1].strip()
                line = f.read().splitlines()
                variant_array = [item.split("\t") for item in line]
    return header, header_array, variant_array


def reading_file(input_vcf_file1, input_vcf_file2, output_vcf_file):
    """This function is the main function in execution of two files from starting of reading
    the file extracting the common and unique variant and then merging to a final
    output file, It has many sub function to execute the task"""

    # The first file is read
    print("** Reading the " + input_vcf_file1 + " file **")
    header_file1, file_info_header_1, file_variant = header_seperation(input_vcf_file1)
    vcf_file_1_dict = variant_file_to_dict(
        file_variant
    )  # Function for converting the variant list to dict with proper key format
    file_name_1 = input_vcf_file1.split("_")[0]

    # The second file is read
    print("** Reading the " + input_vcf_file2 + " file **")
    header_file2, file_info_header_2, file_variant = header_seperation(input_vcf_file2)
    vcf_file_2_dict = variant_file_to_dict(
        file_variant
    )  # Function for converting the variant list to dict with proper key format
    file_name_2 = input_vcf_file2.split("_")[0]

    # Merging the header info tag of both the files
    for x in file_info_header_2:
        if x not in file_info_header_1:
            file_info_header_1.append(x)

    # Merging the sample header of both the files
    for x in header_file2:
        if x not in header_file1:
            header_file1.append(x)

    # Filtering the dict which are common and unique between the two dicts
    common_dict_keys = set(vcf_file_1_dict) & set(vcf_file_2_dict)
    unique_vcf1 = set(vcf_file_1_dict) - set(vcf_file_2_dict)
    unique_vcf2 = set(vcf_file_2_dict) - set(vcf_file_1_dict)
    unique_dict_1 = {key: vcf_file_1_dict[key] for key in unique_vcf1}
    unique_dict_2 = {key: vcf_file_2_dict[key] for key in unique_vcf2}
    common_dict_1 = {key: vcf_file_1_dict[key] for key in common_dict_keys}
    common_dict_2 = {key: vcf_file_2_dict[key] for key in common_dict_keys}

    # Single line for identification of common tag between two files
    list_vcf_dict_1 = list(common_dict_1.values())[1]
    list_vcf_dict_2 = list(common_dict_2.values())[1]
    info_tag, format_tag = identifying_common_info_format_tag(
        list_vcf_dict_1, list_vcf_dict_2
    )  # Function for identification of common tags
    common_dict = merging_all_the_common_variant(
        info_tag, format_tag, common_dict_1, common_dict_2, file_name_1, file_name_2
    )  # Merging the common dict
    final_variant_dict = {**unique_dict_1, **unique_dict_2, **common_dict}
    final_variant_dict = dict(sorted(final_variant_dict.items(), key=lambda k: k[0]))

    print(
        "The unique number of variant" + input_vcf_file1 + " is :", len(unique_dict_1)
    )
    print(
        "The unique number of variant" + input_vcf_file2 + " is :", len(unique_dict_2)
    )
    print(
        "The common number of variant in "
        + input_vcf_file1
        + "and "
        + input_vcf_file2
        + " is:",
        len(common_dict),
    )
    print("The final number of variant in both the files is :", len(final_variant_dict))

    final_list_variant = []
    for key1, value1 in final_variant_dict.items():
        temp_1 = list(key1.split("_"))
        temp_2 = final_variant_dict[key1]
        final_list_variant.append(list(temp_1 + temp_2))

    file_info_header_1.append(" ".join(header_file1))
    with open(output_vcf_file, "w") as output:
        for row in file_info_header_1:
            output.write(row + "\n")
        for row in final_list_variant:
            s = " ".join(map(str, row))
            output.write(s + "\n")

    print("The process is completed")
    print("The file is saved as " + output_vcf_file)


def main(argv):
    if argv != []:
        try:
            opts, args = getopt.getopt(
                argv, "hf:v:o:", ["help", "file1=", "file2=", "ofile="]
            )
        except getopt.GetoptError as err:
            print(
                "Use the script as :merge_vcf.py -f <vcf_file_1> -v <vcf_file2_name> -o <vcf_file_name>"
            )
            sys.exit(2)
        input_vcf_file1 = ""
        output_vcf_file = ""
        for opt, arg in opts:
            if opt == "-h":
                print(
                    "Try : merge_vcf.py -f <vcf_file_1> -v <vcf_file2_name> -o <vcf_file_name>"
                )
                sys.exit()
            if opt in ("-f", "--file1"):
                input_vcf_file1 = arg
            elif opt in ("-v", "--file2"):
                input_vcf_file2 = arg
            elif opt in ("-o", "--ofile"):
                output_vcf_file = arg
            else:
                assert False
        reading_file(
            input_vcf_file1, input_vcf_file2, output_vcf_file
        )  # function to merge two files
    else:
        print(
            "Use the script as: merge_vcf.py -f <vcf_file_1> -v <vcf_file2_name> -o <vcf_file_name> or use python.py -h"
        )


if __name__ == "__main__":
    print("*** The script is developed by : Akshay Dedaniya ***")
    print("The script can be used to merge two vcf ")
    main(sys.argv[1:])
