import os
import collections
import Levenshtein as ls


white_dict = {}
# with open("white_list.tsv") as f:
with open("96lib_gDNA_whitelist_01092023.tsv") as f:
    f.readline()
    for line in f:
        line = line.rstrip().split("\t")
        white_dict[line[0]] = {}
        white_dict[line[0]]["bcseq"] = line[1]
        white_dict[line[0]]["count"] = line[7]
        white_dict[line[0]]["freq"] = line[5]
        white_dict[line[0]]["old_name"] = line[9]

count_dict1 = collections.defaultdict(lambda: collections.defaultdict(
    lambda: {"all": 0, "GTG": 0, "ATG": 0, "others": 0}))
with open("./data/GFP_isolation/project_barcode_count_rep1.tsv") as f:
    for line in f:
        line = line.rstrip().split("\t")
        count_dict1[line[0]][line[2][1:-3]]["all"] += int(line[3])
        if line[2][-3:] == "GTG":
            count_dict1[line[0]][line[2][1:-3]]["GTG"] += int(line[3])
        elif line[2][-3:] == "ATG":
            count_dict1[line[0]][line[2][1:-3]]["ATG"] += int(line[3])
        else:
            count_dict1[line[0]][line[2][1:-3]]["others"] += int(line[3])

white_keys = list(white_dict.keys())
count_dict2 = collections.defaultdict(dict)
for old_name in count_dict1:
    # print(old_name)
    for key in white_keys:
        if old_name == white_dict[key]["old_name"]:
            new_name = key
            break
        else:
            pass

    for key in white_keys:
        count_dict2[new_name][key] = {
            "all": 0, "GTG": 0, "ATG": 0, "others": 0}

    for seq1 in count_dict1[old_name]:
        distances = []
        for key in white_keys:
            distances.append(ls.distance(seq1, white_dict[key]["bcseq"][:-3]))

        min_index = distances.index(min(distances))
        sorted_distances = list(sorted(distances))
        # print(distances[min_index])
        if sorted_distances[0] < 3:
            count_dict2[new_name][white_keys[min_index]
                                  ]["all"] += count_dict1[old_name][seq1]["all"]
            if seq1[-3:] == "ATG":
                count_dict2[new_name][white_keys[min_index]
                                      ]["ATG"] += count_dict1[old_name][seq1]["ATG"]
            elif seq1[-3:] == "GTG":
                count_dict2[new_name][white_keys[min_index]
                                      ]["GTG"] += count_dict1[old_name][seq1]["GTG"]
            else:
                count_dict2[new_name][white_keys[min_index]
                                      ]["others"] += count_dict1[old_name][seq1]["others"]

used_keys = list(count_dict2.keys())
used_keys.sort()
white_keys.sort()


all_keys = used_keys + [key for key in white_keys if key not in used_keys]
print("used_barcode", *all_keys, sep="\t")
for key1 in used_keys:
    values = [count_dict2[key1][key2]["all"] for key2 in all_keys]
    print(key1, *values, sep="\t")
