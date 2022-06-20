import sys
import colorama
colorama.init(autoreset=True, strip=False)

class Main:
    def __init__(self, DEBUG = False, VERBOSE = False):
        self.DEBUG = DEBUG
        self.VERBOSE = VERBOSE
        self.file_path = sys.argv[1]
        
        truth_interval_dict = None
        if len(sys.argv) > 2 and sys.argv[2] != "--verbose":
            self.truth_path = sys.argv[2]
            with open(self.truth_path) as f:
                content = f.read().strip()
            truth_interval_dict = Utils.content_to_dict(content)
            
        sequence = SequenceReader(self.file_path).read_sequence()
        algorithm = ChouFasmanAlgorithm(sequence, truth_interval_dict, self.DEBUG, self.VERBOSE)
        algorithm.compute_structure()

class Utils:
    @staticmethod
    def content_to_dict(content):
        l = [i.strip() for i in content.splitlines() if i.strip()]
        d = {}
        for i in l:
            splitted = i.split()
            d[(int(splitted[1])-1, int(splitted[2]) - 1)] = i[0][0].replace("S", "E")
        return d
    
    @staticmethod
    def generate_position_dict(d, length):
        result = {}
        sorted_keys = sorted(d)
        i = 0
        for interval in sorted_keys:
            ll, ul = interval
            if i < ll:
                for y in range(i, ll):
                    result[y] = 'N'
            for y in range(ll, ul + 1):
                result[y] = d[interval]
            i = ul + 1
        if i < length:
            for y in range(i, length):
                result[y] = 'N'
        return result
    
    @staticmethod
    def generate_text_from_interval_dict(interval_dict):
        string = ""
        sorted_keys = sorted(interval_dict)
        conversion_dict = {"H": "HELIX", "E": "STRAND", "T": "TURN"}
        for key in sorted_keys:
            string += conversion_dict[interval_dict[key]] + "\t" + str(int(key[0]) + 1) + "\t" + str(int(key[1]) + 1) + "\n"
        
        return string
    
    @staticmethod
    def count_for_confusion_matrix(truth_dict, prediction_dict, truth_key, prediction_key):
        start = min(truth_dict.keys())
        end = max(truth_dict.keys())
        counter = 0

        for i in range(start, end + 1):
            if prediction_dict[i] == prediction_key and truth_dict[i] == truth_key:
                counter += 1
        
        return counter
    
    @staticmethod
    def count_individual_confusion_statistics(truth_dict, prediction_dict, key):
        start = min(truth_dict.keys())
        end = max(truth_dict.keys())

        true_positive = 0
        true_negative = 0
        false_positive = 0
        false_negative = 0

        for i in range(start, end + 1):
            if truth_dict[i] == key and prediction_dict[i] == key:
                true_positive += 1
            if truth_dict[i] != key and prediction_dict[i] != key:
                true_negative += 1
            if truth_dict[i] != key and prediction_dict[i] == key:
                false_positive += 1
            if truth_dict[i] == key and prediction_dict[i] != key:
                false_negative += 1
        return true_positive, true_negative, false_positive, false_negative

class ChouFasmanAlgorithm:
    def __init__(self, sequence, truth_interval_dict, DEBUG = False, VERBOSE = False):
        self.INDENTATION = " "*4
        self.DEBUG = DEBUG
        self.VERBOSE = VERBOSE
        self.truth_interval_dict = truth_interval_dict
        self.amino_acid_table = {
            'A': {'Name': 'Alanine', 'P(a)': '1.42', 'P(b)': '0.83', 'P(turn)': '0.66', 'f(i)': '0.06', 'f(i+1)': '0.076', 'f(i+2)': '0.035', 'f(i+3)': '0.058'},
            'R': {'Name': 'Arginine', 'P(a)': '0.98', 'P(b)': '0.93', 'P(turn)': '0.95', 'f(i)': '0.070', 'f(i+1)': '0.106', 'f(i+2)': '0.099', 'f(i+3)': '0.085'},
            'D': {'Name': 'Aspartic Acid', 'P(a)': '1.01', 'P(b)': '0.54', 'P(turn)': '1.46', 'f(i)': '0.147', 'f(i+1)': '0.110', 'f(i+2)': '0.179', 'f(i+3)': '0.081'},
            'N': {'Name': 'Asparagine', 'P(a)': '0.67', 'P(b)': '0.89', 'P(turn)': '1.56', 'f(i)': '0.161', 'f(i+1)': '0.083', 'f(i+2)': '0.191', 'f(i+3)': '0.091'},
            'C': {'Name': 'Cysteine', 'P(a)': '0.70', 'P(b)': '1.19', 'P(turn)': '1.19', 'f(i)': '0.149', 'f(i+1)': '0.050', 'f(i+2)': '0.117', 'f(i+3)': '0.128'}, 
            'E': {'Name': 'Glutamic Acid', 'P(a)': '1.39', 'P(b)': '1.17', 'P(turn)': '0.74', 'f(i)': '0.056', 'f(i+1)': '0.060', 'f(i+2)': '0.077', 'f(i+3)': '0.064'}, 
            'Q': {'Name': 'Glutamine', 'P(a)': '1.11', 'P(b)': '1.10', 'P(turn)': '0.98', 'f(i)': '0.074', 'f(i+1)': '0.098', 'f(i+2)': '0.037', 'f(i+3)': '0.098'}, 
            'G': {'Name': 'Glycine', 'P(a)': '0.57', 'P(b)': '0.75', 'P(turn)': '1.56', 'f(i)': '0.102', 'f(i+1)': '0.085', 'f(i+2)': '0.190', 'f(i+3)': '0.152'}, 
            'H': {'Name': 'Histidine', 'P(a)': '1.00', 'P(b)': '0.87', 'P(turn)': '0.95', 'f(i)': '0.140', 'f(i+1)': '0.047', 'f(i+2)': '0.093', 'f(i+3)': '0.054'}, 
            'I': {'Name': 'Isoleucine', 'P(a)': '1.08', 'P(b)': '1.60', 'P(turn)': '0.47', 'f(i)': '0.043', 'f(i+1)': '0.034', 'f(i+2)': '0.013', 'f(i+3)': '0.056'}, 
            'L': {'Name': 'Leucine', 'P(a)': '1.41', 'P(b)': '1.30', 'P(turn)': '0.59', 'f(i)': '0.061', 'f(i+1)': '0.025', 'f(i+2)': '0.036', 'f(i+3)': '0.070'}, 
            'K': {'Name': 'Lysine', 'P(a)': '1.14', 'P(b)': '0.74', 'P(turn)': '1.01', 'f(i)': '0.055', 'f(i+1)': '0.115', 'f(i+2)': '0.072', 'f(i+3)': '0.095'}, 
            'M': {'Name': 'Methionine', 'P(a)': '1.45', 'P(b)': '1.05', 'P(turn)': '0.60', 'f(i)': '0.068', 'f(i+1)': '0.082', 'f(i+2)': '0.014', 'f(i+3)': '0.055'}, 
            'F': {'Name': 'Phenylalanine', 'P(a)': '1.13', 'P(b)': '1.38', 'P(turn)': '0.60', 'f(i)': '0.059', 'f(i+1)': '0.041', 'f(i+2)': '0.065', 'f(i+3)': '0.065'}, 
            'P': {'Name': 'Proline', 'P(a)': '0.57', 'P(b)': '0.55', 'P(turn)': '1.52', 'f(i)': '0.102', 'f(i+1)': '0.301', 'f(i+2)': '0.034', 'f(i+3)': '0.068'}, 
            'S': {'Name': 'Serine', 'P(a)': '0.77', 'P(b)': '0.75', 'P(turn)': '1.43', 'f(i)': '0.120', 'f(i+1)': '0.139', 'f(i+2)': '0.125', 'f(i+3)': '0.106'}, 
            'T': {'Name': 'Threonine', 'P(a)': '0.83', 'P(b)': '1.19', 'P(turn)': '0.96', 'f(i)': '0.086', 'f(i+1)': '0.108', 'f(i+2)': '0.065', 'f(i+3)': '0.079'}, 
            'W': {'Name': 'Tryptophan', 'P(a)': '1.08', 'P(b)': '1.37', 'P(turn)': '0.96', 'f(i)': '0.077', 'f(i+1)': '0.013', 'f(i+2)': '0.064', 'f(i+3)': '0.167'}, 
            'Y': {'Name': 'Tyrosine', 'P(a)': '0.69', 'P(b)': '1.47', 'P(turn)': '1.14', 'f(i)': '0.082', 'f(i+1)': '0.065', 'f(i+2)': '0.114', 'f(i+3)': '0.125'}, 
            'V': {'Name': 'Valine', 'P(a)': '1.06', 'P(b)': '1.70', 'P(turn)': '0.50', 'f(i)': '0.062', 'f(i+1)': '0.048', 'f(i+2)': '0.028', 'f(i+3)': '0.053'}
        }

        self.sequence = sequence

    def get_propensity_sum_for_region(self, part, propensity_name):
        sum = 0
        for aa in part:
            sum += float(self.amino_acid_table[aa][propensity_name])
        return sum
    
    def get_propensity_average_for_region(self, part, propensity_name):
        return self.get_propensity_sum_for_region(part, propensity_name) / len(part)

    def get_p_bend_for_region(self, part):
        p_bend = 1
        for i in range(len(part)):
            if i == 0:
                string = "f(i)"
            else:
                string = f"f(i+{str(i)})"
            p_bend *= float(self.amino_acid_table[part[i]][string])
        return p_bend


    def detect_overlap(self, tuple1, tuple2):
        return tuple1[1] >= tuple2[0] and tuple2[1] >= tuple1[0]
    
    def merge_overlap(self, tuple1, tuple2):
        return (min(tuple1[0], tuple2[0]), max(tuple1[1], tuple2[1]))
    
    def intersect_overlap(self, tuple1, tuple2):
        return (max(tuple1[0], tuple2[0]), min(tuple1[1], tuple2[1]))

    def resolve_overlaps_inside(self, l):
        i = 0
        while i < len(l) - 1:
            if self.detect_overlap(l[i], l[i + 1]):
                merged = self.merge_overlap(l[i], l[i + 1])
                l.pop(i)
                l[i] = merged
                i = i - 1
            i = i + 1
    
    def resolve_turn_overlaps(self, result_record, alpha_helix_key, beta_sheet_key, turn_key, names_dict):
        for key in (alpha_helix_key, beta_sheet_key):
            if self.DEBUG:
                print(f"\tResolving overlaps between turn and {names_dict[key]} regions")

            overlap = self.detect_overlaps_between_lists(result_record[key], result_record[turn_key])
            self.merge_consecutive_regions(overlap)
            result_record[key] = self.split_from_overlapped_regions(result_record[key], overlap)
            self.merge_consecutive_regions(result_record[key])

        self.delete_less_than_5(result_record, alpha_helix_key, beta_sheet_key, names_dict)

    def resolve_overlaps(self, overlapped_regions, result_record, alpha_helix_key, beta_sheet_key, names_dict):
        for region in overlapped_regions:
            ll, ul = region
            part = self.sequence[ll : ul + 1]
            alpha_sum = self.get_propensity_sum_for_region(part, 'P(a)')
            beta_sum = self.get_propensity_sum_for_region(part, 'P(b)')
            if alpha_sum > beta_sum:
                if self.DEBUG:
                    print(f"\tFor the region {region}, P(a) = {alpha_sum:.4f} > P(b) = {beta_sum:.4f}")
                result_record[alpha_helix_key].append(region)
            else:
                if self.DEBUG:
                    print(f"\tFor the region {region}, P(b) = {beta_sum:.4f} > P(a) = {alpha_sum:.4f}")
                result_record[beta_sheet_key].append(region)
        
        result_record[alpha_helix_key].sort()
        result_record[beta_sheet_key].sort()
        self.merge_consecutive_regions(result_record[alpha_helix_key])
        self.merge_consecutive_regions(result_record[beta_sheet_key])
        self.delete_less_than_5(result_record, alpha_helix_key, beta_sheet_key, names_dict)
    
    def delete_less_than_5(self, result_record, alpha_helix_key, beta_sheet_key, names_dict):
        for key in (alpha_helix_key, beta_sheet_key):
            i = 0
            while i < len(result_record[key]):
                region = result_record[key][i]
                ll, ul = region
                if (ul - ll) + 1 < 5:
                    if self.DEBUG:
                        print(f"\tDeleting {names_dict[key]} assignment for the region {region} since length of the region is {ul - ll + 1} < 5.")
                    result_record[key].remove(region)
                    i = i - 1

                i = i + 1
    
    def merge_consecutive_regions(self, l):
        i = 0
        while i < len(l) - 1:
            if l[i][1] + 1 == l[i + 1][0]:
                merged = (l[i][0], l[i + 1][1])
                l.pop(i)
                l[i] = merged
                i = i - 1
            i = i + 1

    def split_from_overlapped_regions(self, regions, overlapped_regions):
        result = []
        for region in regions:
            result += self.split_region_by(region, overlapped_regions)
        return result

    def split_region_by(self, region, overlapped_regions):
        result_list = []
        to_be_splitted = region
        added = False
        overlap = False
        for overlapped_region in overlapped_regions:
            if self.detect_overlap(to_be_splitted, overlapped_region):
                overlap = True
                result = self.perform_split_between_two_tuples(to_be_splitted, overlapped_region)
                if len(result) < 2:
                    added = False
                    if len(result) == 0: 
                        added = True
                        continue
                    to_be_splitted = result[0]
                else:
                    result_list.append(result[0])
                    added = True
                    to_be_splitted = result[1]
                    overlap = False

        if not added or not overlap:
            result_list.append(to_be_splitted)
        return result_list

    def perform_split_between_two_tuples(self, region, overlapped_region):
        result = []
        if region[0] <= overlapped_region[0] - 1:
            result.append((region[0], overlapped_region[0] - 1))
        if overlapped_region[1] + 1 <= region[1]:
            result.append((overlapped_region[1] + 1, region[1]))
        return result


    def create_position_map(self, region_dict):
        position_map = {}
        for region_name, region_list in region_dict.items():
            for region in region_list:
                position_map[region] = region_name
        return position_map


    def detect_overlaps_between_lists(self, list1, list2):
        result = []
        for i in list1:
            for j in list2:
                if self.detect_overlap(i, j):
                    intersection = self.intersect_overlap(i, j)
                    result.append(intersection)
                    if self.DEBUG:
                        print(f"\tOverlap detected between the regions {i} and {j} as {intersection}")
                    
        self.resolve_overlaps_inside(result)
        return result

    def print_colored_sequence(self, regions, position_map, sequence):
        i = 0
        print("\nColored protein sequence:")
        color_map = {"H": colorama.Fore.GREEN,
                     "E": colorama.Fore.BLUE,
                     "O": colorama.Fore.YELLOW,
                     "T": colorama.Fore.RED}
        for region in regions:
            ll = region[0]
            ul = region[1]
            if i < ll:
                print(sequence[i : ll], end = '')
            print(color_map[position_map[region]] + sequence[ll : ul + 1], end = '')
            i = ul + 1
        
        if i < len(sequence):
            print(sequence[i:], end = '')

        print("\n\nColors:\n")

        print(color_map["H"] + "\tAlpha Helix")
        print(color_map["E"] + "\tBeta Sheet")
        print(color_map["O"] + "\tOverlap")
        print(color_map["T"] + "\tTurn")
        print("\tNo assignment\n")

    def compute_turn_points(self, result_record, turn_key):
        n = 4 # window size for turn computation
        i = 0
        p_bend_threshold = 0.000075
        while i < len(self.sequence) - n + 1:
            window = self.sequence[i : i + n]
            current_interval = [i, i + n - 1]

            if self.DEBUG and self.VERBOSE:
                print(f"\nExamining the window with the index interval [{current_interval[0]}, {current_interval[1]}] for turn detection. The word is {window}:")
                print("-"*50)
            
            p_bend = self.get_p_bend_for_region(window)
            p_avg = self.get_propensity_average_for_region(window, 'P(turn)')
            p_turn = self.get_propensity_sum_for_region(window, 'P(turn)')
            p_alpha = self.get_propensity_sum_for_region(window, 'P(a)')
            p_beta = self.get_propensity_sum_for_region(window, 'P(b)')

            if self.DEBUG and self.VERBOSE:
                print(f"\tP(bend) = {p_bend:.10f}\n\tThe average P(turn) = {p_avg:.4f}\n\tSum of P(turn) = {p_turn:.4f}\n\tSum of P(alpha) = {p_alpha:.4f}\n\tSum of P(beta) = {p_beta:.4f}")

            if p_bend > p_bend_threshold and p_avg > 1.00 and p_turn > p_alpha and p_turn > p_beta:
                if self.DEBUG:
                    if self.VERBOSE:
                        print(f"\tTurn detected at region [{current_interval[0]}, {current_interval[1]}]")
                    else:    
                        print(f"\tTurn detected at region [{current_interval[0]}, {current_interval[1]}]. The word is {window}")

                result_record[turn_key].append((current_interval[0], current_interval[1]))
                i = current_interval[1] + 2
            else:
                i = i + 1
            
            if self.DEBUG and self.VERBOSE:
                print(f"End of the examination for the interval [{current_interval[0]}, {current_interval[1]}]")
                print("-"*50)


    def compute_structure(self):
        n = 6 # sliding window size

        helix_and_sheet = ("alpha_helix", "beta_sheet")     
        alpha_helix, beta_sheet = helix_and_sheet
        turn = "turn"
        result_record = {alpha_helix: [], beta_sheet: [], turn: []}
        names = {alpha_helix: "alpha helix", beta_sheet: "beta sheet", turn: "turn"}
        propensity_names = {alpha_helix: "P(a)", beta_sheet: "P(b)"}

        directions = ("right", "left")

        direction_index_pointer_initialization_functions = {directions[0]: lambda i, n, m: i + n - 1 - 2,
                                                            directions[1]: lambda i, n, m: i + m - 2}

        direction_while_boolean_functions = {directions[0]: lambda j, length, m: j < length - m + 1,
                                             directions[1]: lambda k, length, m: k >= m - 1}
        
        direction_index_interval_functions = {directions[0]: lambda j, m: (j, j + m - 1),
                                              directions[1]: lambda k, m: (k - m + 1, k)}

        direction_index_pointer_new_value_functions = {directions[0]: lambda j: j + 1,
                                                       directions[1]: lambda k: k - 1}
        
        direction_index_pointer_scale_boolean_functions = {directions[0]: lambda j, i, n: j < i + n,
                                                           directions[1]: lambda k, i, n: k > i - 1}

        direction_index_pointer_scale_update_functions = {directions[0]: lambda i, n: i + n,
                                                          directions[1]: lambda i, n: i - 1}
        
        
        direction_index_pointer_scale_str_functions = {
            directions[0]: (lambda j, i, n: f"j({j}) was less than i + n({i + n}).",
                            lambda j, i, n: f"Now, j is set back to i + n({i + n})."),
            directions[1]: (lambda k, i, n: f"k({k}) was greater than i - 1({i - 1}).",
                            lambda k, i, n: f"Now, k is set back to i - 1 ({i - 1}).")
        }

        direction_interval_update_functions = {directions[0]: lambda interval, index_pointer: (interval[0], index_pointer - 1),
                                               directions[1]: lambda interval, index_pointer: (index_pointer + 1, interval[1])}


        for key in helix_and_sheet:
            if self.DEBUG:
                print(f"\n\nStarting {names[key]} computation...")
            
            i = 0
            while i < len(self.sequence) - n + 1:
                window = self.sequence[i : i + n]
                current_interval = [i, i + n - 1]
                if self.DEBUG and self.VERBOSE:
                    print(f"\nExamining the window with the index interval [{i}, {i + n - 1}] for {names[key]} detection. The word is {window}:")
                    print("-"*50)
                less = 0
                greater = 0
                for aa in window:
                    propensity = float(self.amino_acid_table[aa][propensity_names[key]])
                    if self.DEBUG and self.VERBOSE:
                        print(f"{self.INDENTATION*1}Propensity value for {aa} is {propensity}.", end = ' ')
                    if propensity > 1:
                        greater += 1
                        if self.DEBUG and self.VERBOSE:
                            print(f"{propensity_names[key]} > 1. Therefore, greater becomes {greater}")
                    else:
                        less += 1
                        if self.DEBUG and self.VERBOSE:
                            print(f"{propensity_names[key]} <= 1. Therefore, less becomes {less}")

                    
                    
                    if greater >= 4 or less >= 3:
                        if self.DEBUG and self.VERBOSE:
                            if greater >= 4:
                                print(f"{self.INDENTATION*1}There are at least 4 greater. {names[key].capitalize()} is detected.")
                            else:
                                print(f"{self.INDENTATION*1}There are at least 3 less. There cannot be {names[key]} in this region.")
                        break

                if greater >= 4:
                    if self.DEBUG and not self.VERBOSE:
                        print(f"{names[key].capitalize()} is detected at the interval [{current_interval[0]}, {current_interval[1]}].")
                    m = 4 # window size
                    length = len(self.sequence)
                    for direction in directions:
                        if self.DEBUG and self.VERBOSE:
                            print(f"{self.INDENTATION*1}{names[key].capitalize()} detection continues with extension to the {direction}.")
                        index_pointer = direction_index_pointer_initialization_functions[direction](i, n, m)
                        while direction_while_boolean_functions[direction](index_pointer, length, m):
                            ll, ul = direction_index_interval_functions[direction](index_pointer, m)
                            inner_window = self.sequence[ll : ul + 1]
                            if self.DEBUG and self.VERBOSE:
                                print(f"{self.INDENTATION*2}Examining the window with the index interval [{ll}, {ul}]. The word is {inner_window}:")

                            inner_less = 0
                            for inner_aa in inner_window:
                                propensity = float(self.amino_acid_table[inner_aa][propensity_names[key]])
                                if self.DEBUG and self.VERBOSE:
                                    print(f"{self.INDENTATION*3}Propensity value for aa {inner_aa} is {propensity}.", end = ' ')
                                if propensity < 1:
                                    inner_less += 1
                                    if self.DEBUG and self.VERBOSE:
                                        print(f"inner_less is incremented and became {inner_less}")
                                else:
                                    if self.DEBUG and self.VERBOSE:
                                        print(f"There is an aa found to have more propensity than 1. Abort the for loop.")
                                    break

                            if inner_less == 4:
                                if self.DEBUG and self.VERBOSE:
                                    print(f"{self.INDENTATION*3}inner_less is 4. The search is done.")
                                break
                            if self.DEBUG and self.VERBOSE:
                                print(f"{self.INDENTATION*3}inner_less is less than 4. The search continues.")
                            index_pointer = direction_index_pointer_new_value_functions[direction](index_pointer)
                        
                        if direction_index_pointer_scale_boolean_functions[direction](index_pointer, i, n):
                            if self.DEBUG and self.VERBOSE:
                                print(direction_index_pointer_scale_str_functions[direction][0](index_pointer, i, n))
                            index_pointer = direction_index_pointer_scale_update_functions[direction](i, n)
                            if self.DEBUG and self.VERBOSE:
                                print(direction_index_pointer_scale_str_functions[direction][1](index_pointer, i, n))
                        
                        current_interval = direction_interval_update_functions[direction](current_interval, index_pointer)

                        if self.DEBUG:
                            if self.VERBOSE:
                                print(f"{self.INDENTATION*1}Extension to the {direction} is over. The interval after the extension is [{current_interval[0]}, {current_interval[1]}]")
                            else:
                                print(f"{self.INDENTATION*1}The interval after the extension to the {direction} is [{current_interval[0]}, {current_interval[1]}]")

                    if self.DEBUG and self.VERBOSE:
                        print(f"End of the examination for the interval [{i}, {i + n - 1}]")
                        print("-"*50)
                    result_record[key].append(current_interval)
                    i = current_interval[1] + 1
                else:
                    if self.DEBUG and self.VERBOSE:
                        print(f"End of the examination for the interval [{i}, {i + n - 1}]")
                        print("-"*50)

                    i = i + 1    
            self.resolve_overlaps_inside(result_record[key])

        if self.DEBUG:
            print("Starting detection of overlaps...")
        overlap = self.detect_overlaps_between_lists(result_record["alpha_helix"], result_record["beta_sheet"])
        self.merge_consecutive_regions(overlap)


        for key in helix_and_sheet:
            result_record[key] = self.split_from_overlapped_regions(result_record[key], overlap)
            self.merge_consecutive_regions(result_record[key])
        

        regions = sorted(result_record["alpha_helix"] + result_record["beta_sheet"] + overlap)

        position_map = self.create_position_map({
            "H": result_record["alpha_helix"],
            "E": result_record["beta_sheet"],
            "O": overlap
        })

        print("\nAlpha helix and beta sheet regions and overlaps between these regions\n" + "-"*75)
        self.print_colored_sequence(regions, position_map, self.sequence)

        if self.DEBUG:
            print("Resolving overlaps...")

        self.resolve_overlaps(overlap, result_record, alpha_helix, beta_sheet, names)
        regions = sorted(result_record[alpha_helix] + result_record[beta_sheet])
        position_map = self.create_position_map({
            "H": result_record["alpha_helix"],
            "E": result_record["beta_sheet"]
        })
        
        print("\nAlpha helix and beta sheet regions after resolving overlaps\n" + "-"*75)
        self.print_colored_sequence(regions, position_map, self.sequence)

        if self.DEBUG:
            print("Starting turn computation...")
        self.compute_turn_points(result_record, "turn")
        self.merge_consecutive_regions(result_record["turn"])

        if self.DEBUG:
            print("Starting detection of turn overlaps...")
        
        regions = []
        alpha_overlap = self.detect_overlaps_between_lists(result_record[alpha_helix], result_record["turn"])
        self.merge_consecutive_regions(alpha_overlap)
        regions += alpha_overlap

        beta_overlap = self.detect_overlaps_between_lists(result_record[beta_sheet], result_record["turn"])
        self.merge_consecutive_regions(beta_overlap)

        result_record[alpha_helix] = self.split_from_overlapped_regions(result_record[alpha_helix], alpha_overlap)
        self.merge_consecutive_regions(result_record[alpha_helix])

        result_record[beta_sheet] = self.split_from_overlapped_regions(result_record[beta_sheet], beta_overlap)
        self.merge_consecutive_regions(result_record[beta_sheet])

        overlap = sorted(alpha_overlap + beta_overlap)
        self.resolve_overlaps_inside(overlap)
        self.merge_consecutive_regions(overlap)
        regions = sorted(result_record["alpha_helix"] + result_record["beta_sheet"] + overlap)

        position_map = self.create_position_map({
            "H": result_record["alpha_helix"],
            "E": result_record["beta_sheet"],
            "O": overlap
        })
        
        print("\nOverlaps between turn regions and other regions after computing turns\n" + "-"*75)
        self.print_colored_sequence(regions, position_map, self.sequence)

        if self.DEBUG:
            print("Deleting alpha helix and beta sheet assignments where the region is short than 5 residues:")
        self.delete_less_than_5(result_record, alpha_helix, beta_sheet, names)
        regions = sorted(result_record[alpha_helix] + result_record[beta_sheet] + result_record["turn"])
        position_map = self.create_position_map({
            "H": result_record["alpha_helix"],
            "E": result_record["beta_sheet"],
            "T": result_record["turn"]
        })
        
        print("\nFinal secondary structure prediction after resolving turn overlaps\n" + "-"*75)
        self.print_colored_sequence(regions, position_map, self.sequence)

        print(Utils.generate_text_from_interval_dict(position_map))
        
        if self.truth_interval_dict:
            truth_dict = Utils.generate_position_dict(self.truth_interval_dict, len(self.sequence))
            prediction_dict = Utils.generate_position_dict(position_map, len(self.sequence))

            print("3x3 confusion matrix computations:")
            print("True".ljust(10), "Predicted".ljust(10), "Count".ljust(10))
            for key_i in "HET":
                for key_j in "HET":
                    print (key_i.ljust(10), key_j.ljust(10), str(Utils.count_for_confusion_matrix(truth_dict, prediction_dict, key_i, key_j)).ljust(10))
            
            print("Individual confusion matrix computations:")
            for key in "HET":
                print(f"Individual confusion matrix computations for {key}:")
                print("TP".ljust(10), "TN".ljust(10), "FP".ljust(10), "FN".ljust(10))
                tp, tn, fp, fn = Utils.count_individual_confusion_statistics(truth_dict, prediction_dict, key)
                print(str(tp).ljust(10), str(tn).ljust(10), str(fp).ljust(10), str(fn).ljust(10))

    
class SequenceReader:
    def __init__(self, file_path):
        self.file_path = file_path

    def set_file_path(self, file_path):
        self.file_path = file_path
    
    def get_file_path(self):
        return self.file_path

    def read_sequence(self):
        with open(self.file_path) as f:
            lines = f.read().strip().splitlines()

        sequence = None

        for line in lines:
            if not (line.startswith(">") or line.startswith("#")):
                sequence = line
                break

        return sequence.upper()

if __name__ == '__main__':
    VERBOSE = False
    DEBUG = False
    if sys.argv[-1] == "--verbose":
        DEBUG = True
        VERBOSE = True

    main = Main(DEBUG, VERBOSE = VERBOSE)
