from os import path
from json import load

from Bio.Seq import Seq
from Bio import Align


class Binder:
    def __init__(self, mode, seq) -> None:
        if mode == 'nt':
            self.load_seq_nt(seq)
        elif mode == 'aa':
            self.load_seq_aa(seq)

    def load_params(self, params: dict):
        self.start = params["start"]
        self.st_dis = len(self.start) + params["start_dis"]
        self.end = params["end"]
        self.end_dis = params["end_dis"]
        self.limit = params["limit"]
        self.end_pos = -1

    def load_rfs(self, mode='nt'):
        file_dir = path.dirname(path.abspath(__file__))
        file_name = 'randomized_fragments.json'
        file_path = path.join(file_dir, file_name)
        with open(file_path, 'r') as file:
            rf_dict = load(file)
        self.rf_dict = rf_dict[mode][self.lib]

    def load_seq_aa(self, seq):

        pelB = Seq("MKYLLPTAAAGLLLLAAQPAMA")
        ha_tag = Seq("GTKVEIK")

        pelB_i = self.get_position(seq=seq, target=pelB)
        ha_tag_i = self.get_position(seq=seq, target=ha_tag)

        # Add length if sequence found
        if ha_tag_i > -1:
            ha_tag_i += len(ha_tag)

        # Extract binder sequence
        if pelB < ha_tag:
            seq_extract = Seq(seq[pelB_i:ha_tag_i])
        else:
            seq_extract = Seq(seq[pelB_i:])

        self.seq = seq_extract

    def load_seq_nt(self, seq):

        # TODO Create tag list and algorithm

        pelB = Seq("GCGGCCCAGCCGGCCATGGCG")
        his_tag = Seq("CACCATCACCACCATCAT")
        myc_tag = Seq("GAACAAAAACTCATCTCA")

        pelB_i = self.get_position(seq=seq, target=pelB)
        his_tag_i = self.get_position(seq=seq, target=his_tag)
        myc_tag_i = self.get_position(seq=seq, target=myc_tag)

        end_i = max(his_tag_i, myc_tag_i)

        # Add length if sequence found
        if end_i > -1:
            end_i += len(his_tag)  # TODO myctag

        # Extract binder sequence
        seq_extract = Seq(seq[pelB_i:end_i])

        # If any sequence extracted
        if seq_extract:
            self.seq = seq_extract
        else:
            self.seq = None

    def get_position(self, target, seq=None):
        if seq is None:
            seq = self.seq
        
        aligner = Align.PairwiseAligner()
        # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.query_internal_open_gap_score = -100
        aligner.query_internal_extend_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        # aligner.query_left_open_gap_score = 0
        # aligner.query_left_extend_gap_score = 0
        # aligner.query_right_open_gap_score = 0
        # aligner.query_right_extend_gap_score = 0
        aligner.mismatch_score = -0.5
        # aligner.match_score = 1

        target_len = len(target)

        target_alignments = aligner.align(seq, target)
        target_alignment = str(target_alignments[0]).split()

        try:
            target_i = target_alignment[1].replace('.', '|').index("|")
        except ValueError:
            pass

        alignment_len = len(target_alignment[1].strip("-"))
        target_len = len(target)

        mismatches = str(target_alignments[0]).count(".")

        # Check if the seq is revcomp
        #
        # if target_i == -1:
        #     seq_revcomp = seq_nt.reverse_complement()
        #     target_i_revcomp = self.get_position(seq_revcomp, target)
        #     if target_i_revcomp != 0:
        #         target_i = target_i_revcomp
        #         self.revcomp = True

        if alignment_len == target_len and mismatches <= 6:
            return target_i
        else:
            return -1

    def extract_rfs(self, mode):
        self.lib_recognition()
        self.load_rfs(mode)

        # TODO load exact frame

        rf_extracted = []
        last_pos = [0, 0]

        for rf_data in self.rf_dict.values():
            rf = RF()
            rf.load_params(rf_data)
            extract, *pos = rf.extract(self.seq)

            if pos[0] > last_pos[1]:
                rf_extracted.append(extract)
                last_pos = pos
            else:
                rf_extracted.append(Seq("-"))

        self.rfs = rf_extracted

    def translate_rfs(self):
        self.rfs = [rf.translate() for rf in self.rfs]

    # TODO refactor
    def lib_recognition(self):

        lib_patterns = [
            {
            'anchor': 'GCCCAGCCGGCCATG',
            'i_start': 31,
            'i_stop': 34,
            'patterns':
                {
                "GCT": ["SH_VH3-23_Vk1-39",
                        "PureLibra"],
                "GCC": ["AI_VH3-23_Vk3-20",
                        "AI_VH3-23_Vk1-39"],
                "AGG": ["AI_VH1-69_Vk3-20",
                        "AI_VH1-69_Vk1-39"]
                },
            },
            {
             'anchor': 'GCCCAGCCGGCCATG',
             'i_start': 68,
             'i_stop': 72,
             'patterns':
                {
                "TGCA": ["SH_VH3-23_Vk1-39",
                         "PureLibra"],
                "CGCA": ["AI_VH3-23_Vk3-20",
                         "AI_VH3-23_Vk1-39"],
                "CAAG": ["AI_VH1-69_Vk3-20",
                         "AI_VH1-69_Vk1-39"]
                },
            },
            {
            'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
            'i_start': 25,
            'i_stop': 29,
            'patterns':
                {
                "GTGT": ["AI_VH3-23_Vk3-20",
                            "AI_VH1-69_Vk3-20"],
                "CAGA": ["SH_VH3-23_Vk1-39",
                            "AI_VH3-23_Vk1-39",
                            "AI_VH1-69_Vk1-39"],
                "ATTG": ["PureLibra"]
                },
            },
            {
            'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
            'i_start': 69,
            'i_stop': 75,
            'patterns':
                {
                "CCGTGT": ["SH_VH3-23_Vk1-39"],
                "AAGAGC": ["AI_VH3-23_Vk3-20",
                            "AI_VH1-69_Vk3-20"],
                "CAGAGT": ["AI_VH3-23_Vk1-39",
                            "AI_VH1-69_Vk1-39"],
                "GGAAAG": ["PureLibra"]
                },
            },
        ]

        lib_score = {
            "SH_VH3-23_Vk1-39": 0,
            "AI_VH3-23_Vk3-20": 0,
            "AI_VH3-23_Vk1-39": 0,
            "AI_VH1-69_Vk3-20": 0,
            "AI_VH1-69_Vk1-39": 0,
            "PureLibra": 0
        }

        for elem in lib_patterns:
            pos = self.get_position(elem['anchor']) + len(elem['anchor'])
            frag = self.seq[pos+elem['i_start']:pos+elem['i_stop']]
            rec_lib = elem['patterns'].get(frag)
            if rec_lib:
                for match in rec_lib:
                    lib_score[match] += 1

        lib_sorted = sorted(list(lib_score.keys()),
                            key=lambda x: lib_score[x],
                            reverse=True)
        top_lib, second_lib, *rest = lib_sorted
        
        if lib_score[top_lib] > lib_score[second_lib]:
            self.lib = top_lib
        else:
            self.lib = "default"


    def check_unique(self, file):
        # unique = True/False
        # return unique
        pass

        # [binder for binder in binder_list if binder.unqiue]


class RF:
    def __init__(self):
        pass

    def load_params(self, params: dict):
        self.start = params["start"]
        self.st_dis = len(self.start) + params["start_dis"]
        self.end = params["end"]
        self.end_dis = params["end_dis"]
        self.limit = params["limit"]
        self.end_pos = -1

    def get_position(self, target, seq=None):
        if seq is None:
            self.seq
        aligner = Align.PairwiseAligner()
        aligner.query_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        aligner.mismatch_score = -0.5

        target_len = len(target)

        target_alignments = aligner.align(seq, target)
        target_alignment = str(target_alignments[0]).split()

        try:
            target_i = target_alignment[1].replace('.', '|').index("|")
        except ValueError as e:
            print(e)
        except Exception as e:
            print(e)

        alignment_len = len(target_alignment[1].strip("-"))
        target_len = len(target)

        mismatches = str(target_alignments[0]).count(".")

        if alignment_len == target_len and mismatches <= 7:
            return target_i
        else:
            return -1

    # Extract RF sequence
    def extract(self, seq):
        self.st_pos = self.get_position(self.start, seq)
        self.end_pos = self.get_position(self.end, seq)

        # All positions found -> adjust distance
        if self.check_positions():
            self.rf_start = self.st_pos + self.st_dis
            self.rf_end = self.end_pos + self.end_dis
            seq = Seq(seq[self.rf_start:self.rf_end])

            if self.limit < len(seq):
                seq = Seq("-")

            return (seq, self.rf_start, self.rf_end)

        else:
            return (Seq("-"), 0, 0)

    def check_positions(self):
        cond_1 = self.st_pos >= 0
        cond_2 = self.end_pos > self.st_pos
        cond_3 = self.st_pos + self.st_dis < self.end_pos + self.end_dis

        return all([cond_1, cond_2, cond_3])


if __name__ == "__main__":
    seq = """CCAAAACGTCTCAGGCAAGAAGGAGCACCGGCATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAGTGCGCGATGAAGCAGTCGATCGACTATTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCAGCGGCGGTGGCGGGTCGACGGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAGAAGGACCGGCGCATGTAGACTTTTGGCCAGGGGACCAAGCTGGAGATCAAACGAGCGGCCGCAGGGGCCGCAGAACAAAAACTCATCTCAGAAGAGGATCTGGGAGACGCGGGTGGCGGCGGTTCTACTGTTGAAAGTTGTTTAGCAAAACCTCATACAGAAAATTCATTTACTAACGTCTGGAAAGACGACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGCTACAGGCGTTGTGGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACATGGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAAGGGGGGTGGCTCTGAGGGTGGCGGTTCTGAAGGTGGGCGGTTCTGAAGGTGGCGATACAAACCTCCTGAGATAGGGATACACCCATTTCCGGGCAAA"""
    seq = Seq(seq)
    binder = Binder()
    binder.load_seq_nt(seq)
    # binder.load_seq_aa(seq)
    binder.lib_recognition()
    print(binder.seq)
    binder.extract_rfs(mode="nt")
    print(binder.lib)
    print(binder.rfs)
