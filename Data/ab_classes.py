from Bio.Seq import Seq
from Bio import Align
from Data.ab_data import *

class Binder:
    def __init__(self) -> None:
        pass

    def load_params(self, params: dict):
        self.start = params["start"]
        self.st_dis = len(self.start) + params["start_dis"]
        self.end = params["end"]
        self.end_dis = params["end_dis"]
        self.limit = params["limit"]
        self.end_pos = -1

    def load_seq_aa(self, seq):
        
        pelB = Seq("MKYLLPTAAAGLLLLAAQPAMA")
        ha_tag = Seq("GTKVEIK")

        pelB_i = self.get_position(seq_nt=seq, target=pelB)
        ha_tag_i = self.get_position(seq_nt=seq, target=ha_tag)
        
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

        pelB_i = self.get_position(seq_nt=seq, target=pelB)
        his_tag_i = self.get_position(seq_nt=seq, target=his_tag)
        myc_tag_i = self.get_position(seq_nt=seq, target=myc_tag)
        
        end_i = max(his_tag_i, myc_tag_i)
        
        # Add length if sequence found
        if end_i > -1:
            end_i += len(his_tag)

        # Extract binder sequence
        seq_extract = Seq(seq[pelB_i:end_i])
        
        # If any sequence extracted
        if seq_extract:
            self.seq = seq_extract
        else:
            self.seq = None

    def get_position(self, seq_nt, target):
        aligner = Align.PairwiseAligner()
        # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.query_internal_open_gap_score = -100
        aligner.query_internal_extend_gap_score = -100
        aligner.target_internal_open_gap_score= -100
        aligner.target_internal_extend_gap_score = -100
        # aligner.query_left_open_gap_score = 0
        # aligner.query_left_extend_gap_score = 0
        # aligner.query_right_open_gap_score = 0
        # aligner.query_right_extend_gap_score = 0
        aligner.mismatch_score = -0.5
        # aligner.match_score = 1

        target_len = len(target)

        target_alignments = aligner.align(seq_nt, target)
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
        rf_extracted = []
        last_pos = [0, 0]   
        
        for rf_data in rf_dict[mode][self.lib]:
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

    def lib_recognition(self):

        start = "GCCCAGCCGGCCATG"
        linker = "GGTGGCGGTGGATCGGGCGGTGGTGG"

        start_i = self.get_position(self.seq, target=start) + len(start)
        linker_i = self.get_position(self.seq, target=linker) + len(linker)

        lib_pattern_start = [
            ((31, 34), {
                "GCT": ["SH_VH3-23_Vk1-39",
                        "PureLibra"],
                "GCC": ["AI_VH3-23_Vk3-20",
                        "AI_VH3-23_Vk1-39"],
                "AGG": ["AI_VH1-69_Vk3-20",
                        "AI_VH1-69_Vk1-39"]
            }),
            ((68, 72), {
                "TGCA": ["SH_VH3-23_Vk1-39",
                        "PureLibra"],
                "CGCA": ["AI_VH3-23_Vk3-20",
                        "AI_VH3-23_Vk1-39"],
                "CAAG": ["AI_VH1-69_Vk3-20",
                        "AI_VH1-69_Vk1-39"]
            })  
        ]

        lib_pattern_linker = [
            ((25, 29), {
                "CAGA": ["SH_VH3-23_Vk1-39"],
                "GTGT": ["AI_VH3-23_Vk3-20",
                        "AI_VH1-69_Vk3-20"],
                "CAGA": ["AI_VH3-23_Vk1-39",
                        "AI_VH1-69_Vk1-39"],
                "ATTG": ["PureLibra"]
            }),
            ((69, 75), {
                "CCGTGT": ["SH_VH3-23_Vk1-39"],
                "AAGAGC": ["AI_VH3-23_Vk3-20",
                        "AI_VH1-69_Vk3-20"],
                "CAGAGT": ["AI_VH3-23_Vk1-39",
                        "AI_VH1-69_Vk1-39"],
                "GGAAAG": ["PureLibra"]
            })  
        ]

        lib_score = {
            "SH_VH3-23_Vk1-39": 0,
            "AI_VH3-23_Vk3-20": 0,
            "AI_VH3-23_Vk1-39": 0,
            "AI_VH1-69_Vk3-20": 0,
            "AI_VH1-69_Vk1-39": 0,
            "PureLibra": 0
        }
        
        for pattern in lib_pattern_start:
            pos = pattern[0]
            lib = pattern[1]
            frag = self.seq[start_i+pos[0]:start_i+pos[1]]
            rec_lib = lib.get(frag)
            if rec_lib:
                for match in rec_lib:
                    lib_score[match] = lib_score[match] + 1
                    
        for pattern in lib_pattern_linker:
            pos = pattern[0]
            lib = pattern[1]
            frag = self.seq[linker_i+pos[0]:linker_i+pos[1]]
            rec_lib = lib.get(frag)
            if rec_lib:
                for match in rec_lib:
                    lib_score[match] = lib_score[match] + 1
        
        lib_sorted = sorted(list(lib_score.items()), key=lambda x: x[1], reverse=True)
        
        if lib_sorted[0][1] > lib_sorted[1][1]:
            self.lib = lib_sorted[0][0]
        # elif lib_sorted[0][1] > 0:
        #     return lib_sorted[0][0] + " or " + lib_sorted[1][0]
        else:
            self.lib = 0

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
    
    def get_position(self, seq_nt, target):
        aligner = Align.PairwiseAligner()
        aligner.query_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score= -100
        aligner.target_internal_extend_gap_score = -100
        aligner.mismatch_score = -0.5

        target_len = len(target)

        target_alignments = aligner.align(seq_nt, target)
        target_alignment = str(target_alignments[0]).split()
        try:
            target_i = target_alignment[1].replace('.', '|').index("|")
        except ValueError:
            pass       
        
        alignment_len = len(target_alignment[1].strip("-"))
        target_len = len(target)
        
        mismatches = str(target_alignments[0]).count(".")
        
        if alignment_len == target_len and mismatches <= 7:
            return target_i
        else:
            return -1
    
    # Extract RF sequence
    def extract(self, seq_nt):
            self.st_pos = self.get_position(seq_nt, self.start)
            self.end_pos = self.get_position(seq_nt, self.end)
            
            # All positions found -> adjust distance
            if self.check_positions():
                self.rf_start = self.st_pos + self.st_dis
                self.rf_end = self.end_pos + self.end_dis
                seq = Seq(seq_nt[self.rf_start:self.rf_end])
                
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


# if __name__ == "__main__":
#     seq = """AAAAGTAAGTTCCAACCTTACCTTTCCAATTCTTTTCCCTCCCTAAAATAGTCCCGCCCGGCCTCATATCTATATATTATCCCCCATGGTGGGTGGTGATGTTGTTCCCCCCCCCGTCTTTTCTTCTCCCCTTTGTCCCCCCCCCAAAATGGCGGAGGGCAGTGCTACTTTTGTTAACAGAATTAACTACCAAAACCTTTTGGTGCCAAACTGTGGAGGGGGAAATAGAAATCTTACCTACTATCACACCTCCACACGTCAAAGCTGGAGCCCCACCTCTTAATGTTAGATTGATATCACACGTATCGGTATACACTAGAGCTATATGTCCCTTCCTCTGTTTTTCTTCTAGATACATAATTATAGTATCACTGTTCCAATGTATGCACACTCCCCGCGCGAGGAGATGGTGACACGCCTCTCCTACATACACACAGAGGGGAGATGGAGACTGCTTGATCTCATATGCCCGACCCACCTCCTCCACATCCCACCACCCCCCTCTCCCCCCCCACCACGAGGACACGATGCACGCAGTGGCCCCCTGCCCCCCCATAATCTCAAGCCTAACAAAACACCACAATCACGTGCACAGTAATACACTGTCGTGTTCCTCGGCTCTCAGGCTGTTCATTTCCAGATACAGCCTGTTTCTTGAAATTGTCTCTGGACTATGGTAGAACCGCGCCCTCTCATAAGACTCGCACTCTCATAATAGAGACCCCCCACCTCACCACATTATCACGACACTACCCTTCTCCAGCTACCTTCGCCTGTCATCCTGGGATGCACTCCAGCACCATGCACATAGCTAGAAAGACAAGTACCATCACAGTAGCCCTACCACAACGTATTACATCTCACTAGGCACCACCCCCACCACTATACCCAACCCCCCCCCCCCACCACTGCCTAGGACCCTGCCACCTTCCGTCCATACGTCCGGCCTCGGAGCCAGCGAAGTAATCAATCAATTCCAATCGGTCCGGCCCCTTATGGCATATTAAATTTATTTTTCTTTTCTTCCCCTGGTCCTCCTTTGAAATTCAAAACTTCGTTAATCCTTCCTCATCAAATTTCCCTCACCAACAGTAACGGAGCCCGGGCAATTCAAAACACGTTTTAAAATCCCCTATGGGGGGTGCCCTCAATTCCAAGTTGAGCCTTAAAGGGCGAAGGTGTGAGTTGGGCTGTTTGCCGGCTTTTCCCTGGGCCCCCGGGGGATTTCCCGGGGCTGGGGGAGCAAATCTCTTGTTTATTTGTGGCCTTTTCGCGCGCCTTTTTAAATTGCAAATAAGGGGCCGCAACCCGGCTTATCTTGGTGGGAAAAAAAGCGCGCGGGTGATATTGCGCGCGGATGATCATAGGGGGGGGCCGCACTCTATTTAGAGGGGCGCAAAAAAATGCGCCGCCACCCCCAAAAAAAGAAGGCCGCCCCCGCGAAATAGCCGATGATAAATAAATAGATATTTTCTCGTCCGAGGTTGTGTCCCGTCCGCGCGCGCCTGTTACTCTGTCTCCCCAAGAAAAA"""
#     seq = Seq(seq)
#     binder = Binder()
#     binder.load_seq_nt(seq)
#     # binder.load_seq_aa(seq)
#     binder.lib_recognition()
#     print(binder.seq)
#     binder.extract_rfs(mode="nt")
#     print(binder.lib)
#     print(binder.rfs)