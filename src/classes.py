from os import path
from json import load

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align

# Fix AA compatibility


class Binder:
    def __init__(self, seq_rec: SeqRecord) -> None:

        self.mode = self._select_mode(seq_rec)
        if self.mode == 'nt':
            self.mm_lmt = 7
        elif self.mode == 'aa':
            self.mm_lmt = 1
        self.hq_lmt = 40

        if isinstance(seq_rec._seq, str):
            seq_rec._seq = Seq(seq_rec._seq)

        for k, v in seq_rec.__dict__.items():
            if not k.startswith('_'):
                self.__setattr__(k, v)

        seq_extract = self._extract_seq(seq_rec)

        # Check reverse complement nt sequence if forward failed
        if seq_extract is None and self.mode == 'nt':
            seq_rec = seq_rec[::-1]  # Reverse to change phred if present
            seq_rec._seq = seq_rec._seq.complement()
            seq_extract = self._extract_seq(seq_rec)
        
        # Sequence recognized and extracted
        if seq_extract is not None:  
            self.seq = seq_extract._seq
            self.phred = seq_extract._per_letter_annotations.get('phred_quality', [])         
        else:
            self.seq = Seq('-')
            self.phred = []

        self.rfs = self._extract_rfs()
        self.quality = self._get_quality_results()

    @classmethod
    def seq(cls, seq: str):
        return cls(SeqRecord(seq))

    def _select_mode(self, seq_rec: SeqRecord) -> str:
        seq_set = set(seq_rec)

        def _isnucleotide(seq_set) -> bool:
            nucleotides = set('ATGCU')
            return seq_set.issubset(nucleotides)

        def _isaminoacid(seq_set) -> bool:
            aminoacids = set('ARNDCQEGHILKMNFPOSUTWYVBZXJ*')
            return seq_set.issubset(aminoacids)

        if _isnucleotide(seq_set):
            return 'nt'
        elif _isaminoacid(seq_set):
            return 'aa'
        else:
            raise ValueError()

    def _load_rfs(self) -> dict:
        # TODO fix !!!
        def _translate_to_aa(rf_dict):
            for rf, params in rf_dict.items():
                for key, value in params.items():
                    if isinstance(value, str):
                        rf_dict[rf][key] = str(Seq(value).translate())
                    elif isinstance(value, int):
                        rf_dict[rf][key] = int(value / 3)
            return rf_dict
        
        file_dir = path.dirname(path.abspath(__file__))
        file_name = 'randomized_fragments.json'
        file_path = path.join(file_dir, file_name)
        
        with open(file_path, 'r') as file:
            rf_dict = load(file)

        rf_dict = rf_dict.get(self.lib)

        if self.mode == 'aa':
            return _translate_to_aa(rf_dict)
        
        return rf_dict

    def _extract_seq(self, seq_rec: SeqRecord) -> SeqRecord:

        # Remove non-letter characters except * (STOP)
        # def seq_cleanup(seq):
        #     chars = ascii_letters + '*'
        #     return ''.join([x for x in seq if x in chars])

        # seq = seq_cleanup(seq)

        # TODO Extend tag lists/objects and create algorithm
        pelB = Seq('GCGGCCCAGCCGGCCATGGCG')
        his_tag = Seq('GGCCCGGGAGGCCAACACCATCACCACCATCAT')
        myc_tag = Seq('GAACAAAAACTCATCTCAGAAGAGGATCTG')

        if self.mode == 'aa':
            pelB = pelB.translate()
            his_tag = his_tag.translate()
            myc_tag = myc_tag.translate()

        pelB_i = self._get_positions(pelB, seq_rec._seq)[1]
        his_tag_i = self._get_positions(his_tag, seq_rec._seq)[0]
        myc_tag_i = self._get_positions(myc_tag, seq_rec._seq)[0]

        end_i = max(his_tag_i, myc_tag_i)

        # Add length if sequence found
        if end_i > -1:
            end_i += len(his_tag)  # TODO myc_tag

        # Extract binder sequence
        seq_extract = seq_rec[pelB_i:end_i]

        # If any sequence extracted
        if seq_extract._seq:
            return seq_extract
        else:
            return None

    def _get_positions(self, target, seq = None) -> tuple:
        
        if seq is None:
            seq = self.seq

        aligner = Align.PairwiseAligner()
        # aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
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

        alignments = aligner.align(seq, target)
        target_alignment = str(alignments[0]).split()[1]
        
        # mm_condition = target_alignment.count('.') <= self.mm_lmt
        # len_condition = len(target_alignment.strip('-')) == len(target)

        # if len_condition and mm_condition:
        #     try:
        #         target_i = target_alignment.replace('.', '|').index('|')
        #         print(alignments.alignment.aligned, target_i)
        #         return target_i
        #     except ValueError:
        #         return -1
        # return -1

        try:
            alignments = aligner.align(str(seq), target)
            target_alignment = str(alignments[0]).split()
            aligned_i = alignments.alignment.aligned[0][0]

            mm_condition = target_alignment.count('.') <= self.mm_lmt
            len_condition = len(range(*aligned_i)) == len(target)

            if len_condition and mm_condition:
                return aligned_i
            return (-1, -1)
        except AttributeError:
            return (-1, -1)
        except IndexError:
            return (-1, -1)

    # TODO translate cannot be properly used | always auto translate?
    def _extract_rfs(self) -> list[Seq]:
        self.lib = self._lib_recognition()
        self.rf_dict = self._load_rfs()

        if self.seq is not None:

            # TODO load exact frame

            rf_extracted = []
            last_pos = (0, 0)

            # Iterate over randomized fragments
            for rf_data in self.rf_dict.values():
                rf = RF(rf_data, self.mode)
                extract_seq, *pos = rf.extract(self.seq)

                if self.mode == 'nt':
                    extract_seq = extract_seq.translate()

                if pos[0] > last_pos[1]:
                    rf_extracted.append(extract_seq)
                    last_pos = pos
                else:
                    rf_extracted.append(Seq('-'))

            return rf_extracted
            
        else:
            return [Seq('-')] * len(self.rf_dict.values())

    # def translate_rfs(self):
    #     if self.mode == 'nt':
    #         self.rfs = [rf.translate() for rf in self.rfs]

    # TODO refactor
    def _lib_recognition(self) -> str:

        lib_patterns = {
            'nt':
                [
                    {
                        'anchor': 'GCCCAGCCGGCCATG',
                        'i_start': 31,
                        'i_stop': 34,
                        'patterns':
                            {
                            'GCT': ['SH_VH3-23_Vk1-39',
                                    'PureLibra'],
                            'GCC': ['AI_VH3-23_Vk3-20',
                                    'AI_VH3-23_Vk1-39'],
                            'AGG': ['AI_VH1-69_Vk3-20',
                                    'AI_VH1-69_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'GCCCAGCCGGCCATG',
                        'i_start': 68,
                        'i_stop': 72,
                        'patterns':
                            {
                            'TGCA': ['SH_VH3-23_Vk1-39',
                                     'PureLibra'],
                            'CGCA': ['AI_VH3-23_Vk3-20',
                                     'AI_VH3-23_Vk1-39'],
                            'CAAG': ['AI_VH1-69_Vk3-20',
                                     'AI_VH1-69_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
                        'i_start': 25,
                        'i_stop': 29,
                        'patterns':
                            {
                                'GTGT': ['AI_VH3-23_Vk3-20',
                                         'AI_VH1-69_Vk3-20'],
                                'CAGA': ['SH_VH3-23_Vk1-39',
                                         'AI_VH3-23_Vk1-39',
                                         'AI_VH1-69_Vk1-39'],
                                'ATTG': ['PureLibra']
                            },
                    },
                    {
                        'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
                        'i_start': 69,
                        'i_stop': 75,
                        'patterns':
                            {
                                'CCGTGT': ['SH_VH3-23_Vk1-39'],
                                'AAGAGC': ['AI_VH3-23_Vk3-20',
                                           'AI_VH1-69_Vk3-20'],
                                'CAGAGT': ['AI_VH3-23_Vk1-39',
                                           'AI_VH1-69_Vk1-39'],
                                'GGAAAG': ['PureLibra']
                            },
                    },
                ],
            'aa': 
                [
                    {
                        'anchor': 'AAQPAMA',
                        'i_start': 8,
                        'i_stop': 13,
                        'patterns':
                            {
                                'AEVKK': ['AI_VH1-69_Vk3-20',
                                          'AI_VH1-69_Vk1-39'],
                                'GGLVQ': ['AI_VH3-23_Vk3-20',
                                          'AI_VH3-23_Vk1-39',
                                          'SH_VH3-23_Vk1-39',
                                          'PureLibra'],
                            },
                    },
                    {
                        'anchor': 'AAQPAMA',
                        'i_start': 15,
                        'i_stop': 23,
                        'patterns':
                            {
                                'SSVKVSCK': ['AI_VH1-69_Vk3-20',
                                             'AI_VH1-69_Vk1-39'],
                                'GSLRLSCA': ['AI_VH3-23_Vk3-20',
                                             'AI_VH3-23_Vk1-39',
                                             'SH_VH3-23_Vk1-39',
                                             'PureLibra'],
                            },
                    },
                    {
                        'anchor': 'GGGGSGGGGSGGGGS',
                        'i_start': 0,
                        'i_stop': 4,
                        'patterns':
                            {
                                'EIVL': ['AI_VH3-23_Vk3-20',
                                        'AI_VH1-69_Vk3-20'],
                                'DIQM': ['AI_VH1-69_Vk1-39'
                                         'AI_VH3-23_Vk1-39',
                                         'SH_VH3-23_Vk1-39'],
                                'TEIV': ['PureLibra'],
                            },
                    },
                    {
                        'anchor': 'GGGGSGGGGSGGGGS',
                        'i_start': 12,
                        'i_stop': 22,
                        'patterns':
                            {
                                'SLSPGERATL': ['AI_VH3-23_Vk3-20',
                                               'AI_VH1-69_Vk3-20',
                                               'PureLibra'],
                                'ASVGDRVTIT': ['AI_VH1-69_Vk1-39'
                                               'AI_VH3-23_Vk1-39',
                                               'SH_VH3-23_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'WYQQKPG',
                        'i_start': 0,
                        'i_stop': 4,
                        'patterns':
                            {
                                'QAPR': ['AI_VH3-23_Vk3-20',
                                        'AI_VH1-69_Vk3-20',
                                        'PureLibra'],
                                'KAPK': ['AI_VH1-69_Vk1-39'
                                         'AI_VH3-23_Vk1-39',
                                         'SH_VH3-23_Vk1-39'],
                            },
                    },
                ],
            }
        

        lib_score = {
                'SH_VH3-23_Vk1-39': 0,
                'AI_VH3-23_Vk3-20': 0,
                'AI_VH3-23_Vk1-39': 0,
                'AI_VH1-69_Vk3-20': 0,
                'AI_VH1-69_Vk1-39': 0,
                'PureLibra': 0
        }

        for elem in lib_patterns[self.mode]:
            try:
                pos = self._get_positions(elem['anchor'])[1]
                frag = self.seq[pos+elem['i_start']:pos+elem['i_stop']]
                rec_lib = elem['patterns'].get(frag)
                if rec_lib:
                    for match in rec_lib:
                        lib_score[match] += 1
            except TypeError as e:
                print(e)

        lib_sorted = sorted(list(lib_score.keys()),
                            key=lambda x: lib_score[x],
                            reverse=True)
        top_lib, second_lib, *rest = lib_sorted


        if lib_score[top_lib] > lib_score[second_lib]:
            return top_lib
        else:
            return 'default'

    def _get_quality_results(self) -> tuple:
        if self.phred:
            phred_hq = round(sum(x > self.hq_lmt for x in self.phred) / len(self.phred), 3)
            phred_min = min(self.phred)
            phred_min_i = self.phred.index(phred_min)
            return (phred_hq, phred_min, f'({phred_min_i})')
        return ()

    def check_unique(self, file) -> bool:
        # unique = True/False
        # return unique
        raise NotImplementedError

        # [binder for binder in binder_list if binder.unique]


class RF:
    def __init__(self, params: dict, mode) -> None:
    
        self.mode = mode
        self.start = params['start']
        self.st_dis = len(self.start) + params['start_dis']
        self.end = params['end']
        self.end_dis = params['end_dis']
        self.limit = params['limit']
        self.end_pos = -1
        if mode == 'nt':
            self.mm_lmt = 7
        elif mode == 'aa':
            self.mm_lmt = 1

    def _get_positions(self, target, seq) -> tuple:
        aligner = Align.PairwiseAligner()
        aligner.query_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        aligner.mismatch_score = -0.5

        alignments = aligner.align(str(seq), target)
        target_alignment = str(alignments[0]).split()

        # mm_condition = target_alignment.count('.') <= self.mm_lmt
        # len_condition = len(target_alignment[1].strip('-')) == len(target)

        # if len_condition and mm_condition:
        #     try:
        #         target_i = target_alignment[1].replace('.', '|').index('|')
        #         return target_i
        #     except ValueError:
        #         return -1
        # return -1

        try:
            alignments = aligner.align(str(seq), target)
            target_alignment = str(alignments[0]).split()
            aligned_i = alignments.alignment.aligned[0][0]

            mm_condition = target_alignment.count('.') <= self.mm_lmt
            len_condition = len(range(*aligned_i)) == len(target)

            if len_condition and mm_condition:
                return aligned_i
            return (-1, -1)
        
        except AttributeError:
            return (-1, -1)
        except IndexError:
            return (-1, -1)

    # Extract RF sequence
    def extract(self, seq) -> tuple:
        self.st_pos = self._get_positions(self.start, seq)[0]
        self.end_pos = self._get_positions(self.end, seq)[0]

        # All positions found -> adjust distance
        if self._check_positions():
            self.rf_start = self.st_pos + self.st_dis
            self.rf_end = self.end_pos + self.end_dis
            seq = Seq(seq[self.rf_start:self.rf_end])

            if self.limit < len(seq):
                seq = Seq('-')

            return (seq, self.rf_start, self.rf_end)

        else:
            return (Seq('-'), 0, 0)

    def _check_positions(self) -> bool:
        conditions = [
            self.st_pos >= 0,
            self.end_pos > self.st_pos,
            self.st_pos + self.st_dis < self.end_pos + self.end_dis,
        ]
        return all(conditions)