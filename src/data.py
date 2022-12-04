from typing import NamedTuple


class Frames(NamedTuple):
    frame_h1: str
    frame_h2: str
    frame_h3: str
    frame_h4: str
    frame_l1: str
    frame_l2: str
    frame_l3: str
    frame_l4: str


# Dict of FRs of all supported libraries
lib_frame_dict = {
    'SH_VH3-23_Vk1-39':
        Frames(
            'GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGA',  # FRH1
            'TGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTC',  # FRH2
            'GACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTGTATTACTGT',  # FRH3
            'TGGGGCCAGGGCACCCTGGTCACCGTCTCCTCA',  # FRH4
            'GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACCGTGTCACCATCACTTGCCGG',  # FRL1
            'TTAAATTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTAT',  # FRL2
            'GGGGTCCCATCACGCTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAA',  # FRL3
            'TTCGGCGGAGGTACCAAGGTGGAGATCAAA',  # FRL4
        ),
    'AI_VH3-23_Vk3-20':
        Frames(
            'GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCCTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGCGCAGCCTCTGGA',  # FRH1
            'TGGGTCCGCCAGGCTCCTGGGAAGGGGCTGGAGTGGGTC',  # FRH2
            'GACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACTCCAAGAACACACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGT',  # FRH3
            'TGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA',  # FRH4
            'GAAATCGTGTTGACCCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGAGAAAGAGCCACCCTCTCCTGCAGG',  # FRL1
            'TTAGCCTGGTACCAACAGAAACCAGGGCAGGCCCCTAGGCTCCTGATCTAT',  # FRL2
            'GGGATCCCTGACAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGGAGCCTGAAGATTTTGCAGTTTATTACTGTCAG',  # FRL3
            'TTTGGCCAAGGGACCAAGGTGGAGATCAAA'  # FRL4
        ),
    'AI_VH3-23_Vk1-39':
        Frames(
            'GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCCTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGCGCAGCCTCTGGA',  # FRH1
            'TGGGTCCGCCAGGCTCCTGGGAAGGGGCTGGAGTGGGTC',  # FRH2
            'GACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACTCCAAGAACACACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGT',  # FRH3
            'TGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA',  # FRH4
            'GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCAGG',  # FRL1
            'TTAAACTGGTACCAACAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTAT',  # FRL2
            'GGGGTCCCTTCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAGCCTGAAGATTTTGCAACTTATTACTGTCAA',  # FRL3
            'TTTGGCCAAGGGACCAAGGTGGAGATCAAA'  # FRL4
        ),
    'AI_VH1-69_Vk3-20':
        Frames(
            'CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCCGTGAAGGTCTCCTGCAAGGCCTCTGGA',  # FRH1
            'TGGGTCCGCCAGGCTCCTGGGCAAGGGCTGGAGTGGATG',  # FRH2
            'CAGAAGTTCCAGGGCCGAGTCACCATCACCGCCGACGAATCCACGAGCACAGCCTATATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCTGTGTATTACTGT',  # FRH3
            'TGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA',  # FRH4
            'GAAATCGTGTTGACCCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGAGAAAGAGCCACCCTCTCCTGCAGG',  # FRL1
            'TTAGCCTGGTACCAACAGAAACCAGGGCAGGCCCCTAGGCTCCTGATCTAT',  # FRL2
            'GGGATCCCTGACAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGGAGCCTGAAGATTTTGCAGTTTATTACTGTCAA',  # FRL3
            'TTCGGCCAAGGGACCAAGGTGGAGATCAAA'  # FRL4
        ),
    'AI_VH1-69_Vk1-39':
        Frames(
            'CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCCGTGAAGGTCTCCTGCAAGGCCTCTGGA',  # FRH1
            'TGGGTCCGCCAGGCTCCTGGGCAAGGGCTGGAGTGGATG',  # FRH2
            'CAGAAGTTCCAGGGCCGAGTCACCATCACCGCCGACGAATCCACGAGCACAGCCTATATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCTGTGTATTACTGT',  # FRH3
            'TGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA',  # FRH4
            'GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCAGG',  # FRL1
            'TTAAACTGGTACCAACAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTAT',  # FRL2
            'GGGGTCCCTTCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAGCCTGAAGATTTTGCAACTTATTACTGTCAA',  # FRL3
            'TTTGGCCAAGGGACCAAGGTGGAGATCAAA'  # FRL4
        ),
    'PureLibra':
        Frames(
            'GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT',  # FRH1
            'GCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCT',  # FRH2
            'TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCG',  # FRH3
            'GACCCCTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCA',  # FRH4
            'ACGGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGT',  # FRL1
            'CTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTAT',  # FRL2
            'AGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGT',  # FRL3
            'TTTGGCCAGGGGACCAAGCTGGAGATCAAACGAGC'  # FRL4
        )
    }
