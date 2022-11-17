import json
import pathlib

mypath = pathlib.Path(__file__).parent.absolute()

filename = 'reading_frames.json'

rf_dict = {
    'aa': {
        0: {'rf_1':
                {
                'start': 'AQPAM',
                'start_dis': 27,
                'end': 'WVRQA',
                'end_dis': 0,
                'limit': 15
                },
            'rf_2':
                {
                'start': 'GLEWV',
                'start_dis': 0,
                'end': 'DSVKGRFTI',
                'end_dis': 0,
                'limit': 15
                },
            'rf_3':
                {
                'start': 'DTAVYYC',
                'start_dis': 0,
                'end': 'WGQGTLV',
                'end_dis': 0,
                'limit': 20
                },
            'rf_4':
                {
                'start': 'GSGGGGS',
                'start_dis': 24,
                'end': 'LNWYQQK',
                'end_dis': 0,
                'limit': 15
                },
            'rf_5':
                {
                'start': 'APKLLIY',
                'start_dis': 0,
                'end': 'GVPSRFS',
                'end_dis': 0,
                'limit': 15
                },
            'rf_6':
                {
                'start': 'FATYYCQ',
                'start_dis': 0,
                'end': 'FGGGTKV',
                'end_dis': 0,
                'limit': 20
                },
            }
        },
    'nt' : {
        0: {
            'rf_1':
                {
                'start': 'GCCCAGCCGGCCATG',
                'start_dis': 81,
                'end': 'TGGGTCCGCCAGGCT',
                'end_dis': 0,
                'limit': 45
                },
            'rf_2':
                {
                'start': 'GGGCTGGAGTGGGTC',
                'start_dis': 0,
                'end': 'GACTCCGTGAAGGGCCGGTTCACCATC',
                'end_dis': 0,
                'limit': 45
                },
            'rf_3':
                {
                'start': 'GACACGGCCGTGTATTACTGT',
                'start_dis': 0,
                'end': 'TGGGGCCAGGGCACCCTGGTG',
                'end_dis': 0,
                'limit': 75
                },
            'rf_4':
                {
                'start': 'GGATCTGGAGGAGGTGGCTCG',
                'start_dis': 69,
                'end': 'TGGTATCAGCAGAAACCAGGG',
                'end_dis': 0,
                'limit': 45
                },
            'rf_5':
                {
                'start': 'GCCCCTAAGCTCCTGATCTAT',
                'start_dis': 0,
                'end': 'GGGGTCCCATCACGCTTCAGT',
                'end_dis': 0,
                'limit': 45
                },
            'rf_6':
                {
                'start': 'GATTTTGCAACTTACTACTGT',
                'start_dis': 0,
                'end': 'TTCGGCGGAGGTACCAAGGTG',
                'end_dis': 0,
                'limit': 75
                }
            },
        'SH_VH3-23_Vk1-39':
            {
            'rf_1':
                {
                'start': 'GCCCAGCCGGCCATG',
                'start_dis': 78,
                'end': 'TGGGTCCGCCAGGCT',
                'end_dis': 0,
                'limit': 45
                },
            'rf_2':
                {
                'start': 'GGGCTGGAGTGGGTC',
                'start_dis': 0,
                'end': 'GACTCCGTGAAGGGCCGGTTCACCATC',
                'end_dis': 0,
                'limit': 45
                },
            'rf_3':
                {
                'start': 'GACACGGCCGTGTATTACTGT',
                'start_dis': 0,
                'end': 'TGGGGCCAGGGCACCCTGGTG',
                'end_dis': 0,
                'limit': 75
                },
            'rf_4':
                {
                'start': 'GGATCTGGAGGAGGTGGCTCG',
                'start_dis': 69,
                'end': 'TGGTATCAGCAGAAACCAGGG',
                'end_dis': 0,
                'limit': 45
                },
            'rf_5':
                {
                'start': 'GCCCCTAAGCTCCTGATCTAT',
                'start_dis': 0,
                'end': 'GGGGTCCCATCACGCTTCAGT',
                'end_dis': 0,
                'limit': 45
                },
            'rf_6':
                {
                'start': 'GATTTTGCAACTTACTACTGT',
                'start_dis': 0,
                'end': 'TTCGGCGGAGGTACCAAGGTG',
                'end_dis': 0,
                'limit': 75
                }
            },
    'AI_VH3-23_Vk3-20': 
        {
        'rf_1':
            {
            'start': 'GCCCAGCCGGCCATG',
            'start_dis': 81,
            'end': 'TGGGTCCGCCAGGCT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_2':
            {
            'start': 'GGGCTGGAGTGGGTC',
            'start_dis': 0,
            'end': 'GACTCCGTGAAGGGCCGGTTCACCATC',
            'end_dis': 0,
            'limit': 45
            },
        'rf_3':
            {
            'start': 'GACACGGCCGTGTATTACTGT',
            'start_dis': 0,
            'end': 'TGGGGCCAGGGCACCCTGGTG',
            'end_dis': 0,
            'limit': 75
            },
        'rf_4':
            {
            'start': 'GGATCTGGAGGAGGTGGCTCG',
            'start_dis': 72,
            'end': 'TTAAATTGGTATCAGCAGAAA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_5':
            {
            'start': 'GCCCCTAAGCTCCTGATCTAT',
            'start_dis': 0,
            'end': 'GGGGTCCCATCACGCTTCAGT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_6':
            {
            'start': 'TTTGCAACTTACTACTGTCAA',
            'start_dis': 0,
            'end': 'TTCGGCGGAGGTACCAAGGTG',
            'end_dis': 0,
            'limit': 75
            }
        },
    'AI_VH3-23_Vk1-39':
        {
        'rf_1':
            {
            'start': 'GCCCAGCCGGCCATG',
            'start_dis': 81,
            'end': 'TGGGTCCGCCAGGCT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_2':
            {
            'start': 'GGGCTGGAGTGGGTC',
            'start_dis': 0,
            'end': 'GACTCCGTGAAGGGCCGGTTCACCATC',
            'end_dis': 0,
            'limit': 45
            },
        'rf_3':
            {
            'start': 'GACACGGCCGTGTATTACTGT',
            'start_dis': 0,
            'end': 'TGGGGCCAGGGCACCCTGGTG',
            'end_dis': 0,
            'limit': 75
            },
        'rf_4':
            {
            'start': 'GGATCTGGAGGAGGTGGCTCG',
            'start_dis': 72,
            'end': 'TTAAATTGGTATCAGCAGAAA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_5':
            {
            'start': 'GCCCCTAAGCTCCTGATCTAT',
            'start_dis': 0,
            'end': 'GGGGTCCCATCACGCTTCAGT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_6':
            {
            'start': 'TTTGCAACTTACTACTGTCAA',
            'start_dis': 0,
            'end': 'TTCGGCGGAGGTACCAAGGTG',
            'end_dis': 0,
            'limit': 75
            }
        },
    'AI_VH1-69_Vk3-20':
        {
        'rf_1':
            {
            'start': 'GCCCAGCCGGCCATG',
            'start_dis': 81,
            'end': 'TGGGTCCGCCAGGCT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_2':
            {
            'start': 'GGGCTGGAGTGGATG',
            'start_dis': 0,
            'end': 'CAGAAGTTCCAGGGCCGAGTCACCATC',
            'end_dis': 0,
            'limit': 45
            },
        'rf_3':
            {
            'start': 'GACACGGCTGTGTATTACTGT',
            'start_dis': 0,
            'end': 'TGGGGCCAAGGAACCCTGGTC',
            'end_dis': 0,
            'limit': 75
            },
        'rf_4':
            {
            'start': 'GGATCTGGAGGAGGTGGCTCG',
            'start_dis': 72,
            'end': 'TTAGCCTGGTACCAACAGAAA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_5':
            {
            'start': 'GCCCCTAGGCTCCTGATCTAT',
            'start_dis': 0,
            'end': 'GGGATCCCTGACAGGTTCAGT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_6':
            {
            'start': 'TTTGCAGTTTATTACTGTCAA',
            'start_dis': 0,
            'end': 'TTCGGCCAAGGGACCAAGGTG',
            'end_dis': 0,
            'limit': 75
            }
        },
    'AI_VH1-69_Vk1-39':
        {
        'rf_1':
            {
            'start': 'GCCCAGCCGGCCATG',
            'start_dis': 81,
            'end': 'TGGGTCCGCCAGGCT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_2':
            {
            'start': 'GGGCTGGAGTGGATG',
            'start_dis': 0,
            'end': 'CAGAAGTTCCAGGGCCGAGTCACCATC',
            'end_dis': 0,
            'limit': 45
            },
        'rf_3':
            {
            'start': 'GACACGGCTGTGTATTACTGT',
            'start_dis': 0,
            'end': 'TGGGGCCAAGGAACCCTGGTC',
            'end_dis': 0,
            'limit': 75
            },
        'rf_4':
            {
            'start': 'GGATCTGGAGGAGGTGGCTCG',
            'start_dis': 72,
            'end': 'TTAAATTGGTATCAGCAGAAA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_5':
            {
            'start': 'GCCCCTAAGCTCCTGATCTAT',
            'start_dis': 0,
            'end': 'GGGGTCCCATCACGCTTCAGT',
            'end_dis': 0,
            'limit': 45
            },
        'rf_6':
            {
            'start': 'TTTGCAACTTACTACTGTCAA',
            'start_dis': 0,
            'end': 'TTCGGCGGAGGTACCAAGGTG',
            'end_dis': 0,
            'limit': 75
            }
        },
    'PureLibra':
        {
        'rf_1':
            {
            'start': 'GCCCAGCCGGCCATG',
            'start_dis': 78,
            'end': 'GCCATGAGCTGGGTC',
            'end_dis': 0,
            'limit': 45
            },
        'rf_2':
            {
            'start': 'GGGCTGGAGTGGGTGTCAGCT',
            'start_dis': 0,
            'end': 'TACTACGCAGACTCCGTGAAG',
            'end_dis': 0,
            'limit': 45
            },
        'rf_3':
            {
            'start': 'ACGGCCGTATATTACTGTGCG',
            'start_dis': 0,
            'end': 'GACCCCTGGGGCCAGGGAACCCTG',
            'end_dis': 0,
            'limit': 75
            },
        'rf_4':
            {
            'start': 'GGCAGCGGCGGTGGCGGGTCG',
            'start_dis': 81,
            'end': 'CTTAGCCTGGTACCAGCAGAA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_5':
            {
            'start': 'GCTCCCAGGCTCCTCATCTAT',
            'start_dis': 0,
            'end': 'AGCAGGGCCACTGGCATCCCA',
            'end_dis': 0,
            'limit': 45
            },
        'rf_6':
            {
            'start': 'GATTTTGCAGTGTATTACTGT',
            'start_dis': 0,
            'end': 'TTTGGCCAGGGGACCAAGCTG',
            'end_dis': 0,
            'limit': 75
            }
        }    
    }
}

with open(mypath.joinpath(filename), 'w') as file:
    json.dump(rf_dict, file, indent=4)