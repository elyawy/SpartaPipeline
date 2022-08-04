 
length_distributions = ["zipf", "geometric", "poisson"]
indel_models = ["sim","rim"]
PARAMS_LIST = [
    "root_length",
    "length_param_insertion",
    "length_param_deletion",
    "insertion_rate",
    "deletion_rate"
]
SUMSTATS_LIST = [f'SS_{i}' for i in range(0,27)]
SUMSTATS_DEFINITION = {
    'SS_0': "AVG_GAP_SIZE",
    'SS_1': "MSA_LEN",
    'SS_2': "MSA_MAX_LEN",
    'SS_3': "MSA_MIN_LEN",
    'SS_4': "TOT_NUM_GAPS",
    'SS_5': "NUM_GAPS_LEN_ONE",
    'SS_6': "NUM_GAPS_LEN_TWO",
    'SS_7': "NUM_GAPS_LEN_THREE",
    'SS_8': "NUM_GAPS_LEN_AT_LEAST_FOUR",
    'SS_9': "AVG_UNIQUE_GAP_SIZE",
    'SS_10': "TOT_NUM_UNIQUE_GAPS",
    'SS_11': "NUM_GAPS_LEN_ONE\nPOS_1_GAPS",
    'SS_12': "NUM_GAPS_LEN_ONE\nPOS_2_GAPS",
    'SS_13': "NUM_GAPS_LEN_ONE\nPOS_N_MINUS_1_GAPS",
    'SS_14': "NUM_GAPS_LEN_TWO\nPOS_1_GAPS",
    'SS_15': "NUM_GAPS_LEN_TWO\nPOS_2_GAPS",
    'SS_16': "NUM_GAPS_LEN_TWO\nPOS_N_MINUS_1_GAPS",
    'SS_17': "NUM_GAPS_LEN_THREE\nPOS_1_GAPS",
    'SS_18': "NUM_GAPS_LEN_THREE\nPOS_2_GAPS",
    'SS_19': "NUM_GAPS_LEN_THREE\nPOS_N_MINUS_1_GAPS",
    'SS_20': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_1_GAPS",
    'SS_21': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_2_GAPS",
    'SS_22': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_N_MINUS_1_GAPS",
    'SS_23': "MSA_POSITION_WITH_0_GAPS",
    'SS_24': "MSA_POSITION_WITH_1_GAPS",
    'SS_25': "MSA_POSITION_WITH_2_GAPS",
    'SS_26': "MSA_POSITION_WITH_N_MINUS_1_GAPS"
}

