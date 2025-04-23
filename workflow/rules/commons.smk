import pandas as pd

samples = pd.read_table(config['samples'], sep="\t").set_index("sample_id", drop=False)
units = pd.read_table(config['units'], sep="\t").set_index("sample_id", drop=False)
all_files = list(units['path'])
# If there is a single input
input_files = dict(zip(units['sample_id'], units['path']))
condition = dict(zip(samples['sample_id'], samples['condition']))

# Add a new list of sample_id names to the config dictionary to inlude in Snakefile
sample_reads = pd.read_table(config['reads'], sep="\t").set_index("sample_id", drop=False)
reads = pd.read_table(config['reads'], sep="\t").set_index("sample_id", drop=False)

read_1 = dict(zip(reads['sample_id'], reads['read_1']))
read_2 = dict(zip(reads['sample_id'], reads['read_2']))

def run_sample(sample_id):
    '''
    Only process samples based on a specified condition
    '''
    if condition[sample_id] == "treated":
        return(input_files[sample_id])
    else:
        return("Sample is control!")

def get_sample_param(sample_id):
    '''
    Return a sample specific parameter
    '''
    if condition[sample_id] == "treated":
        output = "treated"
    else:
        output = "control"
    
    return(output)
