import argparse
import numpy as np
import pandas as pd
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModelForCausalLM
import transformers,torch


def chunkstring(string, length):

    '''
    Chunks string into sub-strings of specified max length.
    Returns list of chunks.
    '''

    return (list(string[0+i:length+i] for i in range(0, len(string), length)))

def get_surprisal(input_str,model,tokenizer,ws_ind,char_repl):

    '''
    Returns surprisal for the last word of a string.

    input_str (str): input string
    model (object): pre-loaded model instance
    tokenizer (object): pre-loaded tokenizer instance
    ws_ind (char): special character indicating initial whitespace when using convert_ids_to_tokens
    char_repl (bool) : whether character replacement is needed (unicode problems in some GPT-models)

    '''

    model.eval()
    if hasattr(model.config, "max_position_embeddings"): # attribute that contains context size is dependent on model-specific config
        ctx_window = model.config.max_position_embeddings # Llama2 config
    elif hasattr(model.config, "n_positions"):
        ctx_window = model.config.n_positions # GPT config

    chunk_size = int(0.75*ctx_window) # chunk size based on LLM's context window size

    seq_chunks = chunkstring(input_str.split(),chunk_size) # returns chunks with words as items

    words, surprisals = [] , []

    for seq in seq_chunks:

        subword_tokens, subword_surprisals = [] , []
        
        inputs = tokenizer(seq, is_split_into_words=True)

        model_inputs = transformers.BatchEncoding({"input_ids":torch.tensor(inputs.input_ids).unsqueeze(0),
            "attention_mask":torch.tensor(inputs.attention_mask).unsqueeze(0)})

        with torch.no_grad():
            outputs = model(**model_inputs)
        
        output_ids = model_inputs.input_ids.squeeze(0)[1:]
        tokens = tokenizer.convert_ids_to_tokens(model_inputs.input_ids.squeeze(0))[1:]
        index = torch.arange(0, output_ids.shape[0])
        surp = -1 * torch.log2(F.softmax(outputs.logits, dim = -1).squeeze(0)[index, output_ids])

        subword_tokens.extend(tokens)
        subword_surprisals.extend(np.array(surp))

        # Word surprisal
        i = 0
        temp_token = ""
        temp_surprisal = 0
        
        while i <= len(subword_tokens)-1:

            temp_token += subword_tokens[i]
            temp_surprisal += subword_surprisals[i]
            
            if i == len(subword_tokens)-1 or tokens[i+1].startswith(ws_ind):
                # remove start-of-token indicator
                words.append(temp_token[1:])
                surprisals.append(temp_surprisal)
                # reset temp token/surprisal
                temp_surprisal = 0
                temp_token = ""
            i += 1

    if char_repl:
        replace_dict = {'ÃĦ':'Ä','Ã¤':'ä','Ãĸ':'Ö','Ã¶':'ö','Ãľ':'Ü','Ã¼':'ü',
                        'ÃŁ':'ß','âĢľ':'“','âĢŀ':'„','Ãł':'à','ÃĢ':'À','Ã¡':'á',
                        'Ãģ':'Á','Ã¨':'è','ÃĪ':'È','Ã©':'é','Ãī':'É','Ã»':'û',
                        'ÃĽ':'Û','ÃŃ':'í','âĢĵ':'–','âĢĻ':'’'}
        for k in replace_dict.keys():
            words = [w.replace(k,replace_dict[k]) for w in words]

    return surprisals[-1]


##############
#### data ####
##############

def adsbc21_surprisal():
    df = pd.read_csv('../data/adsbc21/adsbc21.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model,tokenizer,ws_ind,char_repl))
    df.to_csv('../data/adsbc21/adsbc21.csv', sep = ';', index = False)

def dbc19_surprisal():
    df = pd.read_csv('../data/dbc19/dbc19.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model,tokenizer,ws_ind,char_repl))
    df.to_csv('../data/dbc19/dbc19.csv', sep = ';', index = False)

def adbc23_surprisal():
    df = pd.read_csv('../data/adbc23/adbc23.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model,tokenizer,ws_ind,char_repl))
    df.to_csv('../data/adbc23/adbc23.csv', sep = ';', index = False)

###########################################################################################
###########################################################################################

if __name__ == '__main__':

    all_model_ids = ['leo13b','secret-gpt-2']

    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--model',help = f'Models:{all_model_ids}')
    args = parser.parse_args()
    if args.model == 'leo13b':
        model_id = 'LeoLM/leo-hessianai-13b'
        surp_id = 'leo13b_surp'
        ws_ind = "▁" # Unicode code point is U+2581, not U+005F
        char_repl = False
        
    elif args.model == 'secret-gpt-2':
        model_id = 'stefan-it/secret-gpt2'
        surp_id = 'secretgpt2_surp'
        ws_ind = 'Ġ'
        char_repl = True

    print(f'Model id: {model_id}')
    tokenizer = AutoTokenizer.from_pretrained(model_id, add_prefix_space=True)
    model = AutoModelForCausalLM.from_pretrained(model_id)

    #test_sent = "Gestern schärfte der Holzfäller, bevor er das Holz stapelte, die Axt"
    #get_surprisal(test_sent,model,tokenizer,ws_ind,char_repl)
    
    print('Collecting adsbc21 surprisal...')
    adsbc21_surprisal()
    print('Collecting dbc19 surprisal...')
    dbc19_surprisal()
    print('Collecting adbc23 surprisal...')
    adbc23_surprisal()    

