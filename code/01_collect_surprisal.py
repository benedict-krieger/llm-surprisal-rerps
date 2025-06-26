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


def get_surprisal(input_str, model, tokenizer, ws_ind, char_repl, bos_pad):

    '''
    Returns surprisal for the last word of a string.

    input_str (str): input string
    model (object): pre-loaded model instance
    tokenizer (object): pre-loaded tokenizer instance
    ws_ind (char): special character indicating initial whitespace when using convert_ids_to_tokens
    char_repl (bool) : whether character replacement is needed (unicode problems in some GPT-models)
    bos_pad (bool) : whether input sequence needs to be padded with the bos token

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

        if bos_pad:
            bos_id = tokenizer.bos_token_id
            model_inputs = transformers.BatchEncoding({"input_ids":torch.tensor([bos_id]+inputs.input_ids).unsqueeze(0),
                "attention_mask":torch.tensor([1]+inputs.attention_mask).unsqueeze(0)})
        else:
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


def bpe_split(word, bos_pad):

    '''
    Test if a given (target) word is split by the tokenizer into multiple subwords.
    
    If the tested tokenizer automatically prepends the BOS token (bos_pad=False), length 2 indicates 
    BOS + single token, while length > 2 indicates multiple subwords.
    
    If the tokenizer does not prepend the BOS token (bos_pad=True), length > 1 indicates multiple subwords.
    '''

    encoded_w = tokenizer.encode(word)

    if len(encoded_w) > 2:
        bpe_split = 1
    elif bos_pad and (len(encoded_w) > 1):
        bpe_split = 1
    else:
        bpe_split = 0

    return bpe_split


##############
#### data ####
##############

def adsbc21_surprisal():
    df = pd.read_csv('../data/adsbc21/adsbc21.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model, tokenizer, ws_ind, char_repl, bos_pad))
    df[bpe_id] = df['Target'].apply(bpe_split, args=(bos_pad,))
    df.to_csv('../data/adsbc21/adsbc21.csv', sep = ';', index = False)

def dbc19_surprisal():
    df = pd.read_csv('../data/dbc19/dbc19.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model, tokenizer, ws_ind, char_repl, bos_pad))
    df[bpe_id] = df['Target'].apply(bpe_split, args=(bos_pad,))
    df.to_csv('../data/dbc19/dbc19.csv', sep = ';', index = False)

def adbc23_surprisal():
    df = pd.read_csv('../data/adbc23/adbc23.csv', sep = ';')
    df[surp_id] = df['Stimulus_tf'].apply(get_surprisal, args=(model,tokenizer,ws_ind,char_repl, bos_pad))
    df[bpe_id] = df['Target'].apply(bpe_split, args=(bos_pad,))
    df.to_csv('../data/adbc23/adbc23.csv', sep = ';', index = False)

###########################################################################################
###########################################################################################

if __name__ == '__main__':

    all_models = ['leo13b','gerpt2','gerpt2-large']

    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--model',help = f'Models:{all_models}')
    args = parser.parse_args()
    if args.model == 'leo13b':
        model_id = 'LeoLM/leo-hessianai-13b'
        surp_id = 'leo13b_surp'
        bpe_id = 'leo13b_bpe_split'
        ws_ind = "▁" # Unicode code point is U+2581, not U+005F
        char_repl = False
        bos_pad = False

    elif args.model == 'gerpt2':
        model_id = 'benjamin/gerpt2'
        surp_id = 'gerpt2_surp'
        bpe_id = 'gerpt2_bpe_split'
        ws_ind = 'Ġ'
        char_repl = True
        bos_pad = True

    elif args.model == 'gerpt2-large':
        model_id = 'benjamin/gerpt2-large'
        surp_id = 'gerpt2large_surp'
        bpe_id = 'gerpt2large_bpe_split'
        ws_ind = 'Ġ'
        char_repl = True
        bos_pad = True

    print(f'Model id: {model_id}')
    tokenizer = AutoTokenizer.from_pretrained(model_id, add_prefix_space=True)
    model = AutoModelForCausalLM.from_pretrained(model_id)

    print('Collecting adsbc21 surprisal...')
    adsbc21_surprisal()
    print('Collecting dbc19 surprisal...')
    dbc19_surprisal()
    print('Collecting adbc23 surprisal...')
    adbc23_surprisal()    

