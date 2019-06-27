import sys
import pdb

def read_hmm(afile):
    afile = open(afile, 'r')
    hmm_list = []
    null_model = []
    transition_freq = []
    local_div = []
    in_hmm = False

    for line in afile:
        line = line.strip()
        if line.startswith('NULL'):
            null_model = line.split()[1:]
            

        if line.startswith('HMM'):
            in_hmm = True
            continue

        if in_hmm:
            line_arr = line.split()
            if len(line_arr) == 10: #Not amino acid frequencies --> transition freq
                transition_freq.append(line_arr[0:7])
		local_div.append(line_arr[7:])
            
	    if len(line_arr) == 23:
	        aa = line_arr[0]
                freq = line_arr[2:22]
                hmm_list.append((aa, freq))
    
    afile.close()
    return hmm_list, null_model, transition_freq, local_div


if __name__ == "__main__":

    afile = sys.argv[1]
    hmm_list, null_model, transition_freq, local_div = read_hmm(afile) #The two first entries of the transition_freq and local_div will belong to the null model
    afile.close()
    pdb.set_trace()
    #print hmm_list
