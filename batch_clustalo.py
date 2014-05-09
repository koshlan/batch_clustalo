# batch_clustalo.py
# (c) Koshlan Mayer-Blackwell
# April 25 - May 05, 2014
# Typical invocation: python blast_clustalo.py

### Module Block ###
def get_list_of_files_ending_in(path, fileextension):
    '''Usage: Return a list of all files in <path_directory> ending with <file_extension>.
    Input: 
        path                str     directory to search
        file_extension      str     file extention to search
    Output:
         list_of_filenames list    a list of filenames'''
    import os
    list_of_filenames = os.listdir(path)
    list_of_filenames = [x for x in list_of_filenames if x.endswith(fileextension)]
    return list_of_filenames

def compile_into_single_file(path, filelist, outputfilename):
    '''Usage: Given a str<path> to a folder and a list<filelist> append all files in order of the list to a single str<outputfile>'''
    oh = open(outputfilename, 'w')
    for filename in filelist:
        fh = open(path + filename ,'r')
        for line in fh:
            oh.write(line)
        fh.close()
    oh.close()

def clean_any_blank_line(file, outputfilename):
    '''Usage: Given a str<file>, return a str<outputfilename> file with all blank lines removed'''
    # ^ is the beginning of string anchor
    #     $ is the end of string anchor
    #     \s is the whitespace character class
    #     * is zero-or-more repetition of
    # ^\s*$
    import re
    fh = open(file, 'r')
    oh = open(outputfilename, 'w')    
    for line in fh:
        if re.match("^\s*$", line):
            pass
        else: 
            oh.write(line)
    fh.close()
    oh.close()

def return_seq_dict_from_fasta(fasta_file):
    '''Usage: Return A Seq Dictionary, keyed on record.id'''
    from Bio import SeqIO
    handle = open(fasta_file, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict

def subfasta_from_dict(list, record_dict, outputfile):
    '''Usage: Write a subfasta to <outputfile> or records from a <list> of record.id looked up in a <record_dict> from a Seq Dictionary'''
    import sys
    oh = open(outputfile, 'w')
    for x in list:
        try:
            oh.write(">%s\n%s\n" %(record_dict[x].id, record_dict[x].seq))
            
        except KeyError:
            sys.stderr.write("%s not found in Records"%(x))
            pass
    oh.close()

def sort_blast(fn_blast_output_name, fn_blast_best_hit):
     import sys
     import os
     os.system('sort -k1,1 -k12,12gr -k11,11g  %s | sort -u -k1,1 --merge > %s' %(fn_blast_output_name, fn_blast_best_hit))

def custom_currate_blastp_output(fn_blast_best_hit):
    fh = open(fn_blast_best_hit, 'r')
    D= {}
    Dc = {}
    Dg = {}
    for line in fh:
        try:
            a = line.split("\t")[1] # HIT IN REFERENCE
            b = line.split("\t")[0] # QUERY GENE
            g = b.split("_")[0]  # GENOME OF QUERRY GENE
            l = line.split("\t")[3] # ALLIGNMENT LENGTH
            ev = line.split("\t")[-2] # EVAL
            bs = line.strip().split("\t")[-1] # BIT SCORE
        except IndexError or ValueError:
            continue
        if int(l) > 100: # CONSIDERS ONLY > 100 AMINO ACID LENGTH PROTEINS
            if a not in D.keys():
              D[a] = {}
            if b not in D[a].keys():
              D[a][b] = True, l,ev, bs
            if a not in Dc.keys():
                Dc[a]= {'count' : 0}
            if g not in Dc[a].keys():
                Dc[a][g] = {'bs': bs , 'query_gene': b, 'alignment_length': l}
                Dc[a]['count'] += 1  
            if bs > Dc[a][g]['bs']: # IF WITHIN THE SAME GENOME A HIGHER BIT SCORE HIT EXISTS UPDATE THE DICTIONARY WITH THAT GENE NAME
                Dc[a][g]={'bs': bs , 'query_gene': b, 'alignment_length': l}           
    fh.close()
    return Dc # THIS DICTIONARY

def clean_list(list, prefixes):
        ''' THIS IS A ONE OFF FUNCTION TO ADD "NA" WHERE MISSING HITS ARE'''
        original_list = [x.split("_")[0] for x in list]
        
        final_list = ["NA"]*len(prefixes)
        for x,y in zip(original_list,list):
            try:
                ind = prefixes.index(x)
                final_list[ind] = y
            except ValueError:
                pass
        return final_list

def output_hit_list_as_matrix(output_array, fn_lookup, prefixes):
    import os
    ''' THIS IS A ONE OF FUNCTION TO OUTPUT A MATRIX RESULT'''
    fh = open(fn_lookup, 'r')
    D = {}
    for line in fh:
        x,y = line.strip().split("\t")
        D[x.replace(">","")] = y 
    fh.close()
    L = []
    for a,b in sorted(output_array):
        cl = clean_list(b,prefixes)
        sys.stdout.write( a  +"\t"+ D[a] +"\t"+ str(len(b))+"\t"+ "\t".join(map(str,cl)) + "\n")
    fh.close()
    '''Usage
    Inputs: tuple array
    Usage: 
    '''



class kdict():# NOT USED CURRENTLY
    '''kdict is a new class meant to deal with nested dictionaries'''
    def __init__(self, D):
        self.D = D
        self.key1 = str() 
        self.key2 = str()
        self.key3 = str()
        self.key4 = str()
    
    def key_1_3_lookup(self):
        #'''Usage: Lookup all elments in a nested dictionary using a primary and tertiary Key. 
        #D{key1:{key2:{key3 : A} key2: {key3: B}}}
                #    -->    A, B   '''
        output_list = []
        for free_key in self.D[self.key1].keys():
            try:
                output_list.append(self.D[self.key1][free_key][self.key3])
            except TypeError:
                pass 
        return output_list
    # def reduce_dictionary_to_k1_with_greater_than_n_genomic_matches(self, special_key = "count"):
    #     out.D
    #     for self.D
    #     return out_D
    





import sys
import os 
import subprocess
blast_result_best_hits = "Example_Inputs/NC_002936_195.gbk_converted.faa.blastresult.best_hit" 
number_of_genomes_found_in = 0 #13
D = custom_currate_blastp_output(blast_result_best_hits) # Makes a Dictionary
filtered_dictionary = {k: v for k, v in D.iteritems() if int(v['count']) > number_of_genomes_found_in } # Refines The Dictionary

def rid(x):
    '''Usage: Removes an entry from a inputed <dictionary> and 
    outpust the resulting <dictionary>'''
    x.pop('count')
    return x 

filtered_dictionary = {k: rid(v) for k, v in filtered_dictionary.iteritems()}

def dump(D,my_key):
    '''Usage: from a dictionary output a list of the values associated with all keys and a secondary <my_key>'''
    L = list()
    for i in sorted(D.keys()):
        L.append(D[i][my_key])
    return L

output_array = [(k,dump(v,'query_gene'))for k,v in filtered_dictionary.iteritems()] # Tuple Array

from datetime import datetime
unique_timestamp = datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
os.system("mkdir %s"%(unique_timestamp))


# ONLY UNQUOTE THIS BLOCK IF YOU WANT A MATRIX OF ALL THE HITS (TO MAKE THIS WORK YOU WILL NEED TO INCLUDE THE PROPER PREFIXES SO THE COLUMNS OF THE MATRIX LINE UP)
fn_lookup = '195.gbk.name_lookup'
prefixes = ["AHS","DET","DehaBAV1","DehalGT","Dehly","DhcVS","DscP1","DscP2","GY50","RBG1351","RBG2","btf","cbdb","dcmb"]
output_hit_list_as_matrix(output_array, fn_lookup, prefixes)


fasta_file = "Example_Inputs/PanGenome_DhcRelatedChloroflexi.fna"
seq_dict = return_seq_dict_from_fasta(fasta_file)
for ref_seq_name , hits_list in output_array:
    subfasta_from_dict(hits_list, seq_dict, "temp.fasta")
    subprocess.call('../util/clustalo -i %s -o %s/%s.fna -v' %("temp.fasta", unique_timestamp, ref_seq_name) , shell=True)