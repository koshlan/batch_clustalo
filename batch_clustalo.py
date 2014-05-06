

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
    


# SCRIPT ONE CREATE A ALL VS. ALL VS. BLAST
# CURATE BLAST FOR

# CONCISE INVOCATION #
import sys
import os 
import subprocess
blast_result_best_hits = "Example_Inputs/NC_002936_195.gbk_converted.faa.blastresult.best_hit" #'/Users/koshlan/Dropbox/Mayer_Kaster_Transfer_Folder/ClusterViz_v1_for_use/Example_Inputs/NC_002936_195.gbk_converted.faa.blastresult.best_hit'
number_of_genomes_found_in = 13

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

fasta_file = "Example_Inputs/PanGenome_DhcRelatedChloroflexi.fna"
seq_dict = return_seq_dict_from_fasta(fasta_file)
for ref_seq_name , hits_list in output_array:
    subfasta_from_dict(hits_list, seq_dict, "temp.fasta")
    subprocess.call('./util/clustalo -i %s -o %s/%s.fna -v' %("temp.fasta", unique_timestamp, ref_seq_name) , shell=True)
    break

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    #filtered_list = filter(lambda x: int(D[x]['count']) > 13, D.keys()) # Gets Keys where D['count'] 
    
    
    # DET0541    Dehly_0636  dcmb_550    DhcVS_482   DscP2_2527150851    DscP1_05    GY50_0467   DET_0541    RBG1351_27_18   AHS_6666666.47752.peg.1339  DehalGT_0480    btf_504 cbdb_A515   DehaBAV1_0517   RBG2_3_123
    #  DET0540  Dehly_0635  dcmb_549    DhcVS_481   DscP2_2527150852    DscP1_06    GY50_0466   DET_0540    RBG1351_27_19   AHS_6666666.47752.peg.1338  DehalGT_0479    btf_503 cbdb_A514   DehaBAV1_0516   RBG2_3_124
    #  DET0543  Dehly_0638  dcmb_552    DhcVS_484   DscP2_2527151626    DscP1_025   GY50_0469   DET_0543    RBG1351_27_15   AHS_6666666.47752.peg.354   DehalGT_0482    btf_506 cbdb_A517   DehaBAV1_0519   RBG2_13_5
    #  DET0542  Dehly_0637  dcmb_551    DhcVS_483   DscP2_2527150850    DscP1_04    GY50_0468   DET_0542    RBG1351_27_17   AHS_6666666.47752.peg.1374  DehalGT_0481    btf_505 cbdb_A516   DehaBAV1_0518   RBG2_3_122
    #  DET0545  Dehly_0640  dcmb_554    DhcVS_486   DscP2_2527151889    DscP1_09    GY50_0471   DET_0545    RBG1351_27_13   AHS_6666666.47752.peg.1344  DehalGT_0484    btf_508 cbdb_A520   DehaBAV1_0521   RBG2_3_121
    #  DET0544  Dehly_0639  dcmb_553    DhcVS_485   DscP2_2527150846    DscP1_03    GY50_0470   DET_0544    RBG1351_27_14   AHS_6666666.47752.peg.803   DehalGT_0483    btf_507 cbdb_A519   DehaBAV1_0520   RBG2_3_56
    #  DET0727  Dehly_0626  dcmb_694    DhcVS_633   DscP2_2527151562    DscP1_0169  GY50_0616   DET_0727    RBG1351_2_8 AHS_6666666.47752.peg.1108  DehalGT_0620    btf_648 cbdb_A682   DehaBAV1_0659   RBG2_17_14
    #  DET0726  Dehly_0627  dcmb_693    DhcVS_632   DscP2_2527150814    DscP1_0168  GY50_0615   DET_0726    RBG1351_2_9 AHS_6666666.47752.peg.494   DehalGT_0619    btf_647 cbdb_A681   DehaBAV1_0658   RBG2_33_3
    #  DET0743  Dehly_0700  dcmb_711    DhcVS_650   DscP2_2527151015    DscP1_0105  GY50_0630   DET_0743    RBG1351_44_11   AHS_6666666.47752.peg.443   DehalGT_0634    btf_665 cbdb_A718   DehaBAV1_0673   RBG2_4_110
    #  DET1268  Dehly_0233  dcmb_1127   DhcVS_1051  DscP2_2527151947    DscP1_0444  GY50_1099   DET_1268    RBG1351_14_12   AHS_6666666.47752.peg.1064  DehalGT_1005    btf_1146    cbdb_A1194  DehaBAV1_1079   RBG2_4_27
    #  DET0527  Dehly_0371  dcmb_536    DhcVS_468   DscP2_2527150188    DscP1_095   GY50_0453   DET_0527    RBG1351_26_8    AHS_6666666.47752.peg.615   DehalGT_0466    btf_490 cbdb_A496   DehaBAV1_0503   RBG2_4_116
    #  DET0526  Dehly_0372  dcmb_535    DhcVS_467   DscP2_2527150189    DscP1_097   GY50_0452   DET_0526    RBG1351_26_7    AHS_6666666.47752.peg.616   DehalGT_0465    btf_489 cbdb_A494   DehaBAV1_0502   RBG2_4_115
    #  DET0525  Dehly_0373  dcmb_534    DhcVS_466   DscP2_2527150191    DscP1_099   GY50_0451   DET_0525    RBG1351_26_6    AHS_6666666.47752.peg.617   DehalGT_0464    btf_488 cbdb_A493   DehaBAV1_0501   RBG2_4_114
    #  DET0564  Dehly_1439  dcmb_571    DhcVS_1189  DscP2_2527150291    DscP1_0111  GY50_0489   DET_0564    RBG1351_17_2    AHS_6666666.47752.peg.242   DehalGT_0501    btf_525 cbdb_A538   DehaBAV1_1215   RBG2_73_1
    #  DET1606  Dehly_1345  dcmb_1492   DhcVS_1488  DscP2_2527151082    DscP1_0211  GY50_1500   DET_1606    RBG1351_24_6    AHS_6666666.47752.peg.1078  DehalGT_1400    btf_1546    cbdb_A1700  DehaBAV1_1352   RBG2_3_25
    #  DET0439  Dehly_0150  dcmb_451    DhcVS_381   DscP2_2527150942    DscP1_0104  GY50_0364   DET_0439    RBG1351_11_35   AHS_6666666.47752.peg.624   DehalGT_0380    btf_405 cbdb_A393   DehaBAV1_0416   RBG2_4_111
    #  DET0976  Dehly_0529  dcmb_911    DhcVS_849   DscP2_2527151353    DscP1_0131  GY50_0858   DET_0976    RBG1351_27_1    AHS_6666666.47752.peg.150   DehalGT_0820    btf_925 cbdb_A939   DehaBAV1_0867   RBG2_2_33
    #  DET1278  Dehly_0304  dcmb_1137   DhcVS_1061  DscP2_2527151133    DscP1_015   GY50_1110   DET_1278    RBG1351_6_15    AHS_6666666.47752.peg.1347  DehalGT_1015    btf_1156    cbdb_A1206  DehaBAV1_1089   RBG2_3_117
    #  DET0009  Dehly_1382  dcmb_9  DhcVS_9 DscP2_2527151246    DscP1_086   GY50_0010   DET_0009    RBG1351_5_7 AHS_6666666.47752.peg.126   DehalGT_0009    btf_9   cbdb_A11    DehaBAV1_0009   RBG2_4_140
    #  DET1607  Dehly_1344  dcmb_1493   DhcVS_1489  DscP2_2527150766    DscP1_0212  GY50_1501   DET_1607    RBG1351_24_7    AHS_6666666.47752.peg.1077  DehalGT_1401    btf_1547    cbdb_A1701  DehaBAV1_1353   RBG2_3_26
    #  DET1372  Dehly_0067  dcmb_1231   DhcVS_1153  DscP2_2527150145    DscP1_0480  GY50_1211   DET_1372    RBG1351_37_9    AHS_6666666.47752.peg.752   DehalGT_1090    btf_1250    cbdb_A1324  DehaBAV1_1183   RBG2_7_74
    #  DET1609  Dehly_0010  dcmb_1495   DhcVS_1491  DscP2_2527151079    DscP1_0318  GY50_1503   DET_1609    RBG1351_24_9    AHS_6666666.47752.peg.679   DehalGT_1403    btf_1549    cbdb_A1703  DehaBAV1_1355   RBG2_3_8
    #  DET1608  Dehly_1343  dcmb_1494   DhcVS_1490  DscP2_2527151080    DscP1_046   GY50_1502   DET_1608    RBG1351_24_8    AHS_6666666.47752.peg.678   DehalGT_1402    btf_1548    cbdb_A1702  DehaBAV1_1354   RBG2_1_92
    #  DET1276  Dehly_0302  dcmb_1135   DhcVS_1059  DscP2_2527151845    DscP1_014   GY50_1108   DET_1276    RBG1351_42_2    AHS_6666666.47752.peg.550   DehalGT_1013    btf_1154    cbdb_A1204  DehaBAV1_1087   RBG2_3_118
    #  DET1211  Dehly_0967  dcmb_1074   DhcVS_994   DscP2_2527151358    DscP1_0321  GY50_1016   DET_1211    RBG1351_2_81    AHS_6666666.47752.peg.1042  DehalGT_0952    btf_1091    cbdb_A1128  DehaBAV1_1021   RBG2_3_10
    #  DET0963  Dehly_0979  dcmb_898    DhcVS_836   DscP2_2527151350    DscP1_0141  GY50_0845   DET_0963    RBG1351_34_11   AHS_6666666.47752.peg.598   DehalGT_0807    btf_912 cbdb_A921   DehaBAV1_0854   RBG2_6_113
    #  DET0539  Dehly_0634  dcmb_548    DhcVS_480   DscP2_2527150853    DscP1_07    GY50_0465   DET_0539    RBG1351_27_20   AHS_6666666.47752.peg.1337  DehalGT_0478    btf_502 cbdb_A513   DehaBAV1_0515   RBG2_3_125
    #  DET1210  Dehly_0966  dcmb_1073   DhcVS_993   DscP2_2527151357    DscP1_043   GY50_1015   DET_1210    RBG1351_2_82    AHS_6666666.47752.peg.1041  DehalGT_0951    btf_1090    cbdb_A1127  DehaBAV1_1020   RBG2_3_11
    #  DET1502  Dehly_0066  dcmb_1356   DhcVS_1276  DscP2_2527151125    DscP1_0166  GY50_1323   DET_1502    RBG1351_5_24    AHS_6666666.47752.peg.1250  DehalGT_1205    btf_1370    cbdb_A1474  DehaBAV1_1292   RBG2_4_12
    #  DET0535  Dehly_0693  dcmb_544    DhcVS_476   DscP2_2527150857    DscP1_0435  GY50_0461   DET_0535    RBG1351_40_2    AHS_6666666.47752.peg.669   DehalGT_0474    btf_498 cbdb_A509   DehaBAV1_0511   RBG2_11_31
    #  DET0536  Dehly_0694  dcmb_545    DhcVS_477   DscP2_2527150856    DscP1_0436  GY50_0462   DET_0536    RBG1351_40_1    AHS_6666666.47752.peg.670   DehalGT_0475    btf_499 cbdb_A510   DehaBAV1_0512   RBG2_11_32
    #  DET1630  Dehly_1374  dcmb_1516   DhcVS_1512  DscP2_2527151571    DscP1_024   GY50_1153   DET_1630    RBG1351_23_13   AHS_6666666.47752.peg.1327  DehalGT_1423    btf_1570    cbdb_A1728  DehaBAV1_1375   RBG2_7_21
    #  DET0437  Dehly_0211  dcmb_449    DhcVS_379   DscP2_2527151459    DscP1_0100  GY50_0362   DET_0437    RBG1351_11_39   AHS_6666666.47752.peg.622   DehalGT_0378    btf_403 cbdb_A391   DehaBAV1_0414   RBG2_78_2
    #  DET0004  Dehly_1362  dcmb_4  DhcVS_4 DscP2_2527151611    DscP1_081   GY50_0005   DET_0004    RBG1351_16_2    AHS_6666666.47752.peg.700   DehalGT_0004    btf_4   cbdb_A4 DehaBAV1_0004   RBG2_4_145
    #  DET0038  Dehly_1208  dcmb_37 DhcVS_36    DscP2_2527151981    DscP1_0222  GY50_0039   DET_0038    RBG1351_4_5 AHS_6666666.47752.peg.143   DehalGT_0037    btf_36  cbdb_A46    DehaBAV1_0035   RBG2_158_1
    #  DET1035  Dehly_0648  dcmb_971    DhcVS_905   DscP2_2527150235    DscP1_0492  GY50_0919   DET_1035    RBG1351_13_38   AHS_6666666.47752.peg.1176  DehalGT_0871    btf_988 cbdb_A1008  DehaBAV1_0917   RBG2_3_52
    #  DET1202  Dehly_0915  dcmb_1065   DhcVS_985   DscP2_2527151325    DscP1_0336  GY50_1006   DET_1202    RBG1351_22_15   AHS_6666666.47752.peg.461   DehalGT_0943    btf_1082    cbdb_A1118  DehaBAV1_1012   RBG2_6_69
    #  DET1204  Dehly_0964  dcmb_1067   DhcVS_987   DscP2_2527150759    DscP1_0245  GY50_1008   DET_1204    RBG1351_16_1    AHS_6666666.47752.peg.576   DehalGT_0945    btf_1084    cbdb_A1120  DehaBAV1_1014   RBG2_2_122
    #  DET1342  Dehly_0142  dcmb_1201   DhcVS_1124  DscP2_2527151941    DscP1_0163  GY50_1181   DET_1342    RBG1351_7_29    AHS_6666666.47752.peg.1328  DehalGT_1078    btf_1220    cbdb_A1292  DehaBAV1_1153   RBG2_1_295
    #  DET1263  Dehly_0935  dcmb_1122   DhcVS_1046  DscP2_2527151351    DscP1_0142  GY50_1094   DET_1263    RBG1351_14_15   AHS_6666666.47752.peg.1408  DehalGT_1000    btf_1141    cbdb_A1187  DehaBAV1_1074   RBG2_4_26
    #  DET0538  Dehly_0633  dcmb_327    DhcVS_479   DscP2_2527150854    DscP1_08    GY50_0464   DET_0538    RBG1351_27_21   AHS_6666666.47752.peg.1336  DehalGT_0477    btf_501 cbdb_A512   DehaBAV1_0514   RBG2_3_126






#for k in sorted(filtered_dictionary.keys()):
#    print my_dump(D[k],'gene_query')
    
    
# 
# print sorted(filtered_dictionary.keys())
# print
# print sorted(filtered_list)
# points={'a':{'counts': 44}, 'b':{'counts': 3}, 'c':{'counts': 24}, 'd':{'counts': 1}}
# XX = {k: v for k, v in points.iteritems() if int(v['counts']) > 13}




# import pprint 
# pp = pprint.PrettyPrinter(indent = 4)
# pp.pprint(D)
# K = kdict(D)
# K.key1 = "DET1642"
# K.key3 = "query_gene"
# print K.key_1_3_lookup()
# 
# 
# 
# 
# print K.D.keys()
# 
# print len(u)
# print len(D.keys())





# Currate a Blastp 
# directory = '/Users/koshlan/Dropbox/Mayer_Kaster_Transfer_Folder/kosh_gbk_files/Converted/'
# my_list = get_list_of_files_ending_in(directory, '.fna')   
# compile_into_single_file(directory, my_list, "COMPILED.fna")
# clean_any_blank_line("COMPILED.fna", "COMPILED_ClEAN.fna" )
# D = return_seq_dict_from_fasta("COMPILED_ClEAN.fna")
# my_list = ["dcmb_998", "AHS_6666666.47752.peg.687"]
# for token in ["1","2","3"]:
#     if not os.path.exists("allignments"):
#         os.mkdir("allignments")
#     if not os.path.exists("allignments/%s/"%(token)):
#         os.mkdir("allignments/%s/"%(token))
#     subfasta(my_list,D,"./allignments/%s/subfasta.test.txt"%(token))




# ## #?# ##
# def key_1_3_lookup(dictionary,key1,key3):
#     output_list = []
#     for free_key in dictionary[key1].keys():
#         try:
#             output_list.append(dictionary[key1][free_key][key3])
#         except TypeError:
#             pass 
#     return output_list
# 
# 
