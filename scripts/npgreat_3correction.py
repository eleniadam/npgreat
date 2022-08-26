########################
########################
# NPGREAT - CORRECTION #
########################
########################
# The local alignments between each REXTAL contig with the Nanopore reads are investigated to identify problematic regions within the REXTAL sequence contig.
# Misassemblies caused by deletions are detected and resolved.
# General correction: Seperation of the contig at the point where the deletion is detected in order to use the Nanopore sequence to fill the missing portion.
# Tandem Repeat (TR) region correction: Removal of the stretched TR-region in the REXTAL contig and replaced by the properly represented TR-region in the Nanopore sequence.

###########################
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import csv

folder_loc_orientation = "orientation_step"
folder_loc_position = "position_step"
folder_loc_correction = 'correction_step'

print("- The NPGREAT Correction step begins...")


# Read input
if (len(sys.argv) > 1):
    eval_choice = sys.argv[1]
else:
    eval_choice = 0

# Flag 0: for evalue!=0.0 (16p, 19q, 20p, 22q)
# Flag 1: for evalue=0.0 (all others)
if(eval_choice in ['0', '0.0']):
    flag_eval = 1
else:
    flag_eval = 0
#print("The flag is:", flag_eval)

# Read nano_ids
nano_ids = []
with open(str(folder_loc_position)+'/nano_ids.csv', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    for elem in (list(reader)[1:]):
        nano_ids.append(''.join(elem))
#print(nano_ids)

# Read contig_aligns
contig_aligns = []
l = 0
for elem in nano_ids:
    contig_aligns.append(pd.read_csv(str(folder_loc_position)+'/caligns_'+str(elem)+'.csv', usecols= ['sseqid','pident','length','qstart','qend','sstart','send','qlen','slen']))
    l = l + 1
    
# Read contig_order
contig_order = []
with open(str(folder_loc_position)+'/contig_order1.csv', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    
    for elem in (list(reader)[1:]):
        contig_order.append((((''.join(elem)[1:-1]).replace("'","")).replace(" ","")).split(","))


####################################################
# Detect split by checking the internal alignments #
####################################################
print("Detecting possible splits by checking the internal alignments...")
df_splits = pd.DataFrame()

# Check each nano
for k in range(0, len(nano_ids), 1):
    
    # Nano ID
    nid = nano_ids[k]
    
    # Alignments
    contig_aligns[k] = contig_aligns[k].sort_values(['sseqid', 'qstart'], ascending=False).reset_index(drop=True)
        
    # Nano length
    nlen = contig_aligns[k].loc[0, 'qlen']
    
    # Check every two alignments
    for i in range(0, len(contig_aligns[k]), 1):
        
        # Last alignment
        if(i == len(contig_aligns[k])-1):
            break
        
        a_rextal_id = contig_aligns[k].loc[i, 'sseqid']
        b_rextal_id = contig_aligns[k].loc[i+1, 'sseqid']
        
        clen = contig_aligns[k].loc[i, 'slen']
        a_nano_coord1 = contig_aligns[k].loc[i, 'qstart']
        b_nano_coord2 = contig_aligns[k].loc[i+1, 'qend']
        a_rextal_coord1 = contig_aligns[k].loc[i, 'sstart']
        b_rextal_coord2 = contig_aligns[k].loc[i+1, 'send']
        
        # Contig ID must be the same
        if(a_rextal_id != b_rextal_id):
            continue
        
        # Detect splits
        nano_gap = a_nano_coord1 - b_nano_coord2
        rextal_gap = a_rextal_coord1 - b_rextal_coord2
        sgap = nano_gap - rextal_gap
        
        if(sgap >= 100):
            contigl = b_rextal_coord2 - 1000
            if(contigl <= 0):
                contigl = 1
            contigr = a_rextal_coord1 + 1000
            if(contigr > clen):
                contigr = clen
            nanol = b_nano_coord2 - 1000
            if(nanol <= 0):
                nanol = 1
            nanor = a_nano_coord1 + 1000
            if(nanor > nlen):
                nanor = nlen
            
            # Store possible splits
            # Columns: contig id, nano id, contig coords +/-1kb (or min/max len), nano coords +/-1kb (or min/max len), size_gap
            df_splits = df_splits.append({'contig_id': a_rextal_id, 'contig_lcoord': contigl, 'contig_rcoord': contigr, 'nano_id': nid, 'nano_lcoord': nanol, 'nano_rcoord': nanor, 'size_gap': sgap}, ignore_index=True)

# Sort detected splits by contig ID
df_splits = df_splits.sort_values(['contig_id', 'contig_lcoord'], ascending=True).reset_index(drop=True)
df_splits = df_splits.astype({"contig_lcoord": int, "contig_rcoord": int, "nano_lcoord": int, "nano_rcoord": int})

# Remove splits that contig_lcoord > contig_rcoord
df_splits = df_splits.drop(df_splits[(df_splits["contig_lcoord"] >= df_splits["contig_rcoord"])].index).reset_index(drop=True)


#####################################################
# Run Zoom alignments (blastn) & Investigate splits #
#####################################################
print("Investigating detected possible splits...")

cols_outfmt6 = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.strip().split(' ')

df_splits_locs = pd.DataFrame()

for k in range(0, len(df_splits.index), 1):
    
    # Blastn the zoomed area
    blast_outfmt6 = subprocess.check_output(f"blastn -outfmt 6 -perc_identity 80 -subject {folder_loc_orientation}/{df_splits['contig_id'][k]}.fa -query {folder_loc_orientation}/{df_splits['nano_id'][k]}.fa -subject_loc {df_splits['contig_lcoord'][k]}-{df_splits['contig_rcoord'][k]} -query_loc {df_splits['nano_lcoord'][k]}-{df_splits['nano_rcoord'][k]}", shell=True, encoding='UTF-8')
    try:
        df_zoom_splits = pd.read_table(StringIO(blast_outfmt6), header=None)  
    except:
        # No alignments
        continue
    
    df_zoom_splits.columns = cols_outfmt6
    
    # Clean alignments: Remove reverse alignments
    df_zoom_splits = df_zoom_splits.loc[(df_zoom_splits['qstart'] <= df_zoom_splits['qend']) & (df_zoom_splits['sstart'] <= df_zoom_splits['send'])]
    
    # Select alignments: Keep only alignments with evalue=0
    df_zoom_splits = df_zoom_splits.loc[(df_zoom_splits['evalue'] == 0.0)]
    df_zoom_splits = df_zoom_splits.reset_index(drop=True)
    
    # Sort by nano coordinates
    df_zoom_splits = df_zoom_splits.sort_values(by=['qstart','sstart'], ascending=True)
    df_zoom_splits = df_zoom_splits.reset_index(drop=True)
    
    # Number of alignments
    numAligns = df_zoom_splits.shape[0]
        
    # Investigate alignments
    if (numAligns == 1):
        # If the alignment contains the entire area in question (i.e. entire area aligned)
        # No split
        continue
        
    if (numAligns == 2 or numAligns > 3):
        contig_alA_end = df_zoom_splits.loc[0, 'send']
        contig_alB_start = df_zoom_splits.loc[1, 'sstart']
        nano_alA_end = df_zoom_splits.loc[0, 'qend']
        nano_alB_start = df_zoom_splits.loc[1, 'qstart']
        
        # (i) - TR split
        if ( (contig_alA_end > contig_alB_start) and (nano_alA_end <= nano_alB_start) ):
            #print(str(contig_alB_start) + " - " + str(contig_alA_end))
            df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alB_start, 'region_coord2': contig_alA_end}, ignore_index=True)
        
        # (ii) - TR split
        if( (contig_alA_end > contig_alB_start) and (nano_alA_end > nano_alB_start) ):
            # Nano rep region
            coord1 = nano_alB_start
            coord2 = nano_alA_end
            
            # Length of contig repeated region
            len_rep = contig_alA_end - contig_alB_start
            
            # Extended nano coordinates
            coord1_ext = coord2 - len_rep
            coord2_ext = coord1 + len_rep
            
            # Nano regions where the contig repeat exists
            if( not(coord1_ext > coord1 and coord2 < coord2_ext) and not(coord1 > coord1_ext and coord2_ext < coord2) ):
                #print(str(contig_alB_start) + " - " + str(contig_alA_end))
                df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alB_start, 'region_coord2': contig_alA_end}, ignore_index=True)
        
        if( (contig_alA_end <= contig_alB_start) and (nano_alA_end < nano_alB_start) ):
            len_close = contig_alB_start - contig_alA_end
            len_nano_close = nano_alB_start - nano_alA_end
            ethresh = 5
            
            if( len_close < len_nano_close ):
                
                # (iii) - General split
                if((len_close < ethresh)):
                    #print(contig_alA_end)
                    df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alA_end, 'region_coord2': contig_alA_end}, ignore_index=True)
                
                # (iv) - TR split
                else:
                    #print(str(contig_alA_end) + " - " + str(contig_alB_start))
                    df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alA_end, 'region_coord2': contig_alB_start}, ignore_index=True)


    if (numAligns == 3):
        contig_alA_end = df_zoom_splits.loc[0, 'send']
        contig_alB_start = df_zoom_splits.loc[1, 'sstart']
        nano_alA_end = df_zoom_splits.loc[0, 'qend']
        nano_alB_start = df_zoom_splits.loc[1, 'qstart']
        
        contig_alB_end = df_zoom_splits.loc[1, 'send']
        contig_alC_start = df_zoom_splits.loc[2, 'sstart']
        nano_alB_end = df_zoom_splits.loc[1, 'qend']
        nano_alC_start = df_zoom_splits.loc[2, 'qstart']
        
        # (v) - TR split
        if ( (contig_alA_end > contig_alB_start) and (nano_alA_end <= nano_alB_start) and (contig_alB_end > contig_alC_start) and (nano_alB_end <= nano_alC_start)):
            #print(str(contig_alB_start) + " - " + str(contig_alB_end))
            df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alB_start, 'region_coord2': contig_alB_end}, ignore_index=True)
    
    ####
        # Same as case above
        elif( (contig_alA_end <= contig_alB_start) and (nano_alA_end < nano_alB_start) ):
            len_close = contig_alB_start - contig_alA_end
            len_nano_close = nano_alB_start - nano_alA_end
            ethresh = 5
            
            if( len_close < len_nano_close ):
                
                # (iii) - General split
                if((len_close < ethresh)):
                    #print(contig_alA_end)
                    df_splits_locs = df_splits_locs.append({'contig_id': df_splits['contig_id'][k], 'region_coord1': contig_alA_end, 'region_coord2': contig_alA_end}, ignore_index=True)
                
    ####

df_splits_locs = df_splits_locs.sort_values(['contig_id', 'region_coord1'], ascending=True).reset_index(drop=True)
df_splits_locs = df_splits_locs.astype({"region_coord1": int, "region_coord2": int})

#print(df_splits_locs)

#########################################
# Identify exact splits for each contig #
#########################################
print("Identifying splits...")

# Region length
df_splits_locs['region_length'] = df_splits_locs['region_coord2'] - df_splits_locs['region_coord1']

# Number of same alignments
df_splits_locs = df_splits_locs.groupby(df_splits_locs.columns.tolist(),as_index=False).size()

# Split ID
list_of_ids = []
id_num = 0
cerr = 10

for k in range(0, len(df_splits_locs.index), 1):
    
    if(k == 0):
        id_num = 0
        list_of_ids.append(id_num)
        continue
    
    # Next contig ID - new split
    elif( df_splits_locs.loc[k, 'contig_id'] != df_splits_locs.loc[k-1, 'contig_id'] ):
        id_num = id_num + 1
        list_of_ids.append(id_num)
        continue
    
    # Overlapping region/cut - same split
    elif(df_splits_locs.loc[k, 'region_length'] == 0):
        if(           # Previous is also cut
            ((df_splits_locs.loc[k-1, 'region_length'] != 0) and \
           (df_splits_locs.loc[k, 'region_coord1'] < df_splits_locs.loc[k-1, 'region_coord2'])) \
           \
           or \
           # Previous is interval
           ((df_splits_locs.loc[k-1, 'region_length'] == 0) and \
               (df_splits_locs.loc[k, 'region_coord1'] - df_splits_locs.loc[k-1, 'region_coord2'] <= cerr))\
          ):
            id_num = id_num
        else:
            id_num = id_num + 1

    # Overlapping region - same split
    elif( df_splits_locs.loc[k, 'region_coord1'] <= df_splits_locs.loc[k-1, 'region_coord2'] ): #<
        id_num = id_num
        
    # New split
    else:
        id_num = id_num + 1

    list_of_ids.append(id_num)

# Splits with IDs
df_splits_locs['split_id'] = list_of_ids
df_splits_locs = df_splits_locs.sort_values(['split_id', 'size', 'region_length'], ascending=True).reset_index(drop=True)
#print(df_splits_locs)

# Select split
# which has the largest size/multiple aligns (if such exists), otherwise the largest region length
df_splits_final_coords = df_splits_locs.groupby(['split_id']).last().reset_index()

#print(df_splits_final_coords)


##################
# Execute splits #
##################

# Function to write fasta split file
def functionSplit(scontigid, scoord1, scoord2, snewcontigid):
    record = SeqIO.read(folder_loc_orientation + "/" + scontigid + ".fa", "fasta")
    sequences_string = record.seq
    sequences_id = record.id

    if(scoord2 == 0):
        part = sequences_string[(scoord1-1):]
    else:
        part = sequences_string[(scoord1-1):scoord2]
    
    seq_split = SeqRecord(
        Seq(part),
        id=snewcontigid,
        #name="",
        description="",
    )
    
    SeqIO.write(seq_split, folder_loc_orientation + "/" + snewcontigid + ".fa", "fasta")
    
    return(snewcontigid)


# Newly created contigs (after splits)
contigs_splitted_list = []

new_contig = 0
for k in range(0, len(df_splits_final_coords), 1):
    
    contig_id_k = df_splits_final_coords.loc[k, "contig_id"]
    split_id_k = str(df_splits_final_coords.loc[k, "split_id"])
    region1_k = df_splits_final_coords.loc[k, "region_coord1"]
    region2_k = df_splits_final_coords.loc[k, "region_coord2"]
    
    # If it's a cut
    if(region1_k == region2_k):
        region2_k = region2_k + 1
    
    # If this is the last split in the list
    if(k == len(df_splits_final_coords)-1):
        if(new_contig == 0):
            contigs_splitted_list.append(functionSplit(contig_id_k, 1, region1_k, contig_id_k + "00"))
        contigs_splitted_list.append(functionSplit(contig_id_k, region2_k, 0, contig_id_k + split_id_k))
        break
    
    # If this is the first split of this contig
    if(new_contig == 0):
        contigs_splitted_list.append(functionSplit(contig_id_k, 1, region1_k, contig_id_k + "00"))
    
    contig_id_knext = df_splits_final_coords.loc[k+1, "contig_id"]
    region1_knext = df_splits_final_coords.loc[k+1, "region_coord1"]
    
    # If the next split is for the same contig
    if(contig_id_k == contig_id_knext):
        contigs_splitted_list.append(functionSplit(contig_id_k, region2_k, region1_knext, contig_id_k + split_id_k))
        new_contig = 1
    else:
        contigs_splitted_list.append(functionSplit(contig_id_k, region2_k, 0, contig_id_k + split_id_k))
        new_contig = 0

#####
# If two splits cut at the same location:
# Select based on size (i.e. how many agree)
df_splits_final_coords = df_splits_final_coords.drop(df_splits_final_coords[(df_splits_final_coords["region_coord2"].shift(1) == df_splits_final_coords["region_coord1"]) & (df_splits_final_coords["size"].shift(1) >= df_splits_final_coords["size"])].index).reset_index(drop=True)
#####

# Contigs that were split are removed and only the splitted parts remain
contigs_remove = df_splits_final_coords['contig_id'].unique().tolist()

#print(df_splits_final_coords)


# Read input file
nano_positions = []
l = 0
for elem in nano_ids:
    nano_positions.append(pd.read_csv('position_step/npositions_'+str(elem)+'.csv', usecols= ['sseqid','pident','length','qstart','qend','sstart','send','qlen','slen']))
    l = l + 1


##################################
# Alignments of splitted contigs #
##################################
print("Updating position information of splitted contigs...")

# Query: nanos (rm, trf)
# Subject: rextal contigs

cols_outfmt6 = 'qseqid sseqid pident length qstart qend sstart send evalue qlen slen'.strip().split(' ')
df_splitted_aligns = pd.DataFrame()
df_new_aligns = pd.DataFrame()

for k in range(0, len(contigs_splitted_list), 1):
    
    # Blastn
    blast_outfmt6 = subprocess.check_output(f"blastn -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen\" -perc_identity 80 -subject {folder_loc_orientation}/{contigs_splitted_list[k]}.fa -query {folder_loc_orientation}/oriented_nanos_rm_trf.fa", shell=True, encoding='UTF-8')
    try:
        df_splitted_aligns = pd.read_table(StringIO(blast_outfmt6), header=None)  
    except:
        # No alignments
        continue
    
    df_splitted_aligns.columns = cols_outfmt6
    
    # Clean alignments: Remove reverse alignments
    df_splitted_aligns = df_splitted_aligns.loc[(df_splitted_aligns['qstart'] <= df_splitted_aligns['qend']) & (df_splitted_aligns['sstart'] <= df_splitted_aligns['send'])]
    
    # Select alignments: Keep only alignments with evalue=0 (if specified)
    if(flag_eval):
        df_splitted_aligns = df_splitted_aligns.loc[(df_splitted_aligns['evalue'] == 0.0)]
        df_splitted_aligns = df_splitted_aligns.reset_index(drop=True)

    # Sort alignments by nano name and nano coordinates
    df_splitted_aligns = df_splitted_aligns.sort_values(['qseqid', 'qstart'], ascending=True).reset_index(drop=True)
        
    # Dataframe that has first and last alignment of each contig
    df_nano_first = df_splitted_aligns.groupby(['qseqid']).first().reset_index()
    df_nano_last = df_splitted_aligns.groupby(['qseqid']).last().reset_index()
    df_nano_full = pd.concat([df_nano_first, df_nano_last]).reset_index(drop=True)
    df_nano_full = df_nano_full.sort_values(['qseqid', 'qstart'], ascending=True).reset_index(drop=True)
    
    # Append first and last of this contig for every nano to the total list
    df_new_aligns = df_new_aligns.append(df_nano_full)

# Sort according to nano id
df_new_aligns = df_new_aligns.sort_values(['qseqid', 'qstart'], ascending=True).reset_index(drop=True)
# Drop evalue column
df_new_aligns = df_new_aligns.drop(columns = ['evalue'])


# Add splitted contigs
for i in range(len(nano_ids)):
    nano_positions[i] = nano_positions[i].append((df_new_aligns.loc[df_new_aligns['qseqid'] == nano_ids[i]]).drop(columns = ['qseqid']))
    nano_positions[i] = nano_positions[i].sort_values(['qstart'], ascending=True).reset_index(drop=True)


# Remove non-splitted contigs and get order of each nano
contig_order2 = []

for i in range(len(nano_ids)):
    # Remove non-splitted
    nano_positions[i] = nano_positions[i].loc[~nano_positions[i]['sseqid'].isin(contigs_remove)]
    nano_positions[i] = nano_positions[i].reset_index(drop=True)


for i in range(len(nano_ids)):
    nano_positions[i] = nano_positions[i].sort_values(['qstart','sseqid'], ascending=True).reset_index(drop=True)
    nano_positions[i] = nano_positions[i].drop_duplicates(subset=['qstart','qend'], inplace=False).reset_index(drop=True)
    
    for cnt_triple_check in range(6): #
        # Remove contained contigs
        nano_positions[i] = nano_positions[i].drop(nano_positions[i][(nano_positions[i]["qend"].shift(1) >= nano_positions[i]["qend"])].index).reset_index(drop=True)

        # Remove extended contained contigs
        nano_positions[i]['extension_val'] = nano_positions[i]["qend"]+nano_positions[i]["slen"]-nano_positions[i]["send"]
        # Compare (i) the extended end coordinates between (ii) different contigs
        nano_positions[i] = nano_positions[i].drop(nano_positions[i][  (nano_positions[i]["extension_val"].shift(1) >= nano_positions[i]["extension_val"]) & (nano_positions[i]["sseqid"].shift(1) != nano_positions[i]["sseqid"]) ].index).reset_index(drop=True)
        nano_positions[i] = nano_positions[i].drop(columns = 'extension_val').reset_index(drop=True)
        
        ###
        # Remove contained backwards (it happens only when qstart is same in both)
        nano_positions[i] = nano_positions[i].drop(nano_positions[i][(nano_positions[i]["qstart"].shift(-1) == nano_positions[i]["qstart"]) & (nano_positions[i]["qend"].shift(-1) >= nano_positions[i]["qend"])].index).reset_index(drop=True)
        
        #Should do the same to remove extended contained backwards
        nano_positions[i]['extension_val'] = nano_positions[i]["qstart"]-nano_positions[i]["sstart"]
        nano_positions[i]['extension_val2'] = nano_positions[i]["qend"]+nano_positions[i]["slen"]-nano_positions[i]["send"]
        nano_positions[i] = nano_positions[i].drop(nano_positions[i][  (nano_positions[i]["extension_val"] >= nano_positions[i]["extension_val"].shift(-1)) & (nano_positions[i]["extension_val2"] <= nano_positions[i]["extension_val2"].shift(-1)) & (nano_positions[i]["sseqid"] != nano_positions[i]["sseqid"].shift(-1)) ].index).reset_index(drop=True)
        nano_positions[i] = nano_positions[i].drop(columns = 'extension_val').reset_index(drop=True)
        nano_positions[i] = nano_positions[i].drop(columns = 'extension_val2').reset_index(drop=True)
        ###
    
    # Make doubles the single ones
    nano_positions[i] = nano_positions[i].sort_values(['qstart'], ascending=True).reset_index(drop=True)
    df_nano1 = nano_positions[i].groupby(['sseqid']).first().reset_index()
    df_nano2 = nano_positions[i].groupby(['sseqid']).last().reset_index()
    nano_positions[i] = pd.concat([df_nano1, df_nano2]).reset_index(drop=True)
    nano_positions[i] = nano_positions[i].sort_values(['qstart'], ascending=True).reset_index(drop=True)

    # Get order
    df_nano_first = nano_positions[i].groupby(['sseqid']).first().reset_index()
    df_nano_first = df_nano_first.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    list_order_i = df_nano_first['sseqid'].tolist()

    contig_order2.append(list_order_i)
    

#######################################
# Repeat steps to confirm final order #
#######################################
##################################################
# Find the order of the contigs for the assembly #
##################################################

# The set of all unique contigs
unique_contigs = {x for y in contig_order2 for x in y}

# Dictionary that maps all unique contigs to integers
dict_unique_contigs = dict([(y,x+1) for x,y in enumerate(sorted(set(unique_contigs)))])
#print(dict_unique_contigs)

# Create the directed graph
n = len(unique_contigs)
order_graph = [[] for i in range(n + 1)]

# Iterate through every nano's contig order
for k in range(0, len(contig_order2), 1):
    #print(contig_order2[k])
    for i in range(0, len(contig_order2[k])-1, 1):
        u = dict_unique_contigs[contig_order2[k][i]]
        v = dict_unique_contigs[contig_order2[k][i+1]]
        #print(str(u) + " - " + str(v))
        # Add an edge only if it doesn't already exist
        if((v in order_graph[u]) == False):
            order_graph[u].append(v)
#print(order_graph)

# Use networkX package for the graph
import networkx as nx

# Create graph
DG = nx.DiGraph()
# Add nodes
DG.add_nodes_from(list(range(1,len(order_graph))))
# Add edges
for i in range(1, len(order_graph)):
    #print(order_graph[i])
    for k in range(0, len(order_graph[i])):
        DG.add_edge(i, order_graph[i][k])

#########################
# Detect and fix cycles #
#########################
# Find cycles
#len(list(nx.simple_cycles(DG)))
#list(nx.simple_cycles(DG))

# Remove erroneous edges to fix a cycle
#DG.remove_edge()

#########################
#########################

# Find the longest path in the graph
assembly_contig_order2 = nx.dag_longest_path(DG)

# Get the contig IDs
inv_dict_unique_contigs = {x: y for y, x in dict_unique_contigs.items()}

# Contig order for the assembly
assembly_contig_order_ids2 = [inv_dict_unique_contigs[x] for x in assembly_contig_order2]
#print(assembly_contig_order_ids2)


# Write to csv (keep all in a folder)
if not os.path.exists(folder_loc_correction):
    os.makedirs(folder_loc_correction)

l = 0
for elem in nano_ids:
    nano_positions[l].to_csv(str(folder_loc_correction)+'/npositions_'+str(elem)+'.csv', index=False)
    l = l + 1


# Write assembly_contig_order_ids2
with open(str(folder_loc_correction)+'/assembly_order.csv', 'w', newline='', encoding='utf-8') as f: 
    write = csv.writer(f)
    write.writerow(['id'])
    for elem in assembly_contig_order_ids2:
        write.writerow([elem])


# Write contig_order2
with open(str(folder_loc_correction)+'/contig_order2.csv', 'w', newline='', encoding='utf-8') as f: 
    write = csv.writer(f)
    write.writerow(['order'])
    for elem in contig_order2:
        write.writerow([elem])


print("The NPGREAT Correction step finished.")


