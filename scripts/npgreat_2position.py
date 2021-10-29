######################
######################
# NPGREAT - POSITION #
######################
######################
# REXTAL sequence contigs are aligned and positioned relative to the repeat-masked nanopore sequence reads.

###########################
import pandas as pd
import subprocess
from io import StringIO
import os
import sys
import csv

folder_loc_orientation = "orientation_step"
folder_loc_position = 'position_step'

print("- The NPGREAT Position step begins...")


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


############################
# Calculate the alignments #
############################
print("Calculating the alignments...")

# Query: nanos (rm, trf)
# Subject: rextal contigs
blast_outfmt6 = subprocess.check_output(f"blastn -subject {folder_loc_orientation}/oriented_rextal_contigs.fa -query {folder_loc_orientation}/oriented_nanos_rm_trf.fa -perc_identity 80 -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen\"", shell=True, encoding='UTF-8')
try:
    # Read blast output in outfmt=6 format (with no header)
    df = pd.read_table(StringIO(blast_outfmt6), header=None)  
except:
    # No alignments
    print("No alignments!")

# Columns of outfmt=6
cols_outfmt6 = 'qseqid sseqid pident length qstart qend sstart send evalue qlen slen'.strip().split(' ')
df.columns = cols_outfmt6

# Clean alignments: 
# Remove reverse alignments
df_aligns = df.loc[(df['qstart'] <= df['qend']) & (df['sstart'] <= df['send'])]

# Keep alignments with evalue = 0.0 (if specified)
if(flag_eval):
    df_aligns = df_aligns.loc[(df['evalue'] == 0.0)]


#######################################################################
# Group the alignments by each nano, sort them, extract position info #
#######################################################################
print("Extracting position information...")

# Sort qseqid column, Group by 'qseqid', Sort within each group by 'qstart' (and then 'sstart')
df_aligns = df_aligns.sort_values(['qseqid'], ascending=True).groupby(['qseqid'], sort=False)     .apply(lambda x: x.sort_values(['qstart','sstart'], ascending=True)).reset_index(drop=True)

# Group by qseqid, i.e. each nano
grouped_aligns = df_aligns.groupby("qseqid")

contig_order = []
nano_positions = []
nano_ids = []
contig_aligns = []
group_iter = -1
# Get each group (nano)
for name, group in grouped_aligns:
    group_iter = group_iter+1
    
    # Nano IDs
    df_nano0 = grouped_aligns.get_group(name)
    nano_ids.append(name)
    df_nano0 = df_nano0.drop(columns=['qseqid', 'evalue'])
    df_nano = df_nano0
        
    # Dataframe that has first and last alignment of each contig
    df_nano_first = df_nano.groupby(['sseqid']).first().reset_index()    
    df_nano_first = df_nano_first.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    
    df_nano_last = df_nano.groupby(['sseqid']).last().reset_index()
    df_nano_last = df_nano_last.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    
    df_nano_full = pd.concat([df_nano_first, df_nano_last]).reset_index(drop=True)
    df_nano_full = df_nano_full.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    
    # Clean alignments
    # Remove same contigs (it also removes internally duplicate alignments)
    df_nano_full = df_nano_full.sort_values(['qstart','sseqid'], ascending=True).reset_index(drop=True)
    df_nano_full = df_nano_full.drop_duplicates(subset=['qstart','qend'], inplace=False).reset_index(drop=True)
    
    # Remove contained contigs
    df_nano_full = df_nano_full.drop(df_nano_full[(df_nano_full["qend"].shift(1) >= df_nano_full["qend"])].index).reset_index(drop=True)
    
    # Remove extended contained contigs
    df_nano_full['extension_val'] = df_nano_full["qend"]+df_nano_full["slen"]-df_nano_full["send"]
    # Compare (i) the extended end coordinates between (ii) different contigs that (iii) are not possible splits
    df_nano_full = df_nano_full.drop(df_nano_full[  (df_nano_full["extension_val"].shift(1) >= df_nano_full["extension_val"]) & (df_nano_full["sseqid"].shift(1) != df_nano_full["sseqid"]) & (df_nano_full["sseqid"].shift(1) == df_nano_full["sseqid"].shift(2))  ].index).reset_index(drop=True)
    df_nano_full = df_nano_full.drop(columns = 'extension_val').reset_index(drop=True)
    
    # Select good aligns
    # Keep the align with the longest length (splits are not affected by this selection - they depend on internal alignments)
    uniq_contigs = df_nano_full['sseqid'].unique()
    for cnt_uniq in range(len(uniq_contigs)):
        df_eachcontig = df_nano_full.drop(df_nano_full[(df_nano_full["sseqid"] != uniq_contigs[cnt_uniq])].index)
        contig_indexes = df_eachcontig.index.tolist()
        if(len(contig_indexes) == 2 and contig_indexes[1] - contig_indexes[0] > 1):
            if(df_eachcontig.loc[contig_indexes[0],"length"] > df_eachcontig.loc[contig_indexes[1],"length"]):
                df_nano_full = df_nano_full.drop(contig_indexes[1])
            else:
                df_nano_full = df_nano_full.drop(contig_indexes[0])
    df_nano_full = df_nano_full.reset_index(drop=True)
       
    # Make doubles the single ones
    df_nano_full = df_nano_full.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    df_nano1 = df_nano_full.groupby(['sseqid']).first().reset_index()
    df_nano2 = df_nano_full.groupby(['sseqid']).last().reset_index()
    df_nano_full = pd.concat([df_nano1, df_nano2]).reset_index(drop=True)
    df_nano_full = df_nano_full.sort_values(['qstart'], ascending=True).reset_index(drop=True)
    
    # Contig positions according to the nano
    nano_positions.append(df_nano_full)
    
    # Contig order according to the nano
    contig_order.append(df_nano_first['sseqid'].tolist())
    
    # Contig alignments internally
    contig_aligns.append(df_nano)


# Write to csv (keep all in a folder)
if not os.path.exists(folder_loc_position):
    os.makedirs(folder_loc_position)

# Write nano_positions
l = 0
for elem in nano_ids:
    nano_positions[l].to_csv(str(folder_loc_position)+'/npositions_'+str(elem)+'.csv', index=False)
    l = l + 1

# Write nano_ids
with open(str(folder_loc_position)+'/nano_ids.csv', 'w', newline='', encoding='utf-8') as f: 
    write = csv.writer(f)
    write.writerow(['id'])
    for elem in nano_ids:
        write.writerow([elem])

# Write contig_aligns
l = 0
for elem in nano_ids:
    contig_aligns[l].to_csv(str(folder_loc_position)+'/caligns_'+str(elem)+'.csv', index=False)
    l = l + 1

# Write contig_order
with open(str(folder_loc_position)+'/contig_order1.csv', 'w', newline='', encoding='utf-8') as f: 
    write = csv.writer(f)
    write.writerow(['order'])
    for elem in contig_order:
        write.writerow([elem])


print("The NPGREAT Position step finished.")


