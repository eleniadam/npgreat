###############################
###############################
# NPGREAT - REGION EXTRACTION #
###############################
###############################
# The alignment of REXTAL contigs with nanopore reads, yields two possibilities for neighboring REXTAL contigs: 
# their overlap or a gap between them. In the region extraction/connector segments step, the overlapping REXTAL contigs are merged 
# and nanopore read segments that can bridge gaps between neighboring contigs, are identified and extracted.

###########################
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

folder_loc_orientation = "orientation_step"
folder_loc_position = "position_step"
folder_loc_correction = "correction_step"

print("- The NPGREAT Region Extraction/Connector Segments step begins...")


# Nanopore IDs
# Read nano_ids
nano_ids = []
with open(str(folder_loc_position)+'/nano_ids.csv', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    for elem in (list(reader)[1:]):
        nano_ids.append(''.join(elem))

# Contig first and last alignments for each nano
# Read nano_positions
nano_positions = []
l = 0
for elem in nano_ids:
    nano_positions.append(pd.read_csv(str(folder_loc_correction)+'/npositions_'+str(elem)+'.csv', usecols= ['sseqid','pident','length','qstart','qend','sstart','send','qlen','slen']))
    l = l + 1

# Total contig order for the assembly
# Read assembly_contig_order_ids2
assembly_contig_order_ids2 = []
with open(str(folder_loc_correction)+'/assembly_order.csv', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    for elem in (list(reader)[1:]):
        assembly_contig_order_ids2.append(''.join(elem))

# Contig order in each nano
# Read contig_order2
contig_order2 = []
with open(str(folder_loc_correction)+'/contig_order2.csv', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)
    
    for elem in (list(reader)[1:]):
        contig_order2.append((((''.join(elem)[1:-1]).replace("'","")).replace(" ","")).split(","))


# Call blastn and obtain info
def run_blastn(nano_list, contigs_list, nano_coords):
    # Query: nano
    # Subject: rextal contigs
    
    cols_outfmt6 = 'qseqid sseqid pident length qstart qend sstart send evalue qlen slen'.strip().split(' ')
    df_curr_aligns = pd.DataFrame()
    df_aligns = pd.DataFrame()
    
    for k in range(0, len(contigs_list), 1):
        # Blastn
        blast_outfmt6 = subprocess.check_output(f"blastn -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen\" -perc_identity 80 -subject {folder_loc_orientation}/{contigs_list[k]}.fa -query {folder_loc_orientation}/{nano_list}.fa -query_loc {nano_coords[0]}-{nano_coords[1]}", shell=True, encoding='UTF-8')
        try:
            df_curr_aligns = pd.read_table(StringIO(blast_outfmt6), header=None)  
        except:
            # No alignments
            continue

        df_curr_aligns.columns = cols_outfmt6

        # Clean alignments: Remove reverse alignments
        df_curr_aligns = df_curr_aligns.loc[(df_curr_aligns['qstart'] <= df_curr_aligns['qend']) & (df_curr_aligns['sstart'] <= df_curr_aligns['send'])]

        # Sort alignments by nano name and nano coordinates
        df_curr_aligns = df_curr_aligns.sort_values(['qseqid', 'evalue'], ascending=True).reset_index(drop=True)

        ###
        df_curr_aligns = df_curr_aligns.sort_values(['qseqid', 'length'], ascending=False).reset_index(drop=True)
        
        # Dataframe that has first and last alignment of each contig
        df_nano_first = df_curr_aligns.groupby(['qseqid']).first().reset_index()
        df_nano_full = pd.concat([df_nano_first]).reset_index(drop=True)
        df_nano_full = df_nano_full.sort_values(['qseqid', 'qstart'], ascending=True).reset_index(drop=True)

        # Append first and last of this contig for every nano to the total list
        df_aligns = df_aligns.append(df_nano_full)

    if(df_aligns.empty):
        return(df_aligns)
    
    # Sort according to nano id
    df_aligns = df_aligns.sort_values(['qseqid', 'qstart'], ascending=True).reset_index(drop=True)
    
    return(df_aligns)


# Region extraction function
# Find: (1) Connectors (incl. extensions), (2) Overlaps (incl. merging), (3) further investigation - {for 3 repeat code to reach[1,2]}
def region_extraction(nano_ids, assembly_contig_order_ids2, nano_positions):
    
    # Dataframe containing all regions' coordinates
    df_regions = pd.DataFrame(columns = ['category', 'contig1', 'contig2', 'contigcoord', 'nano', 'nanocoord1', 'nanocoord2', 'pid'])
    
    for k in range(len(nano_ids)):

        ###############
        #  Extension  #
        ###############
        ###############
        # Extension-1
        i = 0
        if((nano_positions[k].loc[i, "sseqid"] == assembly_contig_order_ids2[0])): 
            #telom: p-end

            nano_id0 = nano_ids[k]
            rextal_id0 = nano_positions[k].loc[i, "sseqid"]
            pid0 = nano_positions[k].loc[i, "pident"]
            coord_nano0 = nano_positions[k].loc[i, "qstart"]
            coord_rextal0 = nano_positions[k].loc[i, "sstart"]

            # Calculate region
            nuc_end0 = coord_nano0 - coord_rextal0
            if(nuc_end0 >= 0):
                df_regions = df_regions.append({'contig1': "x", 'contig2': rextal_id0, 'category': "extension_1", 'nano': nano_id0, 'nanocoord1': 1, 'nanocoord2': nuc_end0, 'contigcoord': 0, 'pid': pid0}, ignore_index=True)

        ###############
        # Extension-2
        i = len(nano_positions[k])-1
        if(nano_positions[k].loc[i, "sseqid"] == assembly_contig_order_ids2[len(assembly_contig_order_ids2)-1]):
            #telom: q-end

            nano_id0 = nano_ids[k]
            rextal_id0 = nano_positions[k].loc[i, "sseqid"]
            pid0 = nano_positions[k].loc[i, "pident"]
            len_nano0 = nano_positions[k].loc[i, "qlen"]
            coord_nano0 = nano_positions[k].loc[i, "qend"]
            len_rextal0 = nano_positions[k].loc[i, "slen"]
            coord_rextal0 = nano_positions[k].loc[i, "send"]

            # Calculate region
            nuc_beg0 = coord_nano0 + len_rextal0 - coord_rextal0
            nuc_end0 = len_nano0
            if(nuc_end0 >= 0):

                # If there is an extension
                if(nuc_end0 > nuc_beg0):
                    df_regions = df_regions.append({'contig1': rextal_id0, 'contig2': "x", 'category': "extension_2", 'nano': nano_id0, 'nanocoord1': nuc_beg0, 'nanocoord2': nuc_end0, 'contigcoord': 0, 'pid': pid0}, ignore_index=True)
                else:
                    df_regions = df_regions.append({'contig1': rextal_id0, 'contig2': "x", 'category': "extension_2", 'nano': "x", 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': 0, 'pid': 0}, ignore_index=True)

        #############################################################################################  
        for i in range(1, len(nano_positions[k])-1, 2):

            nano_id = nano_ids[k]
            len_nano_q = nano_positions[k].loc[i, "qlen"]
            
            rextal_id1 = nano_positions[k].loc[i, "sseqid"]
            rextal_id2 = nano_positions[k].loc[i+1, "sseqid"]
            
            len_align1 = nano_positions[k].loc[i, "length"]
            len_align2 = nano_positions[k].loc[i+1, "length"]

            lpid = nano_positions[k].loc[i, "pident"]
            rpid = nano_positions[k].loc[i+1, "pident"]
            ave_pid = (lpid + rpid)/2

            coord_rextal1 = nano_positions[k].loc[i, "send"]
            coord_rextal2 = nano_positions[k].loc[i+1, "sstart"]
            coord_nano1 = nano_positions[k].loc[i, "qend"]
            coord_nano2 = nano_positions[k].loc[i+1, "qstart"]
            len_r1 = nano_positions[k].loc[i, "slen"]
            len_r2 = nano_positions[k].loc[i+1, "slen"]

            # Calculate extended coordinates of nano
            lcoord = coord_nano1 + len_r1 - coord_rextal1 +1
            rcoord = coord_nano2 - coord_rextal2
            
            # Check coordinates, whether they correspond to a region or they overlap
            
            ###############
            #  Connector  #
            ###############
            if (lcoord <= rcoord):
                #print("connector")
                # If they correspond to a region - extraction
                df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "connector", 'nano': nano_id, 'nanocoord1': lcoord, 'nanocoord2': rcoord, 'contigcoord': 0, 'pid': ave_pid}, ignore_index=True)
                 
            ###############
            #   Overlap   #
            ###############
            elif (coord_nano1 > coord_nano2):
                #print("overlap")
                # If the original coordinates overlap
                # Calculate number of overlapping bases
                ov = coord_nano1 - coord_nano2 + 2

                ###
                ov_denA = nano_positions[k].loc[i+1, "sstart"]
                ov_denB = len_r1 - nano_positions[k].loc[i, "send"]
                ###
                
                # Keep whole the rextal which has the longest alignment 
                # with the nanopore according to blastn output
                # Cut the beginning (or end correspondingly) of the other rextal
                # Merge them
                # Overlap-1
                if (ov_denA >= ov_denB):
                    #print("overlap-1")
                    nuc_beg = coord_rextal2 + ov + (len_r1 - coord_rextal1)
                    
                    df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_1", 'nano': nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_beg, 'pid': ave_pid}, ignore_index=True)
                
                # Overlap-2
                else:
                    #print("overlap-2")
                    nuc_end_orig = coord_rextal1 - ov - coord_rextal2
                    nuc_end = len_r1 - nuc_end_orig
                    nuc_end = -nuc_end
                    
                    df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_2", 'nano': nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_end, 'pid': ave_pid}, ignore_index=True)
            
            #############################################################################################    
            #########################
            # Further investigation #
            #########################
            else:
                # Call blastn with unmasked nano
                #print("further investigation")
                
                fi_nano_id = nano_ids[k]
                fi_rextal_ids = [rextal_id1, rextal_id2]
                
                chcoordnano1 = coord_nano1 - 1000
                chcoordnano2 = coord_nano2 + 1000
                
                if(chcoordnano1 <= 0):
                    chcoordnano1 = 1
                if(chcoordnano2 > len_nano_q):
                    chcoordnano2 = len_nano_q
                
                fi_nano_coords = [chcoordnano1, chcoordnano2]
                
                # Call blastn
                aligns_investigate = run_blastn(fi_nano_id, fi_rextal_ids, fi_nano_coords)

                if (aligns_investigate.empty or len(aligns_investigate['sseqid'].unique()) == 1):
                    #print("No further investigation possible.")
                    ("")
                else:
                    rextal_id1 = aligns_investigate.loc[0, "sseqid"]
                    rextal_id2 = aligns_investigate.loc[1, "sseqid"]

                    len_align1 = aligns_investigate.loc[0, "length"]
                    len_align2 = aligns_investigate.loc[1, "length"]

                    lpid = aligns_investigate.loc[0, "pident"]
                    rpid = aligns_investigate.loc[1, "pident"]
                    ave_pid = (lpid + rpid)/2

                    coord_rextal1 = aligns_investigate.loc[0, "send"]
                    coord_rextal2 = aligns_investigate.loc[1, "sstart"]
                    coord_nano1 = aligns_investigate.loc[0, "qend"]
                    coord_nano2 = aligns_investigate.loc[1, "qstart"]
                    len_r1 = aligns_investigate.loc[0, "slen"]
                    len_r2 = aligns_investigate.loc[1, "slen"]
                    
                    lcoord = coord_nano1 + len_r1 - coord_rextal1 +1
                    rcoord = coord_nano2 - coord_rextal2

                    # Connector
                    if (lcoord <= rcoord):
                        df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "connector", 'nano': fi_nano_id, 'nanocoord1': lcoord, 'nanocoord2': rcoord, 'contigcoord': 0, 'pid': ave_pid}, ignore_index=True)
                    
                    # Overlap
                    elif (coord_nano1 > coord_nano2):
                        ov = coord_nano1 - coord_nano2 + 2

                        ov_denA = aligns_investigate.loc[1, "sstart"]
                        ov_denB = len_r1 - aligns_investigate.loc[0, "send"]

                        # Overlap-1
                        if (ov_denA >= ov_denB):
                            nuc_beg = coord_rextal2 + ov + (len_r1 - coord_rextal1)

                            df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_1", 'nano': fi_nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_beg, 'pid': ave_pid}, ignore_index=True)
                        
                        # Overlap-2
                        else:
                            nuc_end_orig = coord_rextal1 - ov - coord_rextal2
                            nuc_end = len_r1 - nuc_end_orig
                            nuc_end = -nuc_end
                            
                            df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_2", 'nano': fi_nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_end, 'pid': ave_pid}, ignore_index=True)
                    
                    ##########
                    else:
                        coord_nano1 = coord_nano1 +len_r1 -coord_rextal1
                        coord_nano2 = coord_nano2 -coord_rextal2
                        ave_pid = 0 #Set as the least possible option
                        if (coord_nano1 > coord_nano2):
                            ov = coord_nano1 - coord_nano2 + 2

                            ov_denA = aligns_investigate.loc[1, "sstart"]
                            ov_denB = len_r1 - aligns_investigate.loc[0, "send"]

                            # Overlap-1
                            if (ov_denA >= ov_denB):
                                nuc_beg = ov 

                                df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_1", 'nano': fi_nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_beg, 'pid': ave_pid}, ignore_index=True)
                            
                            # Overlap-2
                            else:
                                nuc_end = ov
                                nuc_end = -nuc_end

                                df_regions = df_regions.append({'contig1': rextal_id1, 'contig2': rextal_id2, 'category': "overlap_2", 'nano': fi_nano_id, 'nanocoord1': 0, 'nanocoord2': 0, 'contigcoord': nuc_end, 'pid': ave_pid}, ignore_index=True)

                        else:
                            #print("No further investigation possible.-fi")
                            ("")
                    ##########

    return(df_regions)


# Execute Region Extraction
df_regions_npgreat = region_extraction(nano_ids, assembly_contig_order_ids2, nano_positions)

print("The NPGREAT Region Extraction/Connector Segments step finished.")



#########################
#########################
# NPGREAT - GAP FILLING #
#########################
#########################
# For each gap between REXTAL contigs, several nanopore segments may be available to bridge it. 
# In order to fill the gap, we select the segment that has the highest average percent identity with the flanking contigs.

print("- The NPGREAT Gap Filling step begins...")
#########################

df_regions_npgreat = df_regions_npgreat.sort_values(['contig1', 'contig2', 'pid'], ascending=True).reset_index(drop=True)

df_align_maxpid = df_regions_npgreat.groupby(['contig1', 'contig2']).last().reset_index()
df_aligns_maxpids = pd.concat([df_align_maxpid]).reset_index(drop=True)

print("The NPGREAT Gap Filling step finished.")



#########################
#########################
# NPGREAT - COMBINATION #
#########################
#########################
# In the final step, we combine according to their order the REXTAL contigs, the merged REXTAL contigs and the
# nanopore selected segments that connect as well as extend them. The result is the assembled sequence.

print("- The NPGREAT Combination step begins...")
#########################

data_count = 0
bovl = False
df_combined = pd.DataFrame(columns = ['piece', 'seq_id', 'coord1', 'coord2'])

for i in range(0, len(assembly_contig_order_ids2)):
    
    # First contig extension (extension-1)
    if(i == 0):
        obs = df_aligns_maxpids.loc[df_aligns_maxpids['contig2'] == assembly_contig_order_ids2[i]]
        obs = obs.reset_index(drop=True)
        if( (obs.empty) or (obs.loc[0, "nano"] == 'x') ):
            df_combined = df_combined.append({'piece': data_count, 'seq_id': 'x', 'coord1': 0, 'coord2': 0}, ignore_index=True)
        else:
            df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "nano"], 'coord1': obs.loc[0, "nanocoord1"], 'coord2': obs.loc[0, "nanocoord2"]}, ignore_index=True)    
        data_count = data_count + 1
    
    # Last contig extension (extension-2)
    if(i == len(assembly_contig_order_ids2)-1):
        obs = df_aligns_maxpids.loc[df_aligns_maxpids['contig1'] == assembly_contig_order_ids2[i]]
        obs = obs.reset_index(drop=True)
        if( (obs.empty) or (obs.loc[0, "nano"] == 'x') ):
            if(bovl == False):
                df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig1"], 'coord1': 0, 'coord2': 0}, ignore_index=True)
                data_count = data_count + 1
            df_combined = df_combined.append({'piece': data_count, 'seq_id': 'x', 'coord1': 0, 'coord2': 0}, ignore_index=True)
        else:
            if(bovl == False):
                df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig1"], 'coord1': 0, 'coord2': 0}, ignore_index=True)
                data_count = data_count + 1
            df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "nano"], 'coord1': obs.loc[0, "nanocoord1"], 'coord2': obs.loc[0, "nanocoord2"]}, ignore_index=True)    
        break
    
    # In between: connector or overlap
    obs = df_aligns_maxpids.loc[(df_aligns_maxpids['contig1'] == assembly_contig_order_ids2[i]) & (df_aligns_maxpids['contig2'] == assembly_contig_order_ids2[i+1])]
    obs = obs.reset_index(drop=True)
    obs_category = obs.loc[0, "category"]
    
    # Connector
    if(obs_category == 'connector'):
        if(bovl == False):
            df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig1"], 'coord1': 0, 'coord2': 0}, ignore_index=True) # whole
            data_count = data_count + 1
        df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "nano"], 'coord1': obs.loc[0, "nanocoord1"], 'coord2': obs.loc[0, "nanocoord2"]}, ignore_index=True)
        data_count = data_count + 1
        bovl = False
        
    # Overlap-1
    elif(obs_category == 'overlap_1'):
        
        if(bovl == False):
            df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig1"], 'coord1': 0, 'coord2': 0}, ignore_index=True)
            data_count = data_count + 1
            
        df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig2"], 'coord1': obs.loc[0, "contigcoord"], 'coord2': 0}, ignore_index=True)
        data_count = data_count + 1
        bovl = True
    # Overlap-2
    elif(obs_category == 'overlap_2'):
        if(bovl == False):
            df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig1"], 'coord1': 0, 'coord2': obs.loc[0, "contigcoord"]}, ignore_index=True)
            data_count = data_count + 1
        else:
            ##########
            df_combined.loc[data_count-1, "coord2"] = obs.loc[0, "contigcoord"]
            ###########
        df_combined = df_combined.append({'piece': data_count, 'seq_id': obs.loc[0, "contig2"], 'coord1': 0, 'coord2': 0}, ignore_index=True)

        data_count = data_count + 1
        bovl = True


#######################
# Execute combination #
#######################

# Function to write fasta file
def functionCombine(df_pieces):
    
    assembly_id = "assembly_npgreat"
    res = ""
    
    for i in range(len(df_pieces)):
        part_id = df_pieces.loc[i, "seq_id"]
        part_coord1 = df_pieces.loc[i, "coord1"]
        part_coord2 = df_pieces.loc[i, "coord2"]
        
        if(part_id == 'x'):
            continue
        
        record = SeqIO.read(folder_loc_orientation + "/" + part_id + ".fa", "fasta")
        sequences_string = record.seq
        sequences_id = record.id
        
        if(part_coord1 == 0 and part_coord2 == 0):
            part = sequences_string[:]
        elif(part_coord1 == 0 and part_coord2 != 0):
            part = sequences_string[:part_coord2-1]
        elif(part_coord2 == 0 and part_coord1 != 0):
            part = sequences_string[part_coord1-1:]
        else:
            part = sequences_string[part_coord1-1:part_coord2-1]
            
        res =  res + part
        
        seq_total = SeqRecord(
            Seq(res),
            id=assembly_id,
            #name="",
            description="",
        )

        SeqIO.write(seq_total, assembly_id + ".fasta", "fasta")
    
    return(res)


assembly_seq = functionCombine(df_combined)

len_assembly_seq = len(assembly_seq)
print("The NPGREAT Combination step finished.")
print("*******************************************************************************")
print("*** The NPGREAT assembly has been computed. Total assembly length: "+ str(len_assembly_seq)+"bp *** \nPlease see output file: assembly_npgreat.fasta")
print("*******************************************************************************")

