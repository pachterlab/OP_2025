import os
import re
import pandas as pd
import subprocess 

def read_fastq(fastq_file):
    if fastq_file.endswith(".gz"):
        file = gzip.open(fastq_file, "rt")  # 'rt' mode is for reading text
    else:
        file = open(fastq_file, "r")  # 'r' mode is for reading text

    with file:
        while True:
            header = file.readline().strip()
            sequence = file.readline().strip()
            plus_line = file.readline().strip()
            quality = file.readline().strip()

            if not header:
                break

            yield header, sequence, plus_line, quality

def get_header_set_from_fastq(fastq_file, output_format = "set"):
    if output_format == "set":
        headers = {header[1:].strip() for header, _, _, _ in read_fastq(fastq_file)}
    elif output_format == "list":
        headers = [header[1:].strip() for header, _, _, _ in read_fastq(fastq_file)]
    return headers

def intersect_lists(series):
    return list(set.intersection(*map(set, series)))

def map_transcripts_to_genes(transcript_list, mapping_dict):
    return [mapping_dict.get(transcript, "Unknown") for transcript in transcript_list]

#* TODO: only works when kb count was run with --num (as this means that each row of the BUS file corresponds to exactly one read)
def make_bus_df(kallisto_out, fastq_file, t2g_file, mm = False, union = False, assay = "bulk", bustools = "bustools"):
    with open(f"{kallisto_out}/transcripts.txt") as f:
        transcripts = f.read().splitlines()  # get transcript at index 0 with transcript[0], and index of transcript named "name" with transcript.index("name")

    transcripts.append("dlist")  # add dlist to the end of the list
    
    # important for temp files
    fastq_file = str(fastq_file)

    # TODO: figure out how to handle multiple fastqs - i.e., how to extract barcode from read, and how to handle the case where the barcode is associated with the file as a whole
    if fastq_file.endswith(".fastq") or fastq_file.endswith(".fq") or fastq_file.endswith(".fastq.gz") or fastq_file.endswith(".fq.gz"):
        fastq_header_list = get_header_set_from_fastq(fastq_file, output_format="list")
    elif fastq_file.endswith(".txt"):
        with open(fastq_file) as f:
            fastq_header_list = f.read().splitlines()
    
    # fastq_header_list.insert(0, "")  # because bus file read indices are index 1

    # Get equivalence class that matches to 0-indexed line number of target ID
    ec_df = pd.read_csv(f"{kallisto_out}/matrix.ec", sep="\t", header=None, names=["EC", "transcript_ids"])
    ec_df['transcript_ids'] = ec_df['transcript_ids'].astype(str)
    ec_df['transcript_ids_list'] = ec_df["transcript_ids"].str.split(",")
    ec_df["transcript_ids_list"] = ec_df["transcript_ids_list"].apply(lambda x: list(map(int, x)))
    ec_df["transcript_ids_list"] = ec_df["transcript_ids"].apply(lambda x: list(map(int, x.split(','))))
    ec_df["transcript_names"] = ec_df["transcript_ids_list"].apply(lambda ids: [transcripts[i] for i in ids])
    
    # Get bus output (converted to txt)
    bus_file = f"{kallisto_out}/output.bus"
    bus_text_file = f"{kallisto_out}/output_sorted_bus.txt"
    if not os.path.exists(bus_text_file):
        bus_txt_file_existed_originally = False
        create_bus_txt_file_command = f"{bustools} text -o {bus_text_file} -f {bus_file}"
        subprocess.run(create_bus_txt_file_command, shell=True, check=True)
        # /home/jrich/miniconda3/envs/cartf/lib/python3.10/site-packages/kb_python/bins/linux/bustools/bustools text -p -a -f -d /home/jrich/Desktop/CART_prostate_sc/TEMP_dlist_tests/kb_count_out_delaney/output.bus
    else:
        bus_txt_file_existed_originally = True
    
    bus_df = pd.read_csv(bus_text_file, sep="\t", header=None, names=['barcode', 'UMI', 'EC', 'count', 'read_index'])
    
    if not bus_txt_file_existed_originally:
        os.remove(bus_text_file)

    # TODO: if multiple fastqs, then consider the barcode as well (requires technology parameter)
    bus_df["fastq_header"] = bus_df["read_index"].map(pd.Series(fastq_header_list))

    bus_df = bus_df.merge(ec_df, on="EC", how="left")

    if assay == "sc":
        bus_df_collapsed_1 = bus_df.groupby(["barcode", "UMI", "EC"], as_index=False).agg({
            "count": "sum",  # Sum counts
            "read_index": lambda x: list(x),  # Combine ints in a list
            "fastq_header": lambda x: list(x),  # Combine strings in a list
            "transcript_ids": "first",  # Take the first value for all other columns
            "transcript_ids_list": "first",  # Take the first value for all other columns
            "transcript_names": "first"  # Take the first value for all other columns
        })

        bus_df_collapsed_2 = bus_df_collapsed_1.groupby(["barcode", "UMI"], as_index=False).agg({
            "EC": lambda x: list(x),
            "count": "sum",                                 # Sum the 'count' column
            "read_index": lambda x: sum(x, []),             # Concatenate lists in 'read_index'
            "fastq_header": lambda x: sum(x, []),            # Concatenate lists in 'fastq_header'
            "transcript_ids": lambda x: ",".join(x),    # Join strings in 'transcript_ids_list' with commas  # may contain duplicates indices
            "transcript_ids_list": lambda x: sum(x, []),       # Concatenate lists for 'transcript_ids_list'
            "transcript_names": lambda x: sum(x, []),          # Concatenate lists for 'transcript_names'
        })

        # Add new columns for the intersected lists
        bus_df_collapsed_2["transcript_names_final"] = bus_df_collapsed_1.groupby(["barcode", "UMI"])["transcript_names"].apply(intersect_lists).values
        bus_df_collapsed_2["transcript_ids_list_final"] = bus_df_collapsed_1.groupby(["barcode", "UMI"])["transcript_ids_list"].apply(intersect_lists).values

        bus_df = bus_df_collapsed_2

    elif assay == "bulk":
        # bus_df.rename(columns={"transcript_ids_list": "transcript_ids_list_final", "transcript_names": "transcript_names_final"}, inplace=True)
        bus_df["transcript_ids_list_final"] = bus_df["transcript_ids_list"]
        bus_df["transcript_names_final"] = bus_df["transcript_names"]

    t2g_df = pd.read_csv(t2g_file, sep="\t", header=None, names=["transcript_id", "gene_name"])
    t2g_dict = dict(zip(t2g_df["transcript_id"], t2g_df["gene_name"]))

    # Apply the mapping function to create gene name columns
    bus_df["gene_names"] = bus_df["transcript_names"].apply(lambda x: map_transcripts_to_genes(x, t2g_dict))
    bus_df["gene_names_final"] = bus_df["transcript_names_final"].apply(lambda x: map_transcripts_to_genes(x, t2g_dict))

    bus_df["gene_names_final_set"] = bus_df["gene_names_final"].apply(set)

    if union or mm:
        # union or mm gets added to count matrix as long as dlist is not included in the EC
        bus_df["counted_in_count_matrix"] = bus_df["transcript_names_final"].apply(lambda x: 'dlist' not in x)
    else:
        # only gets added to the count matrix if EC has exactly 1 gene
        bus_df["counted_in_count_matrix"] = bus_df["gene_names_final_set"].apply(lambda x: len(x) == 1)
    
    # adata_path = f"{kallisto_out}/counts_unfiltered/adata.h5ad"
    # adata = sc.read_h5ad(adata_path)
    # barcode_length = len(adata.obs.index[0])
    # bus_df['barcode_without_padding'] = bus_df['barcode'].str[(32 - barcode_length):]

    # so now I can iterate through this dataframe for the columns where counted_in_count_matrix is True - barcode will be the cell/sample (adata row), gene_names_final will be the list of gene name(s) (adata column), and count will be the number added to this entry of the matrix (always 1 for bulk)

    # save bus_df
    bus_df.to_csv(f"{kallisto_out}/bus_df.csv", index=False)
    return bus_df

rootdir = '/home/oakesc/virus_rat'
regex = re.compile('SRR(?P<SRR>\d*)_fastq')
for root, dirs, files in os.walk(rootdir):
	for dir_ in dirs:
		if regex.match(dir_):
			SRR = 'SRR' + re.match(regex, dir_).group(1)
			print(SRR)
			temp_kb_count_out_folder = '/home/oakesc/virus_rat/count_matrices_'+SRR
			temp_fastq_file = '/home/oakesc/virus_rat/'+SRR+'_fastq/'+SRR+'_2.fastq'
			temp_t2g_file = '/home/oakesc/references/virus/palmdb_clustered_t2g.txt'


			bus_df = make_bus_df(kallisto_out = temp_kb_count_out_folder, fastq_file = temp_fastq_file, t2g_file = temp_t2g_file, mm = False, union = False, assay = "bulk", bustools = "/home/oakesc/.local/lib/python3.9/site-packages/kb_python/bins/linux/bustools/bustools")
#read_to_ref_dict = dict(zip(bus_df['fastq_header'], bus_df['gene_names_final']))
