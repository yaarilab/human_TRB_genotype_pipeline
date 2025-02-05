$HOSTNAME = ""
params.outdir = 'results'  

//* params.nproc =  10  //* @input @description:"number of processes cores to use"
params.projectDir="${projectDir}"

params.IgBlastn.num_threads = params.nproc
params.IgBlastn.ig_seqtype = "TCR"
params.IgBlastn.outfmt = "MakeDb"
params.IgBlastn.num_alignments_V = "10"
params.IgBlastn.num_alignments_D = "3"
params.IgBlastn.num_alignments_J = "3"
params.IgBlastn.domain_system = "imgt"
params.IgBlastn.D_penalty = -1

params.MakeDb.failed = "true"
params.MakeDb.format = "airr"
params.MakeDb.regions = "default"
params.MakeDb.extended = "true"
params.MakeDb.asisid = "false"
params.MakeDb.asiscalls = "false"
params.MakeDb.inferjunction = "false"
params.MakeDb.partial = "false"
params.MakeDb.name_alignment = "first"

params.IgBlastn_genotype.num_threads = params.nproc
params.IgBlastn_genotype.ig_seqtype = "TCR"
params.IgBlastn_genotype.outfmt = "MakeDb"
params.IgBlastn_genotype.num_alignments_V = "10"
params.IgBlastn_genotype.num_alignments_D = "3"
params.IgBlastn_genotype.num_alignments_J = "3"
params.IgBlastn_genotype.domain_system = "imgt"
params.IgBlastn_genotype.D_penalty = -1

params.MakeDb_genotype.failed = "true"
params.MakeDb_genotype.format = "airr"
params.MakeDb_genotype.regions = "default"
params.MakeDb_genotype.extended = "true"
params.MakeDb_genotype.asisid = "false"
params.MakeDb_genotype.asiscalls = "false"
params.MakeDb_genotype.inferjunction = "false"
params.MakeDb_genotype.partial = "false"
params.MakeDb_genotype.name_alignment = "Final"

params.trb_deletion.gene_usages_file = "${params.projectDir}/trbv_usage.tsv"

params.ogrdbstats_report.chain = "TRBV"

if (!params.airr_seq){params.airr_seq = ""} 
if (!params.v_germline_file){params.v_germline_file = ""} 
if (!params.d_germline){params.d_germline = ""} 
if (!params.j_germline){params.j_germline = ""} 

Channel.fromPath(params.airr_seq, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_0_fastaFile_g_7;g_0_fastaFile_g_10;g_0_fastaFile_g_8}
Channel.fromPath(params.v_germline_file, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_1_germlineFastaFile_g_4;g_1_germlineFastaFile_g_10;g_1_germlineFastaFile_g_11;g_1_germlineFastaFile_g_8;g_1_germlineFastaFile_g_25}
Channel.fromPath(params.d_germline, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_2_germlineFastaFile_g_5;g_2_germlineFastaFile_g_10;g_2_germlineFastaFile_g_25;g_2_germlineFastaFile_g_37;g_2_germlineFastaFile_g_8}
Channel.fromPath(params.j_germline, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_3_germlineFastaFile_g_6;g_3_germlineFastaFile_g_10;g_3_germlineFastaFile_g_25;g_3_germlineFastaFile_g_8}


process V_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_1_germlineFastaFile_g_4

output:
 file "${db_name}"  into g_4_germlineDb0_g_7

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process D_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_2_germlineFastaFile_g_5

output:
 file "${db_name}"  into g_5_germlineDb0_g_7

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process J_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_3_germlineFastaFile_g_6

output:
 file "${db_name}"  into g_6_germlineDb0_g_7

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process IgBlastn {

input:
 set val(name),file(fastaFile) from g_0_fastaFile_g_7
 file db_v from g_4_germlineDb0_g_7
 file db_d from g_5_germlineDb0_g_7
 file db_j from g_6_germlineDb0_g_7

output:
 set val(name), file("${outfile}") optional true  into g_7_igblastOut0_g_10, g_7_igblastOut0_g_8

script:
num_threads = params.IgBlastn.num_threads
ig_seqtype = params.IgBlastn.ig_seqtype
outfmt = params.IgBlastn.outfmt
num_alignments_V = params.IgBlastn.num_alignments_V
num_alignments_D = params.IgBlastn.num_alignments_D
num_alignments_J = params.IgBlastn.num_alignments_J
domain_system = params.IgBlastn.domain_system
auxiliary_data = params.IgBlastn.auxiliary_data
D_penalty = params.IgBlastn.D_penalty

randomString = org.apache.commons.lang.RandomStringUtils.random(9, true, true)
outname = name + "_" + randomString
outfile = (outfmt=="MakeDb") ? name+"_"+randomString+".out" : name+"_"+randomString+".tsv"
outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : "19"

if(db_v.toString()!="" && db_d.toString()!="" && db_j.toString()!=""){
	"""
	igblastn -query ${fastaFile} \
		-germline_db_V ${db_v}/${db_v} \
		-germline_db_D ${db_d}/${db_d} \
		-germline_db_J ${db_j}/${db_j} \
		-num_alignments_V ${num_alignments_V} \
		-num_alignments_D ${num_alignments_D} \
		-num_alignments_J ${num_alignments_J} \
		-D_penalty ${D_penalty} \
		-domain_system ${domain_system} \
		-ig_seqtype ${ig_seqtype} \
		-auxiliary_data ${auxiliary_data} \
		-outfmt ${outfmt} \
		-num_threads ${num_threads} \
		-out ${outfile}
	"""
}else{
	"""
	"""
}

}


process MakeDb {

input:
 set val(name),file(fastaFile) from g_0_fastaFile_g_8
 set val(name_igblast),file(igblastOut) from g_7_igblastOut0_g_8
 set val(name1), file(v_germline_file) from g_1_germlineFastaFile_g_8
 set val(name2), file(d_germline_file) from g_2_germlineFastaFile_g_8
 set val(name3), file(j_germline_file) from g_3_germlineFastaFile_g_8

output:
 set val(name_igblast),file("*_db-pass.tsv") optional true  into g_8_outputFileTSV0_g_11
 set val("reference_set"), file("${reference_set}") optional true  into g_8_germlineFastaFile11
 set val(name_igblast),file("*_db-fail.tsv") optional true  into g_8_outputFileTSV22

script:

failed = params.MakeDb.failed
format = params.MakeDb.format
regions = params.MakeDb.regions
extended = params.MakeDb.extended
asisid = params.MakeDb.asisid
asiscalls = params.MakeDb.asiscalls
inferjunction = params.MakeDb.inferjunction
partial = params.MakeDb.partial
name_alignment = params.MakeDb.name_alignment

failed = (failed=="true") ? "--failed" : ""
format = (format=="changeo") ? "--format changeo" : ""
extended = (extended=="true") ? "--extended" : ""
regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
asisid = (asisid=="true") ? "--asis-id" : ""
asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
partial = (partial=="true") ? "--partial" : ""

reference_set = "reference_set_makedb_"+name_alignment+".fasta"

outname = name_igblast+'_'+name_alignment
println name_alignment

if(igblastOut.getName().endsWith(".out")){
	"""
	
	cat ${v_germline_file} ${d_germline_file} ${j_germline_file} > ${reference_set}
	
	MakeDb.py igblast \
		-s ${fastaFile} \
		-i ${igblastOut} \
		-r ${v_germline_file} ${d_germline_file} ${j_germline_file} \
		--log MD_${name}.log \
		--outname ${outname}\
		${extended} \
		${failed} \
		${format} \
		${regions} \
		${asisid} \
		${asiscalls} \
		${inferjunction} \
		${partial}
	"""
}else{
	"""
	
	"""
}

}


process igdiscover_igblast {

input:
 set val(name_igblast),file(igblastOut) from g_7_igblastOut0_g_10
 set val(name1), file(v_germline_file) from g_1_germlineFastaFile_g_10
 set val(name1), file(d_germline_file) from g_2_germlineFastaFile_g_10
 set val(name1), file(j_germline_file) from g_3_germlineFastaFile_g_10
 set val(name),file(fastaFile) from g_0_fastaFile_g_10

output:
 set val(name_igblast),file("*_igdiscover-pass.tsv") optional true  into g_10_outputFileTSV0_g_11

script:

"""

mkdir database

# move the germline files to the database folder
awk '/^>/ {print; next} {gsub(/[.-]/, ""); print}' ${v_germline_file} > database/V.fasta
cp ${j_germline_file} database/J.fasta
cp ${d_germline_file} database/D.fasta

igdiscover igblast --threads 1 --no-cache --pathig ${igblastOut} database ${fastaFile} > ${name_igblast}_igdiscover-pass.tsv

"""
}


process tcr_data_to_igdiscover {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_makedb-pass_mut.tsv$/) "reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_igdiscover-pass_mut.tsv$/) "reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapsed.fasta$/) "reads/$filename"}
input:
 set val(name),file(makedb) from g_8_outputFileTSV0_g_11
 set val(name1), file(v_germline_file) from g_1_germlineFastaFile_g_11
 set val(name2),file(igdiscover) from g_10_outputFileTSV0_g_11

output:
 set val(name),file("*_makedb-pass_mut.tsv") optional true  into g_11_outputFileTSV0_g_25
 set val(name),file("*_igdiscover-pass_mut.tsv") optional true  into g_11_outputFileTSV11
 set val(name),file("*_collapsed.fasta") optional true  into g_11_germlineFastaFile2_g_29, g_11_germlineFastaFile2_g_30

script:
"""
#!/usr/bin/env Rscript
library(piglet)
library(stringr)

data <- read.delim("${makedb}", header=T, sep = "\t", stringsAsFactors = F)
germline <- tigger::readIgFasta("${v_germline_file}")

# get the v start and sequence
data[["v_start"]] <- stringi::stri_locate(data[["sequence_alignment"]],regex = "[ATCG]")
data[["v_seq"]] <- sapply(1:nrow(data),function(i) substr(data[["sequence_alignment"]][i],1,data[["v_germline_end"]][i]))
# get the mutation count in region for each sequence
data[["v_mut"]] <- sapply(1:nrow(data), function(i){
  allele <- data[["v_call"]][i]
  # get the first call
  allele <- strsplit(allele, ",", fixed = T)[[1]][1]
  # get the mutated positions
  idx <- piglet::allele_diff(germline[[allele]], data[["v_seq"]][i])
  # find the minimal v start position + 5, and only consider mutation above it
  v_min <- min(data[["v_start"]][grep(allele, data[["v_call"]],fixed=T)])+5
  # sum the number of mutation and check if below or equal to 3.
  sum(idx>v_min & idx<=316)<=3;
})
# filter the data
data_filter <- data[data[["v_mut"]],]

# read igdiscover data
data_igdiscover <- read.delim("${igdiscover}", header=T, sep = "\t", stringsAsFactors = F)
# sort the names
data_igdiscover[["sequence_id"]] <- sapply(data_igdiscover[["name"]], function(x) strsplit(x,"|",fixed=T)[[1]][1])
data_igdiscover <- data_igdiscover[data_igdiscover[["sequence_id"]] %in% data[["sequence_id"]],]

# add cdr3 information
ig_ind <- 1
for (seq_id in data_igdiscover[["sequence_id"]]) {
  seq_ind <- which(data[["sequence_id"]]==seq_id)
  cdr3 <- data[["junction"]][[seq_ind]]
  if (str_length(cdr3) > 6) {
    cdr3 <- substr(cdr3, 4, str_length(cdr3)-3)
    cdr3_aa <- data[["junction_aa"]][[seq_ind]]
    cdr3_aa <- substr(cdr3_aa, 2, str_length(cdr3_aa)-1)
  } else {
    cdr3 <- ""
    cdr3_aa <- ""
  }
  data_igdiscover[["CDR3_nt"]][ig_ind] <- cdr3
  data_igdiscover[["CDR3_aa"]][ig_ind] <- cdr3_aa
  ig_ind <- ig_ind + 1
}

## write both tables

write.table(data, paste0("${name}", "_makedb-pass_mut.tsv"), sep = "\t")
write.table(data_igdiscover, paste0("${name}", "_igdiscover-pass_mut.tsv"), sep = "\t")

# collapse identical sequences
library(dplyr)

if (!"consensus_count" %in% names(data)) {
  if ("reads" %in% names(data)) {
    if ("templates" %in% names(data)) {
      data[["templates"]][!unlist(lapply(data[["templates"]], is.integer))] <- 1
      data[["reads"]][!unlist(lapply(data[["reads"]], is.integer))] <- data[["templates"]][!unlist(lapply(data[["reads"]], is.integer))]
    } else {
      data[["reads"]][!unlist(lapply(data[["reads"]], is.integer))] <- 1
      data[["templates"]] <- 1
    }
    data[["consensus_count"]] <- data[["reads"]]
    data[["duplicate_count"]] <- data[["templates"]]
  } else {
    data[["consensus_count"]] <- 1
    data[["duplicate_count"]] <- 1
  }
} else {
  data[["duplicate_count"]][!unlist(lapply(data[["duplicate_count"]], is.integer))] <- 1
  data[["consensus_count"]][!unlist(lapply(data[["consensus_count"]], is.integer))] <- data[["duplicate_count"]][!unlist(lapply(data[["consensus_count"]], is.integer))]
}

data <- data %>% select(sequence, sequence_alignment, sequence_id,  consensus_count, duplicate_count)
data[["sequence_alignment"]] <- gsub(".", "", data[["sequence_alignment"]], fixed = TRUE)

vdj_seqs <- unique(data[["sequence_alignment"]])
for (vdj_seq in vdj_seqs) {
  vdj_indexes <- which(data[["sequence_alignment"]] == vdj_seq)
  if (length(vdj_indexes) > 1) {
    temp <- data[vdj_indexes, ]
    temp[["consensus_count"]] <- as.numeric(temp[["consensus_count"]])
    temp[["duplicate_count"]] <- as.numeric(temp[["duplicate_count"]])
    
    data <- data[-vdj_indexes, ]
    
    seq_id_ind <- which.max(temp[["consensus_count"]])
    seq_id <- temp[["sequence_id"]][seq_id_ind]
    seq_id_ind <- vdj_indexes[seq_id_ind]
    sequence <- temp[["sequence"]][seq_id_ind]
    
    conscount <- sum(temp[["consensus_count"]])
    dupcount <- sum(temp[["duplicate_count"]])
    data[nrow(data) + 1, ] <- c(sequence, vdj_seq, seq_id, conscount, dupcount)
  }
}

seq.names <- sapply(1:nrow(data), function(x){ 
  paste0(names(data)[3:ncol(data)], rep('=', length(3:ncol(data))), data[x, 3:ncol(data)], collapse = '|')
})
seq.names <- gsub('sequence_id=', '', seq.names, fixed = TRUE)

tigger::writeFasta(setNames(as.list(data[["sequence"]]), seq.names), file = paste0("${name}","_collapsed.fasta"))

"""
}


process trb_genotype_inference {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*v_genotype.tsv$/) "genotype/$filename"}
input:
 set val(name), file(airrseq) from g_11_outputFileTSV0_g_25
 set val(namev), file(germline_v) from g_1_germlineFastaFile_g_25
 set val(named), file(germline_d) from g_2_germlineFastaFile_g_25
 set val(namej), file(germline_j) from g_3_germlineFastaFile_g_25

output:
 set val(name),file("*v_genotype.tsv")  into g_25_outputFileTSV0_g_34
 set val(namev),file("V_gapped_personal.fasta")  into g_25_germlineFastaFile1_g_26, g_25_germlineFastaFile1_g_40, g_25_germlineFastaFile1_g_30
 set val(named),file("D_personal.fasta")  into g_25_germlineFastaFile2_g_27, g_25_germlineFastaFile2_g_30
 set val(namej),file("J_personal.fasta")  into g_25_germlineFastaFile3_g_28, g_25_germlineFastaFile3_g_30
 set val(name),file("*d_genotype.tsv")  into g_25_outputFileTSV4_g_35
 set val(name),file("*j_genotype.tsv")  into g_25_outputFileTSV5_g_35

script:

min_consensus_count = params.trb_genotype_inference.min_consensus_count
filter_chimera = params.trb_genotype_inference.filter_chimera

"""
#!/usr/bin/env Rscript

library(tigger)
library(dplyr)
library(stringi)
library(stringr)

# max V position to look for novel SNPs 
max_snp_position <- 316

# undocumented_alleles_2_ignore <- c()
undocumented_alleles_2_ignore <- c("TRBV13*01_A170T", "TRBV13*01_T158C", "TRBV10-3*02_C225G", "TRBV20-1*01_C142A", "TRBV30*01_A113C", "TRBV6-6*01_C261T",
                                   "TRBV7-9*05_A19G_C256T",
                                   "TRBV15*bp02_A316C", "TRBV5-4*bp01_C159T", "TRBV6-6*bp03_G216C", "TRBV6-6*bp03_T201C_A202C_G216C", "TRBV6-6*bp03_T231C_C261T",
                                   "TRBV15*bp02_G153T", "TRBV19*bp01_T310C_G311C_C314T", "TRBV5-4*bp01_G205A", "TRBV5-5*bp01_G232A", 
                                   "TRBV7-9*bp04_T312A", "TRBV6-6*bp01_C261T", "TRBV10-2*bp01_C214T", 
                                   "TRBV30*bp01_T316G", "TRBV19*bp01_G313T_C315T_A316C", "TRBV5-6*bp01_C223G", "TRBV6-1*bp01_C278A", "TRBV15*bp02_C147A", 
                                   "TRBV18*bp01_G289A", "TRBV11-1*bp01_A164G", "TRBV30*bp01_C168A", "TRBV10-2*bp02_C154A", "TRBV10-1*bp01_C284G", 
                                   "TRBV7-9*bp01_T312C", "TRBV19*bp01_C293T_G294A", "TRBV15*bp02_G275A", "TRBV27*bp01_A155C", "TRBV30*bp01_G169A",
                                   "TRBV5-6*bp01_G233C_A236G", "TRBV11-2*bp01_A238T", "TRBV10-1*bp02_G156A_G274T", "TRBV24-1*bp01_A316C",
                                   "TRBV10-1*bp01_C190A_C195T_A199G", "TRBV5-5*bp01_T284G_G303C", "TRBV6-9*bp01_G155T_C156G_A303G", "TRBV7-4*bp01_T306C_C307T",
                                   "TRBV10-1*bp01_G274T", 
                                   "TRBV20-1*ap02_T310G", "TRBV7-8*ap01_T295C", "TRBV7-4*ap01_G291C_A297G", "TRBV7-4*ap01_G291C_A297G_C314T", "TRBV7-9*ap01_G313T",
                                   "TRBV4-3*ap01_A305C_T306C", "TRBV4-3*ap01_A305C_T306C_T308C", "TRBV4-3*ap01_G311C_G313C", "TRBV4-3*ap01_T308C_G311C", "TRBV4-3*ap01_T308C_T310C_G311C")

filter_chimera_bool <- as.logical("${filter_chimera}")

TRBV_GERM <- readIgFasta("${germline_v}")
TRBD_GERM <- readIgFasta("${germline_d}")
TRBJ_GERM <- readIgFasta("${germline_j}")

DATA <- read.delim("${airrseq}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# filter by the selected minimal constcount.
DATA <- DATA[DATA[["consensus_count"]] >= ${min_consensus_count}, ]

# filter by zero mutations over the V
DATA[["v_seq"]] <- substr(DATA[["sequence_alignment"]], 1, sapply(DATA[["v_germline_end"]], min, max_snp_position)) 
DATA[["v_mut"]] <- sapply(tigger::getMutCount(DATA[["v_seq"]], DATA[["v_call"]], germline_db = TRBV_GERM), function(x){x[[1]]}) 
DATA <- DATA[DATA[["v_mut"]] <= 1, ]

DATA_V_SA <- DATA[!grepl(pattern = ',', DATA[["v_call"]]), ]
DATA_V_SA <- DATA_V_SA[!DATA_V_SA[["v_call"]] %in% undocumented_alleles_2_ignore, ]

geno_BV <- inferGenotypeBayesian(DATA_V_SA, germline_db = TRBV_GERM, find_unmutated = FALSE, novel = new_novel_df_H, v_call = 'v_call')
names(geno_BV) <- names(geno_BV)
geno_BV[["genotyped_alleles"]] <- apply(geno_BV[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})

if (filter_chimera_bool) {
  trizygous <- geno_BV[str_count(geno_BV[["genotyped_alleles"]], ",") >= 2, ]
  trizygous <- trizygous %>% tidyr::separate_rows(genotyped_alleles, sep = ",")
  trizygous[["full_name"]] <- paste0(trizygous[["gene"]], "*", trizygous[["genotyped_alleles"]])
  novel_names <- trizygous[["full_name"]][grepl("_[A-Z]", trizygous[["genotyped_alleles"]])]
  
  if (length(novel_names)){
    novel_allele_mismatches <- list()
    for (novel in novel_names) {
      novel_seq <- unlist(strsplit(as.character(TRBV_GERM[[novel]]), ""))
      gene <- sapply(strsplit(novel, "*", fixed = TRUE), "[", 1)
      gene_alleles <- trizygous[["full_name"]][trizygous[["gene"]] == gene]
      gene_alleles <- TRBV_GERM[gene_alleles]
      
      for (allele in names(gene_alleles)) {
        if (allele == novel) {next}
        allele_seq <- unlist(strsplit(as.character(gene_alleles[allele]), ""))
        mismatches_counter <- c(0)
        for (pos in 2:max_snp_position) {
          mismatches_counter[pos] <- mismatches_counter[pos - 1]
          if (pos > min(str_length(novel_seq), str_length(allele_seq))) {
            next
          }
          if ((novel_seq[[pos]] != ".") & (allele_seq[[pos]] != ".") & (novel_seq[[pos]] != allele_seq[[pos]])) {
            mismatches_counter[pos] <- mismatches_counter[pos] + 1
          }
        }
        novel_allele_mismatches[[novel]][[allele]] <- list()
        novel_allele_mismatches[[novel]][[allele]][["prefix"]] <- mismatches_counter
        novel_allele_mismatches[[novel]][[allele]][["suffix"]] <- mismatches_counter[max_snp_position] - mismatches_counter
      }
    }
    
    chimera_alleles <- c()
    prefix_alleles <- c()
    suffix_alleles <- c()
    min_mismatches <- c()
    position <- c()
    for (novel_allele in names(novel_allele_mismatches)) {
      for (prefix_allele in names(novel_allele_mismatches[[novel_allele]])) {
        for (suffix_allele in names(novel_allele_mismatches[[novel_allele]])) {
          if (suffix_allele != prefix_allele) {
            chimera_alleles <- c(chimera_alleles, novel_allele)
            prefix_alleles <- c(prefix_alleles, prefix_allele)
            suffix_alleles <- c(suffix_alleles, suffix_allele)
            mismatches <- novel_allele_mismatches[[novel_allele]][[prefix_allele]][["prefix"]] + novel_allele_mismatches[[novel_allele]][[suffix_allele]][["suffix"]]
            min_mismatch <- min(mismatches)
            min_mismatches <- c(min_mismatches, min_mismatch)
            position <- c(position, which.min(mismatches))
          }
        }
      }
    }
    
    chimera_df <- data.frame(chimera_alleles, prefix_alleles, suffix_alleles, min_mismatches, position, stringsAsFactors = FALSE)
    if (nrow(chimera_df)) {
      chimera_df[["snp_count"]] <- str_count(chimera_df[["chimera_alleles"]], "_")
      chimera_df <- chimera_df[chimera_df[["min_mismatches"]] < chimera_df[["snp_count"]], ]
    }
    
    if (nrow(chimera_df)) {
      pos_chimera_df <- do.call(rbind, unname(by(chimera_df, chimera_df[["chimera_alleles"]], function(x) x[x[["min_mismatches"]] == min(x[["min_mismatches"]]), ])))
      pos_chimera_df[["prefix_gene"]] <- sapply(strsplit(as.character(pos_chimera_df[["prefix_alleles"]]), "*", fixed = TRUE), '[', 1)
      pos_chimera_df[["suffix_gene"]] <- sapply(strsplit(as.character(pos_chimera_df[["suffix_alleles"]]), "*", fixed = TRUE), '[', 1)
      
      prob_chimera <- unique(pos_chimera_df[["chimera_alleles"]][pos_chimera_df[["min_mismatches"]] == 0])
    } 
    
    if (length(prob_chimera)) {
      DATA_V_SA <- DATA_V_SA[!DATA_V_SA[["v_call"]] %in% prob_chimera, ]
      geno_BV <- inferGenotypeBayesian(DATA_V_SA, germline_db = TRBV_GERM, find_unmutated = FALSE, novel = new_novel_df_H, v_call = 'v_call')
      geno_BV[["genotyped_alleles"]] <- apply(geno_BV[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:5]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})
    }
  }
}

DATA_D_geno <- DATA[(!grepl(pattern = ',', DATA[["d_call"]]) & DATA[["d_call"]] != 'None') & (DATA[["d_sequence_end"]] - DATA[["d_sequence_start"]] >= 8), ]
DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno[["sequence_id"]]), ]

# extract d sequence in the direct orientation
DATA_D_reg <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] < DATA_D_geno[["d_germline_end"]], ]
DATA_D_reg[["d_seq"]] <- substr(DATA_D_reg[["sequence"]], DATA_D_reg[["d_sequence_start"]], DATA_D_reg[["d_sequence_end"]])

# extract convert d sequence in the inverted orientation to the direct orientation
DATA_D_inv <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] > DATA_D_geno[["d_germline_end"]], ]
DATA_D_inv[["d_seq"]] <- substr(DATA_D_inv[["sequence"]], DATA_D_inv[["d_sequence_start"]], DATA_D_inv[["d_sequence_end"]])
DATA_D_inv[["d_seq"]] <- stringi::stri_reverse(DATA_D_inv[["d_seq"]])
DATA_D_inv[["d_seq"]] <- gsub("A", "t", DATA_D_inv[["d_seq"]])
DATA_D_inv[["d_seq"]] <- gsub("T", "a", DATA_D_inv[["d_seq"]])
DATA_D_inv[["d_seq"]] <- gsub("G", "c", DATA_D_inv[["d_seq"]])
DATA_D_inv[["d_seq"]] <- gsub("C", "g", DATA_D_inv[["d_seq"]])
DATA_D_inv[["d_seq"]] <- toupper(DATA_D_inv[["d_seq"]])

d_germ_end <- DATA_D_inv[["d_germline_start"]]
DATA_D_inv[["d_germline_start"]] <- DATA_D_inv[["d_germline_end"]]
DATA_D_inv[["d_germline_end"]] <- d_germ_end

DATA_D_geno <- rbind(DATA_D_reg, DATA_D_inv)

# filter by zero mutations over the D segment
DATA_D_geno[["mut_d"]] <- unlist(lapply(1:nrow(DATA_D_geno), function(i) {
  mut <- 0
  row_seq <- unlist(strsplit(DATA_D_geno[["d_seq"]][[i]], ""))
  allele_seq <- unlist(strsplit(TRBD_GERM[[DATA_D_geno[["d_call"]][[i]]]], ""))
  for (pos in DATA_D_geno[["d_germline_start"]][[i]]:DATA_D_geno[["d_germline_end"]][[i]]) {
    if (row_seq[pos - (DATA_D_geno[["d_germline_start"]][[i]] - 1)] != allele_seq[pos]) {
      mut <- mut + 1
    }
  }
  mut
}))

DATA_D_geno <- DATA_D_geno[DATA_D_geno[["mut_d"]] == 0, ]

geno_BD <- inferGenotypeBayesian(DATA_D_geno, find_unmutated = FALSE, germline_db = TRBD_GERM, v_call = 'd_call')
geno_BD[["genotyped_alleles"]] <- apply(geno_BD[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})

D2_total <- nrow(DATA_D_geno[grepl("TRBD2", DATA_D_geno[["d_call"]]), ])
D2_01_count <- nrow(DATA_D_geno[DATA_D_geno[["d_call"]] == "TRBD2*01", ])
D2_01_freq <- D2_01_count / D2_total

if (D2_01_freq < 0.2066) {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "02"
} else if (D2_01_freq > 0.8969) {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "01"
} else if (geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] == "01") {
  geno_BD[["genotyped_alleles"]][geno_BD[["gene"]] == "TRBD2"] <- "01,02"
}

DATA_J_SA <- DATA[!grepl(pattern = ',', DATA[["j_call"]]), ]
geno_BJ <- inferGenotypeBayesian(DATA, germline_db = TRBJ_GERM, find_unmutated = FALSE, v_call = 'j_call')
geno_BJ[["genotyped_alleles"]] <- apply(geno_BJ[, c(2, 6:9)], 1, function(y){m <- which.max(as.numeric(y[2:3]));paste0(unlist(strsplit((y[1]), ','))[1:m], collapse = ",")})


## Remove from TRBV_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBV_GERM),'*',fixed=T),'[',1) %in%  geno_BV[["gene"]])
TRBV_GERM.NEW <- TRBV_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BV)){
  gene <- geno_BV[["gene"]][i]
  
  alleles <- geno_BV[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBV_GERM) %in%  paste(gene,alleles,sep='*')
  TRBV_GERM.NEW <- c(TRBV_GERM.NEW,TRBV_GERM[IND])
}


## Remove from TRBD_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBD_GERM),'*',fixed=T),'[',1) %in%  geno_BD[["gene"]])
TRBD_GERM.NEW <- TRBD_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BD)){
  gene <- geno_BD[["gene"]][i]
  alleles <- geno_BD[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBD_GERM) %in%  paste(gene,alleles,sep='*')
  TRBD_GERM.NEW <- c(TRBD_GERM.NEW,TRBD_GERM[IND])
}

## Remove from TRBJ_GERM irrelevant alleles
NOTGENO.IND <- !(sapply(strsplit(names(TRBJ_GERM),'*',fixed=T),'[',1) %in%  geno_BJ[["gene"]])
TRBJ_GERM.NEW <- TRBJ_GERM[NOTGENO.IND]

for(i in 1:nrow(geno_BJ)){
  gene <- geno_BJ[["gene"]][i]
  alleles <- geno_BJ[["genotyped_alleles"]][i]
  alleles <- unlist(strsplit(alleles,','))
  IND <- names(TRBJ_GERM) %in%  paste(gene,alleles,sep='*')
  TRBJ_GERM.NEW <- c(TRBJ_GERM.NEW,TRBJ_GERM[IND])
}


### CHECK IF THE REPLACEMENT IS CORRECT

## Combine the genotyped and others and write to a fasta file for reference

writeFasta(TRBV_GERM.NEW, file = "V_gapped_personal.fasta")
writeFasta(TRBD_GERM.NEW, file = "D_personal.fasta")
writeFasta(TRBJ_GERM.NEW, file = "J_personal.fasta")

## save the genotype data
write.table(geno_BV, file = paste0("${name}","_v_genotype.tsv"), quote = F, row.names = F, sep = "\t")
write.table(geno_BD, file = paste0("${name}","_d_genotype.tsv"), quote = F, row.names = F, sep = "\t")
write.table(geno_BJ, file = paste0("${name}","_j_genotype.tsv"), quote = F, row.names = F, sep = "\t")
"""

}


process J_MakeBlastDb_genotype {

input:
 set val(db_name), file(germlineFile) from g_25_germlineFastaFile3_g_28

output:
 file "${db_name}"  into g_28_germlineDb0_g_29

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process D_MakeBlastDb_genotype {

input:
 set val(db_name), file(germlineFile) from g_25_germlineFastaFile2_g_27

output:
 file "${db_name}"  into g_27_germlineDb0_g_29

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process V_MakeBlastDb_genotype {

input:
 set val(db_name), file(germlineFile) from g_25_germlineFastaFile1_g_26

output:
 file "${db_name}"  into g_26_germlineDb0_g_29

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process IgBlastn_genotype {

input:
 set val(name),file(fastaFile) from g_11_germlineFastaFile2_g_29
 file db_v from g_26_germlineDb0_g_29
 file db_d from g_27_germlineDb0_g_29
 file db_j from g_28_germlineDb0_g_29

output:
 set val(name), file("${outfile}") optional true  into g_29_igblastOut0_g_30

script:
num_threads = params.IgBlastn_genotype.num_threads
ig_seqtype = params.IgBlastn_genotype.ig_seqtype
outfmt = params.IgBlastn_genotype.outfmt
num_alignments_V = params.IgBlastn_genotype.num_alignments_V
num_alignments_D = params.IgBlastn_genotype.num_alignments_D
num_alignments_J = params.IgBlastn_genotype.num_alignments_J
domain_system = params.IgBlastn_genotype.domain_system
auxiliary_data = params.IgBlastn_genotype.auxiliary_data
D_penalty = params.IgBlastn_genotype.D_penalty

randomString = org.apache.commons.lang.RandomStringUtils.random(9, true, true)
outname = name + "_" + randomString
outfile = (outfmt=="MakeDb") ? name+"_"+randomString+".out" : name+"_"+randomString+".tsv"
outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : "19"

if(db_v.toString()!="" && db_d.toString()!="" && db_j.toString()!=""){
	"""
	igblastn -query ${fastaFile} \
		-germline_db_V ${db_v}/${db_v} \
		-germline_db_D ${db_d}/${db_d} \
		-germline_db_J ${db_j}/${db_j} \
		-num_alignments_V ${num_alignments_V} \
		-num_alignments_D ${num_alignments_D} \
		-num_alignments_J ${num_alignments_J} \
		-D_penalty ${D_penalty} \
		-domain_system ${domain_system} \
		-ig_seqtype ${ig_seqtype} \
		-auxiliary_data ${auxiliary_data} \
		-outfmt ${outfmt} \
		-num_threads ${num_threads} \
		-out ${outfile}
	"""
}else{
	"""
	"""
}

}


process MakeDb_genotype {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-pass.tsv$/) "rearrangements/$filename"}
input:
 set val(name),file(fastaFile) from g_11_germlineFastaFile2_g_30
 set val(name_igblast),file(igblastOut) from g_29_igblastOut0_g_30
 set val(name1), file(v_germline_file) from g_25_germlineFastaFile1_g_30
 set val(name2), file(d_germline_file) from g_25_germlineFastaFile2_g_30
 set val(name3), file(j_germline_file) from g_25_germlineFastaFile3_g_30

output:
 set val(name_igblast),file("*_db-pass.tsv") optional true  into g_30_outputFileTSV0_g_34, g_30_outputFileTSV0_g_35, g_30_outputFileTSV0_g_37, g_30_outputFileTSV0_g_40
 set val("reference_set"), file("${reference_set}") optional true  into g_30_germlineFastaFile1_g_40
 set val(name_igblast),file("*_db-fail.tsv") optional true  into g_30_outputFileTSV22

script:

failed = params.MakeDb_genotype.failed
format = params.MakeDb_genotype.format
regions = params.MakeDb_genotype.regions
extended = params.MakeDb_genotype.extended
asisid = params.MakeDb_genotype.asisid
asiscalls = params.MakeDb_genotype.asiscalls
inferjunction = params.MakeDb_genotype.inferjunction
partial = params.MakeDb_genotype.partial
name_alignment = params.MakeDb_genotype.name_alignment

failed = (failed=="true") ? "--failed" : ""
format = (format=="changeo") ? "--format changeo" : ""
extended = (extended=="true") ? "--extended" : ""
regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
asisid = (asisid=="true") ? "--asis-id" : ""
asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
partial = (partial=="true") ? "--partial" : ""

reference_set = "reference_set_makedb_"+name_alignment+".fasta"

outname = name_igblast+'_'+name_alignment
println name_alignment

if(igblastOut.getName().endsWith(".out")){
	"""
	
	cat ${v_germline_file} ${d_germline_file} ${j_germline_file} > ${reference_set}
	
	MakeDb.py igblast \
		-s ${fastaFile} \
		-i ${igblastOut} \
		-r ${v_germline_file} ${d_germline_file} ${j_germline_file} \
		--log MD_${name}.log \
		--outname ${outname}\
		${extended} \
		${failed} \
		${format} \
		${regions} \
		${asisid} \
		${asiscalls} \
		${inferjunction} \
		${partial}
	"""
}else{
	"""
	
	"""
}

}


process ogrdbstats_report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "ogrdbstats/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*csv$/) "ogrdbstats/$filename"}
input:
 set val(name),file(airrFile) from g_30_outputFileTSV0_g_40
 set val(name1), file(germline_file) from g_30_germlineFastaFile1_g_40
 set val(name2), file(v_germline_file) from g_25_germlineFastaFile1_g_40

output:
 file "*pdf"  into g_40_outputFilePdf00
 file "*csv"  into g_40_outputFileCSV11

script:

// general params
chain = params.ogrdbstats_report.chain
outname = airrFile.name.toString().substring(0, airrFile.name.toString().indexOf("_db-pass"))

"""

germline_file_path=\$(realpath ${germline_file})

novel=""

if grep -q "_[A-Z][0-9]" ${v_germline_file}; then
	awk '/^>/{f=0} \$0 ~ /_[A-Z][0-9]/ {f=1} f' ${v_germline_file} > novel_sequences.fasta
	novel=\$(realpath novel_sequences.fasta)
	diff \$germline_file_path \$novel | grep '^<' | sed 's/^< //' > personal_germline.fasta
	germline_file_path=\$(realpath personal_germline.fasta)
	novel="--inf_file \$novel"
fi

IFS='\t' read -a var < ${airrFile}

airrfile=${airrFile}

if [[ ! "\${var[*]}" =~ "v_call_genotyped" ]]; then
    awk -F'\t' '{col=\$5;gsub("call", "call_genotyped", col); print \$0 "\t" col}' ${airrFile} > ${outname}_genotyped.tsv
    airrfile=${outname}_genotyped.tsv
fi

airrFile_path=\$(realpath \$airrfile)


run_ogrdbstats \
	\$germline_file_path \
	"Homosapiens" \
	\$airrFile_path \
	${chain} \
	\$novel 

"""

}


process trb_deletion {

input:
 set val(name),file(airrseq) from g_30_outputFileTSV0_g_34
 set val(name1), file(genotype_v) from g_25_outputFileTSV0_g_34

output:
 set val(name), file("*_v_genotype_deletion.tsv")  into g_34_outputFileTSV0_g_35

script:
	
gene_usages_file = params.trb_deletion.gene_usages_file
min_consensus_count = params.trb_deletion.min_consensus_count

"""
#!/usr/bin/env Rscript

library(tigger)
library(dplyr)
library(stringi)

gene_usages <- read.delim("${gene_usages_file}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

TRBV_GERM <- readIgFasta("${germline_v}")

geno_BV <- read.delim("${genotype_v}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

DATA <- read.delim("${airrseq}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

DATA <- DATA[DATA[["consensus_count"]] >= ${min_consensus_count}, ]

cutoff <- 0.0005
p_val_cutoff <- 0.001

DATA[["v_gene"]] <- unlist(lapply(DATA[["v_call"]], function(Vcall){
  assignments <- unlist(strsplit(Vcall, ",", fixed = TRUE))
  v_genes <- unique(sapply(strsplit(assignments, "*", fixed = TRUE), "[", 1))
  v_genes <- v_genes[order(v_genes)]
  paste(v_genes, collapse = ",")
}))

DATA <- DATA[!grepl(",", DATA[["v_gene"]]), ]
DATA <- DATA[grepl("TRBV", DATA[["v_gene"]]), ]

gene_usages[["N"]] <- unlist(lapply(gene_usages[["GENE"]], function(gene){nrow(DATA[DATA[["v_gene"]] == gene, ])}))
gene_usages[["TOTAL"]] <- nrow(DATA)
gene_usages[["USAGE"]] <- gene_usages[["N"]] / gene_usages[["TOTAL"]]

gene_usages <- gene_usages[gene_usages[["MIN_FREQ"]] != Inf & gene_usages[["AVG_USAGE"]] > 1.5 * cutoff, ]

# calculate the p value for each gene in case that the gene is deleted by binom test
gene_usages[["PVAL"]] <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages[["USAGE"]][i] < cutoff) {
    return(binom.test(x = gene_usages[["N"]][i], n = gene_usages[["TOTAL"]][i], p = gene_usages[["MIN_FREQ"]][i])\$p.value)
  } else {
    return(1)
  }
})

# Detect according to the p values if there are deleted genes
gene_usages[["DELETED"]] <- sapply(1:nrow(gene_usages), function(i) {
  if (gene_usages[["PVAL"]][i] <= p_val_cutoff) {
    if ((gene_usages[["USAGE"]][i] < cutoff) & gene_usages[["MIN_FREQ"]][i] != Inf) {
      return(TRUE)
    }
  }
  return(FALSE)
})

gene_usages <- gene_usages[gene_usages[["DELETED"]], ]

if (nrow(gene_usages) > 0) {
  deleted_genes <- gene_usages[["GENE"]]
  geno_BV <- geno_BV[!geno_BV[["gene"]] %in% deleted_genes, ]
  
  for (gene in deleted_genes) {
    geno_BV[nrow(geno_BV) + 1, ] <- c(gene, NA, NA, NA, NA, NA, NA, NA, NA, 1000, "Deletion")
  }
}

write.table(geno_BV, file = paste0("${name}","_v_genotype_deletion.tsv"), quote = F, row.names = F, sep = "\t")

"""

}


process trb_genotype_report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outname}_genotype.tsv$/) "genotype_report/$filename"}
input:
 set val(name), file(airrseq) from g_30_outputFileTSV0_g_35
 set val(namev), file(genotype_v) from g_34_outputFileTSV0_g_35
 set val(named), file(genotype_d) from g_25_outputFileTSV4_g_35
 set val(namej), file(genotype_j) from g_25_outputFileTSV5_g_35

output:
 set val(name), file("${outname}_genotype.tsv")  into g_35_outputFileTSV0_g_37

script:
	
outname = airrseq.name.substring(0, airrseq.name.indexOf("_db-pass"))

"""
#!/usr/bin/env Rscript

library(tigger)
library(dplyr)
library(stringi)

DATA <- read.delim("${airrseq}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BV <- read.delim("${genotype_v}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BD <- read.delim("${genotype_d}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
geno_BJ <- read.delim("${genotype_j}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)



v_table <- DATA[!grepl(",", DATA[["v_call"]]), ]
v_table <- v_table %>% dplyr::group_by(v_call) %>% dplyr::summarise(N = n())
v_table[["gene"]] <- sapply(strsplit(v_table[["v_call"]], "*", fixed = TRUE), "[", 1)
v_table[["allele"]] <- sapply(strsplit(v_table[["v_call"]], "*", fixed = TRUE), "[", 2)
geno_BV[["Freq_by_Seq"]] <- unlist(lapply(geno_BV[["gene"]], function(g) {
  if (!g %in% v_table[["gene"]]) {
    return("0")
  }
  counts <- v_table[["N"]][v_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

d_table <- DATA[(!grepl(pattern = ',', DATA[["d_call"]]) & DATA[["d_call"]] != 'None') & (DATA[["d_sequence_end"]] - DATA[["d_sequence_start"]] >= 8), ]
d_table <- d_table[complete.cases(d_table[["sequence_id"]]), ]
d_table <- d_table %>% dplyr::group_by(d_call) %>% dplyr::summarise(N = n())
d_table[["gene"]] <- sapply(strsplit(d_table[["d_call"]], "*", fixed = TRUE), "[", 1)
d_table[["allele"]] <- sapply(strsplit(d_table[["d_call"]], "*", fixed = TRUE), "[", 2)
geno_BD[["Freq_by_Seq"]] <- unlist(lapply(geno_BD[["gene"]], function(g) {
  if (!g %in% d_table[["gene"]]) {
    return("0")
  }
  counts <- d_table[["N"]][d_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

j_table <- DATA[!grepl(",", DATA[["j_call"]]), ]
j_table <- j_table %>% dplyr::group_by(j_call) %>% dplyr::summarise(N = n())
j_table[["gene"]] <- sapply(strsplit(j_table[["j_call"]], "*", fixed = TRUE), "[", 1)
j_table[["allele"]] <- sapply(strsplit(j_table[["j_call"]], "*", fixed = TRUE), "[", 2)
geno_BJ[["Freq_by_Seq"]] <- unlist(lapply(geno_BJ[["gene"]], function(g) {
  if (!g %in% j_table[["gene"]]) {
    return("0")
  }
  counts <- j_table[["N"]][j_table[["gene"]] == g]
  counts <- sort(counts, decreasing = TRUE)
  return(paste(counts, collapse = ";"))
}))

geno <- rbind(geno_BV, geno_BD, geno_BJ)
geno[["Freq_by_Clone"]] <- geno[["Freq_by_Seq"]]
write.table(geno, file = paste0("${outname}","_genotype.tsv"), quote = FALSE, row.names = F, sep = "\t")
"""
}


process trb_haplotype {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_haplo_D2.tsv$/) "haplotype/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_haplo_J1_6.tsv$/) "haplotype/$filename"}
input:
 set val(name),file(airrseq) from g_30_outputFileTSV0_g_37
 set val(name1),file(genos) from g_35_outputFileTSV0_g_37
 set val(name2),file(germline_d) from g_2_germlineFastaFile_g_37

output:
 set val(name),file("*_haplo_D2.tsv") optional true  into g_37_outputFileTSV00
 set val(name),file("*_haplo_J1_6.tsv") optional true  into g_37_outputFileTSV11

script:
	
"""
#!/usr/bin/env Rscript

library(tigger)
library(dplyr)
library(stringi)
library(rabhit)

TRBV_GERM <- readIgFasta("${germline_v}")
TRBD_GERM <- readIgFasta("${germline_d}")
geno <- read.delim("${genos}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DATA <- read.delim("${airrseq}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DATA[["subject"]] <- "${name}"

hetero_j16 <- "TRBJ1-6" %in% geno[["gene"]]
if (hetero_j16) {
  hetero_j16 <- grepl(",", geno[["genotyped_alleles"]][geno[["gene"]] == "TRBJ1-6"])
}

hetero_d2 <- "TRBD2" %in% geno[["gene"]]
if (hetero_d2) {
  hetero_d2 <- grepl(",", geno[["genotyped_alleles"]][geno[["gene"]] == "TRBD2"])
}

# Filtering
DATA <- DATA[!grepl(",", DATA[["v_call"]]), ] # V single assignment

# zero mutations over the V
v_seqs <- sapply(1:nrow(DATA), function(x) substr(DATA[["sequence_alignment"]][x], 1, DATA[["v_germline_end"]][x]))
DATA[["v_mut"]] <- unlist(tigger::getMutCount(v_seqs, DATA[["v_call"]], germline_db = TRBV_GERM))
DATA <- DATA[DATA[["v_mut"]] <= 1, ]

del_genes <- geno[["gene"]][grepl("Del", geno[["genotyped_alleles"]])]

if (hetero_j16) {
  # Only rearrangements with single assignment of TRBJ1-6
  DATA_J <- DATA[grepl("J1-6", DATA[["j_call"]]), ]
  DATA_J <- DATA_J[!grepl(",", DATA_J[["j_call"]]), ]
  DATA_J[["J_SEQ_LENGTH"]] <- DATA_J[["j_sequence_end"]] - DATA_J[["j_sequence_start"]] + 1
  DATA_J <- DATA_J[DATA_J[["J_SEQ_LENGTH"]] > 10, ]
  
  TRBJ1_6_01_REA <- nrow(DATA_J[DATA_J[["j_call"]] == "TRBJ1-6*01", ])
  total <- nrow(DATA_J)
  
  if ((TRBJ1_6_01_REA / total > 0.3) & (TRBJ1_6_01_REA / total < 0.7)) {
    haplo_j1_6 <- createFullHaplotype(DATA_J, toHap_col = "v_call", hapBy_col = "j_call", chain = "TRB", hapBy = "TRBJ1-6", toHap_GERM = TRBV_GERM, rmPseudo = FALSE)
    haplo_j1_6[["TOTAL"]] <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_j1_6[haplo_j1_6[["gene"]] == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_j1_6, file = paste0("${name}","_haplo_J1_6.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

if (hetero_d2) {
  # Only rearrangements with single assignment of TRBD2
  DATA[["d_germline_length"]] <- as.integer(DATA[["d_germline_end"]]) - as.integer(DATA[["d_germline_start"]]) + 1
  DATA_D_geno <- DATA[(!grepl(pattern = ',', DATA[["d_call"]]) & DATA[["d_call"]] != 'None') & (DATA[["d_germline_length"]] >= 9), ]
  DATA_D_geno <- DATA_D_geno[complete.cases(DATA_D_geno[["sequence_id"]]), ]
  DATA_D_geno <- DATA_D_geno[grepl("D2", DATA_D_geno[["d_call"]]), ]
  
  # filter by zero mutations over the D segment
  # extract d sequence in the direct orientation
  DATA_D_reg <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] < DATA_D_geno[["d_germline_end"]], ]
  DATA_D_reg[["d_seq"]] <- substr(DATA_D_reg[["sequence"]], DATA_D_reg[["d_sequence_start"]], DATA_D_reg[["d_sequence_end"]])
  
  # extract convert d sequence in the inverted orientation to the direct orientation
  DATA_D_inv <- DATA_D_geno[DATA_D_geno[["d_germline_start"]] > DATA_D_geno[["d_germline_end"]], ]
  DATA_D_inv[["d_seq"]] <- substr(DATA_D_inv[["sequence"]], DATA_D_inv[["d_sequence_start"]], DATA_D_inv[["d_sequence_end"]])
  DATA_D_inv[["d_seq"]] <- stringi::stri_reverse(DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("A", "t", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("T", "a", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("G", "c", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- gsub("C", "g", DATA_D_inv[["d_seq"]])
  DATA_D_inv[["d_seq"]] <- toupper(DATA_D_inv[["d_seq"]])
  
  d_germ_end <- DATA_D_inv[["d_germline_start"]]
  DATA_D_inv[["d_germline_start"]] <- DATA_D_inv[["d_germline_end"]]
  DATA_D_inv[["d_germline_end"]] <- d_germ_end
  
  DATA_D_geno <- rbind(DATA_D_reg, DATA_D_inv)
  
  DATA_D_geno[["mut_d"]] <- unlist(lapply(1:nrow(DATA_D_geno), function(i) {
    mut <- 0
    row_seq <- unlist(strsplit(DATA_D_geno[["d_seq"]][[i]], ""))
    allele_seq <- unlist(strsplit(TRBD_GERM[[DATA_D_geno[["d_call"]][[i]]]], ""))
    for (pos in DATA_D_geno[["d_germline_start"]][[i]]:(DATA_D_geno[["d_germline_end"]][[i]])) {
      if (row_seq[pos - (DATA_D_geno[["d_germline_start"]][[i]] - 1)] != allele_seq[pos]) {
        mut <- mut + 1
      }
    }
    mut
  }))
  
  DATA_D_geno <- DATA_D_geno[DATA_D_geno[["mut_d"]] == 0, ]
  TRBD2_01_REA <- nrow(DATA_D_geno[DATA_D_geno[["d_call"]] == "TRBD2*01", ])
  total <- nrow(DATA_D_geno)
  
  if ((TRBD2_01_REA / total > 0.3) & (TRBD2_01_REA / total < 0.7)) {
    haplo_d2 <- createFullHaplotype(DATA_D_geno, toHap_col = "v_call", hapBy_col = "d_call", chain = "TRB", hapBy = "TRBD2", toHap_GERM = TRBV_GERM, rmPseudo = FALSE, kThreshDel = 5)
    haplo_d2[["TOTAL"]] <- total
    
    if (length(del_genes) > 0) {
      for (gene in del_genes) {
        haplo_d2[haplo_d2[["gene"]] == gene, c(3:5, 9)] <- c("Del", "Del", "Del", 1000)
      }
    }
    
    write.table(haplo_d2, file =  paste0("${name}","_haplo_D2.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
