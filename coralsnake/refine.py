#!/usr/bin/env python
import gzip, re, os, subprocess, logging
from collections import defaultdict, OrderedDict

def run_cmd(cmd):
    logging.info("RUN CMD: %s" % cmd)
    subprocess.check_call(cmd, shell=True)

def load_seqname_mapper(path, is_ncbi_report=False):
    rename_mapper = dict()
    for line in open(path):
        if line.startswith("#"):
            continue
        row = line.strip().split("\t")
        if is_ncbi_report:
                rename_mapper[row[6]] = row[0].replace(" ", "_")
        else:
            rename_mapper[row[0]] = row[1].replace(" ", "_")
    return rename_mapper


def load_canonical_transcripts(path, is_ucsc_known_canonical=False):
    canonical_transcripts = set()
    if is_ucsc_known_canonical:
        for i, line in enumerate(open(path)):
            row = line.strip().split("\t")
            if i == 0:
                assert row[4] == "transcript"
            else:
                v = row[4]
                if v.startswith("ENS"):
                    v = v.split(".")[0]
                canonical_transcripts.add(v)
    else:
        for line in open(path): 
            row = line.strip().split("\t")
            canonical_transcripts.add(row[0])
    return canonical_transcripts
    

class FastaRefiner(object):
    def __init__(self, input_fasta, output_prefix, seqname_mapper=None, seqname_pattern=None):
        self.input_fasta = input_fasta
        self.output_fasta = output_prefix + ".genome.fasta"
        self.output_sizes = output_prefix + ".genome.sizes"
        self.seqname_mapper = seqname_mapper
        self.seqname_pattern = seqname_pattern
        
    def run(self):
        logging.info("Refine genome fasta file.")
        f = gzip.open(self.input_fasta, "rt") if self.input_fasta.endswith(".gz") else open(self.input_fasta, "r")
        with open(self.output_fasta, "w+") as fw:
            b = False
            for line in f:
                if line.startswith(">"):
                    s = line.strip("\n")[1:].strip()
                    i = s.find(" ")
                    description = None
                    if i == -1:
                        seqname_old = s
                    else:
                        seqname_old = s[:i]
                        description = s[i + 1:].strip()
                       
                    seqname_new = seqname_old # new seqname
                    if self.seqname_mapper is not None:
                        if seqname_old in self.seqname_mapper:
                            seqname_new = self.seqname_mapper[seqname_old]
                        else:
                            logging.warning(f"{seqname_old} is not in seqname mapper.")
                            seqname_new = seqname_old
                        
                    if self.seqname_pattern is not None:
                        if re.search(self.seqname_pattern, seqname_new) is None:
                            b = False
                            continue
                        
                    if description is None:
                        line = ">%s %s\n" % (seqname_new, seqname_old)
                    else:
                        line = ">%s %s %s\n" % (seqname_new, seqname_old, description)
                        
                    fw.write(line)
                    b = True
                else:                
                    if b:
                        fw.write(line)
        f.close()
        run_cmd("samtools faidx %s" % self.output_fasta)
        run_cmd("cut -f 1,2 %s.fai > %s" % (self.output_fasta, self.output_sizes))


class GtfUtils(object):
    @staticmethod
    def group_rows_by_gene_id(rows):
        data = defaultdict(list)
        for row in rows:
            data[row[-1]["gene_id"]].append(row)
        return data

    @staticmethod
    def group_rows_by_transcript_id(rows):
        data = defaultdict(list)
        for row in rows:
            if row[2] == "gene":
                continue
            data[row[-1]["transcript_id"]].append(row)
        return data

    @staticmethod
    def group_rows_by_feature(rows):
        data = defaultdict(list)
        for row in rows:
            data[row[2]].append(row)
        return data

    @staticmethod
    def check_seqname_and_strand_consistency(rows):
        row1 = rows[0]
        gid = row1[8]["gene_id"]
        b = True
        for i in range(1, len(rows)):
            row2 = rows[i]
            if row1[0] != row2[0]:
                logging.error(f"Gene ID {gid} has rows on different seqnames: {row1[0]} and {row2[0]}")
                b = False
                break
            if row1[6] != row2[6]:
                logging.error(f"Gene ID {gid} has rows on different strands: {row1[6]} and {row2[6]}")
                b = False
                break
        return b
    
    @staticmethod
    def parse_gtf_attribute_string(s):
        attributes = OrderedDict()
        key_x, key_y = None, None
        value_x, value_y = None, None
        value_has_quotes = False
        i = 0
        while i < len(s):
            c = s[i]
            if key_x is None:
                if c != " " and c != ";":
                    key_x = i
            else:
                if key_y is None:
                    if c == " ":
                        key_y = i
                else:
                    if value_x is None:
                        if c != " ":
                            if c == "\"":
                                value_x = i + 1
                                value_has_quotes = True
                            else:
                                value_x = i
                                value_has_quotes = False
                    else:
                        if value_has_quotes:
                            # if c == "\"" and s[i - 1] != "\\":
                            #     assert s[i + 1] == ";"
                            if c == "\"" and s[i - 1] != "\\" and s[i + 1] == ";":
                                value_y = i
                                attributes[s[key_x:key_y]] = s[value_x:value_y]
                                key_x, key_y = None, None
                                value_x, value_y = None, None
                        else:
                            if c == ";":
                                value_y = i
                                attributes[s[key_x:key_y]] = s[value_x:value_y]
                                key_x, key_y = None, None
                                value_x, value_y = None, None
            i += 1
        return attributes
    
    @staticmethod
    def load_gtf(path, seqname_mapper=None, seqname_pattern=None):
        header = []
        rows = []
        f = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
        for line in f:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                row = line.strip().split("\t")
                assert len(row) == 9
                if seqname_mapper is not None:
                    if row[0] in seqname_mapper:
                        row[0] = seqname_mapper[row[0]]
                    else:
                        logging.warning(f"{row[0]} is not in seqname mapper.")
                if seqname_pattern is not None:
                    if re.search(seqname_pattern, row[0]) is None:
                        continue
                row[1] = row[1].replace(" ", "_")
                row[3] = int(row[3])
                row[4] = int(row[4])
                assert row[3] >= 1
                assert row[4] >= row[3]
                row[8] = GtfUtils.parse_gtf_attribute_string(row[8])
                assert "gene_id" in row[8]
                if row[2] == "gene" and "transcript_id" in row[8]:
                    del row[8]["transcript_id"]
                if row[2] != "gene":
                    assert "transcript_id" in row[8]
                if "exon_number" in row[8]:
                    del row[8]["exon_number"]
                rows.append(row)
        f.close()
        return header, rows

    @staticmethod
    def create_gene_row_from_gene_rows(rows):
        gene_id = rows[0][8]["gene_id"]
        data_features = GtfUtils.group_rows_by_feature(rows)
        if "transcript" in data_features or "exon" in data_features:
            logging.info(f"Create gene row from transcript and exon rows for {gene_id}")
            gene_row = None
            attris = None
            if "transcript" in data_features:
                for row in data_features["transcript"]:
                    if attris is None:
                        attris = row[8].copy()
                        del attris["transcript_id"]
                    if gene_row is None:
                        gene_row = row[:8]
                        gene_row[2] = "gene"
                        gene_row.append(attris)
                    else:
                        gene_row[3] = min(gene_row[3], row[3])
                        gene_row[4] = max(gene_row[4], row[4])
            if "exon" in data_features:
                for row in data_features["transcript"]:
                    if attris is None:
                        attris = row[8].copy()
                        del attris["transcript_id"]
                    if gene_row is None:
                        gene_row = row[:8]
                        gene_row[2] = "gene"
                        gene_row.append(attris)
                    else:
                        gene_row[3] = min(gene_row[3], row[3])
                        gene_row[4] = max(gene_row[4], row[4])
            return gene_row
        elif "CDS" in data_features:
            logging.info(f"Create gene row from CDS rows for {gene_id}")
            gene_row = None
            attris = None
            for row in data_features["CDS"]:
                if attris is None:
                    attris = row[8].copy()
                    del attris["transcript_id"]
                if gene_row is None:
                    gene_row = row[:8]
                    gene_row[2] = "gene"
                    gene_row[7] = "."
                    gene_row.append(attris)
                else:
                    gene_row[3] = min(gene_row[3], row[3])
                    gene_row[4] = max(gene_row[4], row[4])
            strand = gene_row[6]
            if strand == "-":
                gene_row[3] -= 3
                assert gene_row[3] >= 1
            else:
                gene_row[4] += 3
            return gene_row
        else:
            return None
            
    @staticmethod
    def infer_gene_type(rows):
        gid = rows[0][8]["gene_id"]
        vs: list[str] = []
        for row in rows:
            attris = row[8]
            if "gene_type" in attris:
                vs.append(attris["gene_type"])
            if "gene_biotype" in attris:
                vs.append(attris["gene_biotype"])
        vs = list(set(vs))
        if len(vs) == 1:
            return vs[0]
        elif len(vs) > 1:
            logging.info(f"Gene ID {gid} has multiple gene types: {vs}")
            return None
        else:
            for row in rows:
                if row[2] == "CDS":
                    return "protein_coding"
            return "unknown"

    @staticmethod
    def infer_gene_name(rows):
        gid = rows[0][8]["gene_id"]
        vs: list[str] = []
        for row in rows:
            attris = row[8]
            if "gene_name" in attris:
                vs.append(attris["gene_name"])
            if "gene" in attris:
                vs.append(attris["gene"])
        vs = list(set(vs))
        if len(vs) == 1:
            return vs[0]
        elif len(vs) > 1:
            logging.info(f"Gene ID {gid} has multiple gene names: {vs}")
            return None
        else:
            return rows[0][8]["gene_id"]
        
    @staticmethod
    def infer_transcript_type(rows):
        tid = rows[0][8]["transcript_id"]
        vs = []
        for row in rows:
            attris = row[8]
            if "transcript_type" in attris:
                vs.append(attris["transcript_type"])
            if "transcript_biotype" in attris:
                vs.append(attris["transcript_biotype"])
        vs = list(set(vs))
        if len(vs) == 1:
            return vs[0]
        elif len(vs) > 1:
            logging.info(f"Transcript ID {tid} has multiple transcript types: {vs}")
            return None
        else:
            for row in rows:
                if row[2] == "CDS":
                    return "protein_coding"
            return "unknown"
        
    @staticmethod
    def infer_transcript_name(rows):
        tid = rows[0][8]["transcript_id"]
        vs: list[str] = []
        for row in rows:
            attris = row[8]
            if "transcript_name" in attris:
                vs.append(attris["transcript_name"])
        vs = list(set(vs))
        if len(vs) == 1:
            return vs[0]
        elif len(vs) > 1:
            logging.info(f"Transcript ID {tid} has multiple transcript names: {vs}")
            return None
        else:
            return tid
        
    @staticmethod
    def output_row(fw, row, check=False):
        assert len(row) == 9
        attris = row[8]
        if check:
            assert "gene_id" in attris
            assert "gene_name" in attris
            assert "gene_type" in attris
            if row[2] == "gene":
                assert "transcript_id" not in attris
                assert "exon_number" not in attris
            else:
                assert "transcript_id" in attris
                assert "transcript_name" in attris
                assert "transcript_type" in attris
                if row[2] == "transcript":
                    assert "exon_number" not in attris
                elif row[2] == "exon" or row[2] == "CDS":
                    assert "exon_number" in attris
        items = []
        for k, v in attris.items():
            items.append('%s "%s";' % (k, v))
        s = " ".join(items)
        fw.write("%s\t%s\n" % ("\t".join([str(x) for x in row[:8]]), s))


    @staticmethod
    def output_rows(fw, rows, check=False):
        for row in rows:
            GtfUtils.output_row(fw, row, check)


    @staticmethod
    def create_transcript_row_from_transcript_rows(rows):
        data_features = GtfUtils.group_rows_by_feature(rows)
        if "exon" in data_features:
            transcript_row = None
            attris = None
            for row in data_features["exon"]:
                if attris is None:
                    attris = row[8].copy()
                    if "exon_number" in attris:
                        del attris["exon_number"]
                if transcript_row is None:
                    transcript_row = row[:8]
                    transcript_row[2] = "transcript"
                    transcript_row.append(attris)
                else:
                    transcript_row[3] = min(transcript_row[3], row[3])
                    transcript_row[4] = max(transcript_row[4], row[4])
            return transcript_row
        elif "CDS" in data_features:
            transcript_row = None
            attris = None
            for row in data_features["CDS"]:
                if attris is None:
                    attris = row[8].copy()
                    if "exon_number" in attris:
                        del attris["exon_number"]
                if transcript_row is None:
                    transcript_row = row[:8]
                    transcript_row[2] = "transcript"
                    transcript_row[7] = "."
                    transcript_row.append(attris)
                else:
                    transcript_row[3] = min(transcript_row[3], row[3])
                    transcript_row[4] = max(transcript_row[4], row[4])
            # Extend to include stop codon
            strand = transcript_row[6]
            if strand == "-":
                transcript_row[3] -= 3
            else:
                transcript_row[4] += 3
            return transcript_row
        else:
            assert False
    @staticmethod
    def create_transcript_row_from_gene_row(gene_row):
        gene_attris = gene_row[8]
        transcript_id = gene_attris["gene_id"] + "_transcript"
        transcript_type = gene_attris["gene_type"]
        transcript_row = gene_row[:8]
        transcript_row[2] = "transcript"
        transcript_attris = gene_attris.copy()
        transcript_attris["transcript_id"] = transcript_id
        transcript_attris["transcript_name"] = transcript_id
        transcript_attris["transcript_type"] = transcript_type
        transcript_row.append(transcript_attris)
        return transcript_row
    @staticmethod
    def create_exon_row_from_transcript_row(transcript_row):
        transcript_attris = transcript_row[8]
        exon_row = transcript_row[:8]
        exon_row[2] = "exon"
        exon_attris = transcript_attris.copy()
        exon_attris["exon_number"] = 1
        exon_row.append(exon_attris)
        return exon_row
    
    @staticmethod
    def create_exon_rows_from_cds_rows(cds_rows):
        exon_rows = []
        strand = cds_rows[0][6]
        for i, cds_row in enumerate(cds_rows):
            exon_row = cds_row[:8]
            exon_row[2] = "exon"
            exon_row[7] = "."
            exon_attris = cds_row[8].copy()
            exon_row.append(exon_attris)
            if strand == "-":
                if i == 0:
                    exon_row[3] -= 3
            else:
                if i == len(cds_rows) - 1:
                    exon_row[4] += 3
            exon_rows.append(exon_row)
        return exon_rows
    
    @staticmethod    
    def create_start_codon_from_cds_rows(exon_rows, cds_rows):
        # TODO
        raise NotImplementedError()
        # if strand == "-":
        #     cds_row = cds_row[-1]
        #     codon_row = cds_row[:8]
        #     codon_row[3] = codon_row[4] - 2
        # else:
        #     cds_row = cds_rows[0]
        #     codon_row = cds_row[:8]
        # codon_row[4] = codon_row[3] + 2
        # codon_row[2] = "start_codon"
        # codon_attris = cds_row[8].copy()
        # del codon_attris["exon_number"]
        # codon_row.append(codon_attris)
        # output_row(fw, codon_row, check=True)
        
    @staticmethod
    def create_stop_codon_from_cds_rows(exon_rows, cds_rows):
        # TODO
        raise NotImplementedError()
        # if strand == "-":
        #     cds_row = cds_row[0]
        #     codon_row = cds_row[:8]
        #     codon_row[3] = codon_row[3] - 3
        # else:
        #     cds_row = cds_rows[-1]
        #     codon_row = cds_row[:8]
        #     codon_row[3] = codon_row[4] + 1
        # codon_row[4] = codon_row[3] + 2
        # codon_row[2] = "stop_codon"
        # codon_attris = cds_row[8].copy()
        # del codon_attris["exon_number"]
        # codon_row.append(codon_attris)
        

    @staticmethod
    def merge_overlapping_exon_rows(exon_rows):
        tid = exon_rows[0][8]["transcript_id"]
        i = 0
        while i < len(exon_rows) - 1:
            row1, row2 = exon_rows[i], exon_rows[i + 1]
            if row1[4] + 1 >= row2[3]:
                row1[4] = max(row1[4], row2[4])
                exon_rows.pop(i + 1)
                logging.warning(f"Transcript ID {tid} has overlap exons. {row1[3]}-{row1[4]} and {row2[3]}-{row2[4]}. Automatic merge them.")
            else:
                i += 1
        return exon_rows
    @staticmethod
    def merge_overlapping_cds_rows(cds_rows):
        tid = cds_rows[0][8]["transcript_id"]
        i = 0
        while i < len(cds_rows) - 1:
            row1, row2 = cds_rows[i], cds_rows[i + 1]
            if row1[4] + 1 >= row2[3]:
                logging.warning(f"Transcript ID {tid} has overlap CDS. {row1[3]}-{row1[4]} and {row2[3]}-{row2[4]}. Automatic merge them.")
                row1[4] = max(row1[4], row2[4])
                cds_rows.pop(i + 1)
            else:
                i += 1
        return cds_rows
    @staticmethod
    def check_cds_consistency_with_exon(exon_rows, cds_rows):
        offset = None
        for ci, cds_row in enumerate(cds_rows):
            # find matched exon
            if ci == 0:
                for ei, exon_row in enumerate(exon_rows):
                    if exon_row[4] < cds_row[3]:
                        continue
                    elif exon_row[3] > cds_row[4]:
                        return False
                    else:
                        offset = ei
                        break
            
            exon_row = exon_rows[offset + ci]
            if ci == 0:
                if ci == len(cds_rows) - 1:
                    if not (cds_row[3] >= exon_row[3] and cds_row[4] <= exon_row[4]):
                        return False
                else:
                    if not (cds_row[3] >= exon_row[3] and cds_row[4] == exon_row[4]):
                        return False
            else:
                if ci == len(cds_rows) - 1:
                    if not (cds_row[3] == exon_row[3] and cds_row[4] <= exon_row[4]):
                        return False
                else:
                    if not (cds_row[3] == exon_row[3] and cds_row[4] == exon_row[4]):
                        return False
        return True
     
    @staticmethod
    def assign_exon_number_for_exon_rows(exon_rows):
        strand = exon_rows[0][6]
        for i, exon_row in enumerate(exon_rows):
            if strand == "-":
                exon_row[8]["exon_number"] = len(exon_rows) - i
            else:
                exon_row[8]["exon_number"] = i + 1 
    @staticmethod
    def assign_exon_number_for_cds_rows(exon_rows, cds_rows):
        transcript_id = exon_rows[0][8]["transcript_id"]
        offset = None
        for ci, cds_row in enumerate(cds_rows):
            # find matched exon
            if ci == 0:
                for ei, exon_row in enumerate(exon_rows):
                    if exon_row[4] < cds_row[3]:
                        continue
                    elif exon_row[3] > cds_row[4]:
                        logging.error(f"Transcript ID {transcript_id} occurs CDS don't match exon!")
                        return False
                    else:
                        offset = ei
                        break
            exon_row = exon_rows[offset + ci]
            cds_row[8]["exon_number"] = exon_row[8]["exon_number"]

    @staticmethod
    def get_canonical_transcript_id(final_gene_rows, canonicals=None):
        transcript_lengths = defaultdict(int)
        for row in final_gene_rows:
            if row[2] == "exon":
                transcript_lengths[row[8]["transcript_id"]] += (row[4] - row[3] + 1)
            
        tids1, tids2, tids3, tids4 = [], [], [], []
        for row in final_gene_rows:
            if row[2] == "transcript":
                transcript_id = row[8]["transcript_id"]
                transcript_type = row[8]["transcript_type"]
                if canonicals is not None and (transcript_id in canonicals or transcript_id.split(".")[0] in canonicals):
                    if transcript_type in ["protein_coding", "mRNA"]:
                        tids1.append(transcript_id)
                    else:
                        tids2.append(transcript_id)
                if transcript_type in ["protein_coding", "mRNA"]:
                    tids3.append(transcript_id)
                else:
                    tids4.append(transcript_id)
        if len(tids1) > 0:
            vs = tids1
        elif len(tids2) > 0:
            vs = tids2
        elif len(tids3) > 0:
            vs = tids3
        else:
            vs = tids4
        # print(tids1)
        vs = list(sorted(vs, key=lambda x: transcript_lengths[x], reverse=True))
        return vs[0]
    

class GtfRefiner(object):
    def __init__(self, input_gtf, output_prefix, seqname_mapper=None, seqname_pattern=None, canonicals=None):
        self.input_gtf = input_gtf
        self.output_gtf = output_prefix + ".annotation.gtf"
        self.output_skip_gtf = output_prefix + ".skip.gtf"
        self.feature_summary_txt = output_prefix + ".gene_features_summary.txt"
        self.seqname_mapper = seqname_mapper
        self.seqname_pattern = seqname_pattern
        self.canonicals = canonicals
        
        self.exist_gene_names = set()
        self.exist_transcript_ids = set()
        self.exist_transcript_names = set()
        
        self.fw = None
        self.fw_skip = None
        
            
    def run(self):
        self.fw = open(self.output_gtf, "w+")
        self.fw_skip = open(self.output_skip_gtf, "w+")
        
        logging.info(f"Loding rows from {self.input_gtf}")
        header, rows = GtfUtils.load_gtf(self.input_gtf, seqname_mapper=self.seqname_mapper, seqname_pattern=self.seqname_pattern)
        logging.info(f"Loaded {len(rows)} rows from GTF file")
        for line in header:
            self.fw.write(line + "\n")
        self.fw.write("#!refined gtf\n")
        
        data_gene = GtfUtils.group_rows_by_gene_id(rows)
        total_genes = len(data_gene)
        logging.info(f"Grouped rows into {total_genes} genes")    
        
        
        features_counter = defaultdict(int)
        succeed_genes = 0
        for _, rows in data_gene.items():
            features = tuple(sorted(set([row[2] for row in rows])))
            features_counter[features] += 1
            if self.process_gene_rows(rows):
                succeed_genes += 1
        logging.info(f"Successfully refined {succeed_genes}/{total_genes} genes (%.2f%%)" % (succeed_genes * 100.0 / len(data_gene)))
        
        self.fw.close()
        self.fw_skip.close()
        
        assert self.output_gtf.endswith(".gtf")
        prefix = self.output_gtf[:-4]
        run_cmd("sort -k1,1 -k4,4n -k5,5n %s.gtf | grep -v '#' > %s.sorted.gtf" % (prefix, prefix))
        run_cmd("bgzip -c %s.sorted.gtf > %s.sorted.gtf.gz" % (prefix, prefix))
        run_cmd("tabix -p gff %s.sorted.gtf.gz" % prefix)
        
        
        with open(self.feature_summary_txt, "w+") as fw:
            total = sum(features_counter.values())
            fw.write("Total\tRatio\tFeatures\n")
            for features, count in features_counter.items():
                fw.write("%d\t%.6f\t%s\n" % (count, count / total, ",".join(features)))
                
    def process_gene_rows(self, rows):
        final_gene_rows = []
        
        gene_id = rows[0][8]["gene_id"]
        
        if not GtfUtils.check_seqname_and_strand_consistency(rows):
            GtfUtils.output_rows(self.fw_skip, rows, check=False)
            logging.warning(f"Gene ID {gene_id} has inconsistent rows. Skipping this gene!")
            return False
        
        gene_type = GtfUtils.infer_gene_type(rows)
        if gene_type is None:
            logging.warning(f"Gene ID {gene_id} has ambiguous gene types. Skipping this gene!")
            GtfUtils.output_rows(self.fw_skip, rows, check=False)
            return False
        
        gene_name = GtfUtils.infer_gene_name(rows)
        if gene_name is None:
            logging.warning(f"Gene ID {gene_id} has ambiguous gene names. Skipping this gene!")
            GtfUtils.output_rows(self.fw_skip, rows, check=False)
            return False
        
        if gene_name in self.exist_gene_names:
            gene_name_new = gene_name + "_" + gene_id
            assert gene_name_new not in self.exist_gene_names
            logging.warning(f"Gene ID {gene_id} has duplicate gene names. rename {gene_name} to {gene_name_new}")
            gene_name = gene_name_new
        self.exist_gene_names.add(gene_name)
        
        for row in rows:
            attris = row[8]
            attris["gene_type"] = gene_type
            attris["gene_name"] = gene_name
        
        data_features = GtfUtils.group_rows_by_feature(rows)
        gene_row = None
        if "gene" in data_features:
            if len(data_features["gene"]) > 1:
                logging.warning(f"Gene ID {gene_id} has multiple gene features. Skipping this gene!")
                GtfUtils.output_rows(self.fw_skip, rows, check=False)
                return False
            gene_row = data_features["gene"][0]
        else:
            logging.warning(f"Gene ID {gene_id} don't have gene feature, try to create gene from transcript, exon and CDS.")
            gene_row = GtfUtils.create_gene_row_from_gene_rows(rows)
        final_gene_rows.append(gene_row)
        
        for row in rows:
            assert row[3] >= gene_row[3]
            assert row[4] <= gene_row[4]
                
        data_transcript = GtfUtils.group_rows_by_transcript_id(rows)
        if len(data_transcript) == 0:
            logging.warning(f"Gene ID {gene_id} don't have any transcript, try to create transcript, exon from gene.")
            transcript_row = GtfUtils.create_transcript_row_from_gene_row(gene_row)
            exon_row = GtfUtils.create_exon_row_from_transcript_row(transcript_row)
            final_gene_rows.append(transcript_row)
            final_gene_rows.append(exon_row)
        else:
            for _, rows in data_transcript.items():
                final_transcript_rows = self.process_transcript_rows(rows)
                if final_transcript_rows is None:
                    return False
                else:
                    final_gene_rows.extend(final_transcript_rows)
        
        canonical_transcript_id = GtfUtils.get_canonical_transcript_id(final_gene_rows, self.canonicals)
        
        # Output rows

        for row in final_gene_rows:
            if row[2] != "gene":
                row[8]["is_canonical"] = row[8]["transcript_id"] == canonical_transcript_id
            GtfUtils.output_row(self.fw, row, check=True)
            
        return True
    
    def process_transcript_rows(self, rows):
        final_transcript_rows = []
        tid = rows[0][8]["transcript_id"]
        
        # Make transcript_id unique
        if tid in self.exist_transcript_ids:
            gid = rows[0][8]["gene_id"]
            tid_new = tid + "_" + gid
            assert tid_new not in self.exist_transcript_ids
            logging.warning(f"Transcript ID {tid} has duplicate transcript id. rename {tid} to {tid_new}")
            tid = tid_new
            for row in rows:
                row[8]["transcript_id"] = tid
        self.exist_transcript_ids.add(tid)
        
        # Infer transcript_type
        transcript_type = GtfUtils.infer_transcript_type(rows)
        if transcript_type is None:
            logging.error(f"Transcript ID {tid} has ambiguous transcript types. Please check your GTF file!")
            return None
            
        # Infer transcript_name. Must run after maing transcript_id unique, 
        # because transcript_name is based on transcript_id when 'transcript_name' does not exists.
        transcript_name = GtfUtils.infer_transcript_name(rows)
        if transcript_name is None:
            logging.error(f"Transcript ID {tid} has ambiguous transcript names. Please check your GTF file!")
            return None
        
        # Make transcript_name unique
        if transcript_name in self.exist_transcript_names:
            transcript_name_new = transcript_name + "_" + tid
            assert transcript_name_new not in self.exist_transcript_names
            logging.warning(f"Transcript ID {tid} has duplicate transcript names. rename {transcript_name} to {transcript_name_new}")
            transcript_name = transcript_name_new
        self.exist_transcript_names.add(transcript_name)           

        # Update attributes
        for row in rows:
            attris = row[8]
            attris["transcript_type"] = transcript_type
            attris["transcript_name"] = transcript_name
                        
        data_features = GtfUtils.group_rows_by_feature(rows)
        
        # Fetch transcript row
        if "transcript" in data_features:
            if len(data_features["transcript"]) > 1:
                logging.error(f"Transcript ID {tid} has multiple transcript features. Please check your GTF file!")
                return None
            transcript_row = data_features["transcript"][0]
        else:
            logging.warning(f"Transcript ID {tid} don't have transcript feature, try to create transcript from exon and CDS.")
            transcript_row = GtfUtils.create_transcript_row_from_transcript_rows(rows)
        final_transcript_rows.append(transcript_row)
        
        # Check coordinates
        for row in rows:
            if not (row[3] >= transcript_row[3] and row[4] <= transcript_row[4]):
                logging.error(f"Transcript ID {tid} has rows outside of transcript. Please check your GTF file!")
                return None
        
        # Check exon and CDS
        cds_rows = None
        if "CDS" in data_features:
            cds_rows = list(sorted(data_features["CDS"], key=lambda row: [row[3], row[4]]))
            cds_rows = GtfUtils.merge_overlapping_cds_rows(cds_rows)
                    
        exon_rows = None
        if "exon" in data_features:
            exon_rows = list(sorted(data_features["exon"], key=lambda row: [row[3], row[4]]))
            exon_rows = GtfUtils.merge_overlapping_exon_rows(exon_rows)
        else:     
            if "CDS" in data_features:
                exon_rows = GtfUtils.create_exon_rows_from_cds_rows(cds_rows)
            else:
                exon_rows = [GtfUtils.create_exon_row_from_transcript_row(transcript_row)]
        GtfUtils.assign_exon_number_for_exon_rows(exon_rows)
        final_transcript_rows.extend(exon_rows)
                
        if "CDS" in data_features:
            if not GtfUtils.check_cds_consistency_with_exon(exon_rows, cds_rows):
                logging.error(f"Transcript ID {tid} occurs CDS don't match exon!")
                return None
            GtfUtils.assign_exon_number_for_cds_rows(exon_rows, cds_rows)
                    
            # if "start_codon" in data_features:
            #     final_transcript_rows.extend(data_features["start_codon"])
            # else:
            #     final_transcript_rows.extend(GtfUtils.create_start_codon_from_cds_rows(exon_rows, cds_rows))
            
            # if "stop_codon" in data_features:
            #     final_transcript_rows.extend(data_features["stop_codon"])
            # else:
            #     final_transcript_rows.extend(GtfUtils.create_stop_codon_from_cds_rows(exon_rows, cds_rows))
                
            final_transcript_rows.extend(cds_rows)
            
        return final_transcript_rows
                            
def refine_genome_references(input_fasta, input_gtf, outdir, name, rename_mapper=None, seqname_pattern=None, canonical_transcripts=None, is_ucsc_known_canonical=False):
    name = os.path.basename(outdir) if name is None else name
    seqname_mapper = None if rename_mapper is None else load_seqname_mapper(rename_mapper)
    canonical_transcripts = None if canonical_transcripts is None else load_canonical_transcripts(canonical_transcripts, is_ucsc_known_canonical=is_ucsc_known_canonical)
    output_prefix = os.path.join(outdir, name)
    os.makedirs(outdir, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    with open(output_prefix + ".seqname_mapper.tsv", "w+") as fw:
        if seqname_mapper is not None:
            for k, v in seqname_mapper.items():
                fw.write("%s\t%s\n" % (k, v))
                
    with open(output_prefix + ".seqname_pattern.tsv", "w+") as fw:
        if seqname_pattern is not None:
            fw.write(seqname_pattern + "\n")
    
    if input_fasta is not None:
        FastaRefiner(
            input_fasta=input_fasta, 
            output_prefix=output_prefix, 
            seqname_mapper=seqname_mapper, 
            seqname_pattern=seqname_pattern).run()

    if input_gtf is not None:
        GtfRefiner(
            input_gtf=input_gtf, 
            output_prefix=output_prefix, 
            seqname_mapper=seqname_mapper, 
            seqname_pattern=seqname_pattern, 
            canonicals=canonical_transcripts).run()