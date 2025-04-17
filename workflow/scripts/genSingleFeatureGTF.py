# Generate hybrid single feature GTF by combining all input fasta contigs
# Author: Moe

import gffutils
import sys
import argparse
from Bio import SeqIO
import singleFeature
import gzip


parser = argparse.ArgumentParser(description="Process multiple input files.")
parser.add_argument("-o", "--output", type=str, required=True, help="Onput file path")
parser.add_argument("input_files", type=str, nargs="+", help="Input file paths")
args = parser.parse_args()

db=""
for fastafile in args.input_files:
    with gzip.open(fastafile, 'rt') as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            seqid = record.id
            print(f"seqid: {seqid}")
            sequence = str(record.seq)
            print(f"sequence: {sequence}")
            start = 1
            end = len(sequence)
            print(f"sequence length: {len(sequence)}")
            new_feature = singleFeature.single_feature(seqid=seqid, end=end)
            if db:
                db.update(new_feature)
            else:
                db = gffutils.create_db(new_feature, ":memory:")

output_gtf = args.output
with open(output_gtf, 'w') as out_gtf:
    for feature in db.all_features(order_by=('seqid', 'start')):
        attr_str = singleFeature.format_attributes(dict(feature.attributes))
        gtf_line = '\t'.join([
            feature.seqid,
            feature.source,
            feature.featuretype,
            str(feature.start),
            str(feature.end),
            str(feature.score) if feature.score is not None else '.',
            feature.strand if feature.strand is not None else '.',
            str(feature.frame) if feature.frame is not None else '.',
            attr_str
        ])
        out_gtf.write(gtf_line + '\n')
