#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import csv
import gzip

from Bio import SeqIO


def load_mapping(tsv_file):
    """
    Load mapping of UniProtKB proteins to Rhea IDs from a TSV file.
    The TSV file should have a header.
    Returns a dictionary with Protein ID as the key and RheaIDs as the value.
    """
    mapping = {}
    with open(tsv_file) as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader) # skip the header
        for row in reader:
            prot_id, rhea_id = row
            rhea_id = " ".join(sorted(set(rhea_id.split())))
            mapping[prot_id] = rhea_id
    return mapping


def filter_fasta(in_handle, out_handle, mapping):
    """
    Filter the FASTA file based on the mapping and add RheaID to the FASTA header.
    Raise an exception if the header is not in the format UniRef90_<prot_id>.
    """
    class InvalidProteinIDException(Exception):
        pass

    for record in SeqIO.parse(in_handle, "fasta"):
        if not record.id.startswith("UniRef90_") or len(record.id.split("_")) != 2:
            raise InvalidProteinIDException(
                f"Invalid protein ID format: {record.id}"
                )
        prot_id = record.id.split("_")[1]
        if prot_id in mapping:
            rhea_id = mapping[prot_id]
            record.description += f' RheaID="{rhea_id}"'
            SeqIO.write(record, out_handle, "fasta")


def processing_handle(input_fasta, output_fasta, mapping):
    """
    Enable processing of both regular and gzipped FASTA files.
    """
    with open(output_fasta, 'w') as out_handle:
        if input_fasta.endswith('.gz'):
            with gzip.open(input_fasta, 'rt') as in_handle:
                filter_fasta(in_handle, out_handle, mapping)
        else:
            filter_fasta(input_fasta, out_handle, mapping)

def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA keeping only proteins that have RheaIDs and add RheaID to the header."
        )
    parser.add_argument('input_fasta', type=str, help='Input FASTA file to be filtered (can be .fasta or .fasta.gz).')
    parser.add_argument('mapping_file', type=str, help='TSV file mapping UniProtKB proteins to Rhea IDs.')
    parser.add_argument('output_fasta', type=str, help='Output FASTA file with filtered records.')
    
    args = parser.parse_args()
    
    mapping = load_mapping(args.mapping_file)
    
    processing_handle(args.input_fasta, args.output_fasta, mapping)

if __name__ == '__main__':
    main()
