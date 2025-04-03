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
import fileinput

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


def filter_fasta(in_handle, mapping):
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
            yield record


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA keeping only proteins that have RheaIDs and add RheaID to the header."
        )
    parser.add_argument('input_fasta', type=str, help='Input FASTA file to be filtered (can be .fasta or .fasta.gz).')
    parser.add_argument('mapping_file', type=str, help='TSV file mapping UniProtKB proteins to Rhea IDs.')
    parser.add_argument('output_fasta', type=str, help='Output FASTA file with filtered records.')
    
    args = parser.parse_args()
    
    mapping = load_mapping(args.mapping_file)
    
    with (
        open(args.output_fasta, 'w') as out_handle,
        fileinput.hook_compressed(args.input_fasta, "r") as in_handle
        ):
        for record in filter_fasta(in_handle, mapping):
            SeqIO.write(record, out_handle, "fasta")

if __name__ == '__main__':
    main()
