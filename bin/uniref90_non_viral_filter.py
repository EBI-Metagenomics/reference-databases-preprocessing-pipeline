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

from Bio import Entrez
import pyfastx
import taxoniq

# Global cache to store TaxID search results
taxid_cache = {}

UNWANTED_TAXIDS = ["Viruses", "unclassified entries"]

def is_unwanted_taxoniq(tax_id):
    """
    Check if taxa id is among UNWANTED_TAXIDS using faster taxoniq search in local, 
    indexed, compressed copy of the NCBI taxonomy database
    """
    taxon = taxoniq.Taxon(tax_id)
    if taxon.ranked_lineage:
        highest_taxon = taxon.ranked_lineage[-1]
        if highest_taxon.rank.name == "superkingdom":
            return highest_taxon.scientific_name in UNWANTED_TAXIDS
    raise ValueError("No information about the lineage in taxoniq")

def is_unwanted_entrez(taxid):
    """
    Check if taxa id is among UNWANTED_TAXIDS using slower approach with Entrez querying
    """
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    lineage = records[0]["Lineage"]
    highest_taxon = lineage.split("; ")[0]
    return highest_taxon in UNWANTED_TAXIDS

def is_unwanted(tax_id):
    """
    Wrapper function that checks cache before querying taxoniq or Entrez.
    Caches the result to avoid redundant queries.
    """
    if tax_id in taxid_cache:
        return taxid_cache[tax_id]

    is_unwanted = False
    try:
        is_unwanted = is_unwanted_taxoniq(tax_id)
    except (KeyError, ValueError):
        # If taxoniq fails, fallback to Entrez
        is_unwanted = is_unwanted_entrez(tax_id)

    # Cache the result
    taxid_cache[tax_id] = is_unwanted
    return is_unwanted

def filter_fasta(fasta, err_handle):
    """
    Remove proteins with UNWANTED_TAXIDS from FASTA using TaxID field.
    Returns a list of filtered records.
    """
    for seq in fasta:
        try:
            tax_id = seq.description.split("TaxID=")[1].split()[0]
            if not is_unwanted(tax_id):
                yield seq, tax_id
        except Exception as e:
            err_handle.write(f"Error while processing {tax_id}: {e}\n")

def processing_handle(input_fasta, output_prefix):
    """
    Filter the proteins in the input fasta and write output.
    Also write a mapping of protein id to taxid.
    """
    with (
        open("failed.txt", "w") as err_handle,
        open(f"{output_prefix}.filtered.fasta", 'w') as out_handle,
        open(f"{output_prefix}.protid2taxid", 'w') as mapping_handle
        ):
        mapping_handle.write("accession.version\ttaxid\n")

        fasta = pyfastx.Fasta(input_fasta)
        for seq, tax_id in filter_fasta(fasta, err_handle):
            out_handle.write(seq.raw)

            mapping_handle.write(f"{seq.name}\t{tax_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Remove viral proteins and proteins that taxonomically classified as 'metagenomes' from FASTA using TaxID field")
    parser.add_argument('input_fasta', type=str, help='Input FASTA file to be cleaned (can be .fasta or .fasta.gz)')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output FASTA file with filtered records and mapping of proteins to taxids')

    args = parser.parse_args()
    
    processing_handle(args.input_fasta, args.output_prefix)

if __name__ == '__main__':
    main()
