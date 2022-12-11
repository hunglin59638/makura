#!/usr/bin/env python3
from pathlib import Path

from ete3 import NCBITaxa

TAXDUMP_FILE = None


class TaxonomyParser(NCBITaxa):
    def __init__(self, dbfile=None, taxdump_file=None):
        super().__init__(dbfile, taxdump_file)

    def convert_name_to_taxid(self, names):
        pass

    def get_lineage_translator(self, taxids):
        return super().get_lineage_translator(taxids)

    def get_lineage(self, taxid):
        try:
            return super().get_lineage(taxid)
        except ValueError as e:
            return []
