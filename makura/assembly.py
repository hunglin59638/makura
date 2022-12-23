#!/usr/bin/env python3
import argparse
import concurrent.futures
import hashlib
import itertools
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from time import sleep

import pandas as pd
import requests
from tqdm import tqdm

import makura
from makura import PKG_PATH, __version__
from makura.taxonparse import TaxonomyParser

AGLEA_TAXIDS = [
    2836,
    304574,
    29197,
    131213,
    3041,
    2825,
    304573,
    3027,
    39119,
    2864,
    3035,
    5747,
    38254,
    2830,
    131220,
    96475,
    2870,
    82162,
    569578,
    38410,
    2763,
    2833,
]


class AssemblySummary:

    _REFSEQ_GROUPS = [
        "archaea",
        "bacteria",
        "fungi",
        "invertebrate",
        "plant",
        "protozoa",
        "vertebrate_mammalian",
        "vertebrate_other",
        "viral",
    ]

    _GENBANK_GROUPS = [
        "archaea",
        "bacteria",
        "fungi",
        "invertebrate",
        "metagenomes",
        "other",
        "plant",
        "protozoa",
        "unknown",
        "vertebrate_mammalian",
        "vertebrate_other",
        "viral",
    ]

    _REFSEQ_CATEGORY = ["na", "reference genome", "representative genome"]

    _ASSEMBLY_LEVEL = ["Chromosome", "Complete Genome", "Contig", "Scaffold"]

    microbial_grops = ["archaea", "bacteria", "fungi", "viral", "protozoa", "aglea"]

    def __init__(self, assembly_summary=None, db_type="refseq"):
        self.db_type = db_type
        if self.db_type not in ("refseq", "genbank"):
            raise Exception("db_type must be refseq or genbank")
        self.assembly_summary = (
            PKG_PATH / f"assembly_summary_{db_type}.txt"
            if assembly_summary is None
            else Path(assembly_summary)
        )
        if not self.assembly_summary.is_file():
            self.update_assembly_summary()

        self.read_assembly_summary()

    def read_assembly_summary(self):
        self.assembly_summary_df = pd.read_csv(
            self.assembly_summary, sep="\t", low_memory=False
        )
        self.assembly_summary_df = self.assembly_summary_df.astype(
            {"taxid": "int64", "species_taxid": "int64"}
        )

    def update_assembly_summary(self):
        groups = (
            self._REFSEQ_GROUPS if self.db_type == "refseq" else self._GENBANK_GROUPS
        )

        rows = []
        for group in groups:
            url = f"https://ftp.ncbi.nlm.nih.gov/genomes/{self.db_type}/{group}/assembly_summary.txt"
            response = requests.get(url)
            for row in response.text.split("\n"):
                if row.startswith("#   See ") or not row:
                    continue
                elif row.startswith("#"):
                    cols = row.strip("# ").strip().split("\t")
                else:
                    row = f"{row}\t{group}".split("\t")
                    rows.append(row)
            print(f"{group} of assembly summary is fetched")

        cols.append("group")
        df = pd.DataFrame(rows, columns=cols)

        taxparser = TaxonomyParser()
        taxparser.update_taxonomy_database()
        taxids = [str(taxid) for taxid in df["taxid"].to_list() if taxid]

        taxid2lineage_dct = dict(
            [(taxid, taxparser.get_lineage(taxid)) for taxid in taxids]
        )

        all_rank_taxids = all_rank_taxids = list(
            set(list(itertools.chain.from_iterable(taxid2lineage_dct.values()))).union(
                set(taxids)
            )
        )
        taxid2name_dct = taxparser.get_taxid_translator(all_rank_taxids)
        taxid2rank_dct = taxparser.get_rank(all_rank_taxids)
        taxid_rows = []
        for taxid in taxids:
            lineage = taxid2lineage_dct[taxid]
            row = ";".join(
                [
                    "|".join([str(t), taxid2rank_dct[t], taxid2name_dct[t]])
                    for t in lineage
                ]
            )
            taxid_rows.append(row)
        taxid_df = pd.DataFrame(taxid_rows, columns=["taxid_lineage"])

        final_df = pd.concat([df, taxid_df], axis=1)
        final_df.to_csv(self.assembly_summary, index=False, sep="\t")

    def filter_accession_by_group(self, groups=[]):
        df = self.assembly_summary_df
        filter_df = df[df["group"].isin(groups)]
        accessions = filter_df["assembly_accession"].to_list()
        return accessions

    def download_by_group(self, out_dir, groups=[]):
        accessions = self.filter_accession_by_group(groups=groups)
        self.download_by_accession(out_dir, accessions)

    def filter_accession_by_taxid(self, taxids=[]):
        df = self.assembly_summary_df.dropna(subset=["taxid_lineage", "taxid"])
        accessions = [
            a
            for a in df.apply(
                lambda rec: rec["assembly_accession"]
                if set([str(t) for t in taxids])
                & set(
                    [value.split("|")[0] for value in rec["taxid_lineage"].split(";")]
                )
                else None,
                axis=1,
            ).to_list()
            if a
        ]
        return accessions

    def download_by_taxid(self, out_dir, taxids=[]):
        accessions = self.filter_accession_by_taxid(taxids=taxids)
        self.download_by_accession(
            out_dir=out_dir, accessions=accessions, use_rsync=True, parallel=4
        )

    def _download_job(self, ftp_path, out_dir, use_rsync=True):

        prefix = ftp_path.split("/")[-1]
        local_genome_file = Path(out_dir) / f"{prefix}_genomic.fna.gz"
        genome_url = f"{ftp_path}/{prefix}_genomic.fna.gz"
        md5checksums_url = f"{ftp_path}/md5checksums.txt"

        if use_rsync:
            genome_url = re.sub(r"^https", "rsync", genome_url)
            returncode = os.system(
                f"rsync --copy-links --times --quiet {genome_url} {out_dir}"
            )
            if returncode:
                print(f"CMD: rsync --copy-links --times --quiet {genome_url} {out_dir}")
                data = None
            else:
                data = local_genome_file

            with tempfile.NamedTemporaryFile("w", delete=False) as fp:
                tmp_file = Path(fp.name)
            md5checksums_url = re.sub(r"^https", "rsync", md5checksums_url)
            returncode = os.system(
                f"rsync --copy-links --times --quiet {md5checksums_url} {tmp_file}"
            )
            if returncode:
                md5checksums_content = None
            else:
                md5checksums_content = tmp_file.read_text()
                tmp_file.unlink()

        else:

            with requests.Session() as session:
                data = b""
                r = session.get(genome_url, stream=True)
                with open(local_genome_file, "wb") as f:
                    for chunk in r.raw.stream(4096, decode_content=False):
                        if chunk:
                            f.write(chunk)
                            data += chunk

                r = session.get(md5checksums_url)
                md5checksums_content = r.text

        original_md5 = ""
        if md5checksums_content:
            for md5_row in md5checksums_content.split("\n"):
                if md5_row:
                    md5_value, file = md5_row.split("  ")
                    if file == f"./{prefix}_genomic.fna.gz":
                        original_md5 = md5_value
                        break

        return (data, original_md5)

    def download_by_accession(self, out_dir, accessions=[], use_rsync=True, parallel=4):
        if isinstance(accessions, str):
            accessions = [accessions]
        out_dir = Path(out_dir)
        df = pd.read_csv(
            self.assembly_summary,
            sep="\t",
            usecols=["assembly_accession", "organism_name", "ftp_path"],
        )
        filter_df = df[
            (df["assembly_accession"].isin(accessions)) & (df["ftp_path"] != "na")
        ]

        def download_iter():
            for index, row in filter_df.iterrows():
                org_name = row["organism_name"]
                ftp_path = row["ftp_path"]
                yield org_name, ftp_path

        def _job(org_name, ftp_path):
            """
            returncode: return download status
            1: the genome file is existed, and skill it
            2: download task is success, and append the md5 of the file to md5checksums.txt
            3: download task is failed and remove the incomplete file
            """
            prefix = ftp_path.split("/")[-1]
            local_genome_file = out_dir / f"{prefix}_genomic.fna.gz"
            if local_genome_file.is_file():
                returncode = 1
            else:
                data, md5_value = self._download_job(ftp_path, out_dir, use_rsync)
                success = self.validate_md5(data, md5_value)
                if success:
                    file = f"./{local_genome_file.name}"
                    with open(out_dir / "md5checksums.txt", "a") as md5_h:
                        md5_h.write(f"{md5_value}\t{file}\n")
                    returncode = 2
                else:
                    local_genome_file.unlink(missing_ok=True)
                    returncode = 3

            return returncode

        with tqdm(total=filter_df.shape[0], desc="Download NCBI genomes") as pbar:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=parallel
            ) as executor:
                futures = {
                    executor.submit(
                        _job,
                        *args,
                    ): args
                    for args in download_iter()
                }

                results = {}
                count = 0
                for future in concurrent.futures.as_completed(futures):
                    org_name, ftp_path = futures[future]
                    pbar.set_description(f"{org_name} is downloaded")
                    returncode = future.result()
                    results[ftp_path] = returncode
                    pbar.update(1)

                    if returncode != 1:
                        count += 1
                    if count == 500:
                        pbar.set_description(f"sleep 30s...")
                        sleep(30)
                        count = 0

        success = True
        download_status_f = out_dir / "download_status.txt"
        with open(download_status_f, "w") as f:
            for ftp_path, returncode in results.items():
                f.write(f"{ftp_path}\t{returncode}\n")
                if returncode == 3:
                    success = False
        if not success:
            download_status_f.unlink()

        return results

    def filter_refseq_category(self, accessions=[], categories=[]):
        if not categories:
            categories = self._REFSEQ_CATEGORY

        df = self.assembly_summary_df
        filter_regex = "|".join([c.split(" ")[0] for c in categories])
        filter_df = df[df["refseq_category"].str.match(filter_regex)]

        if accessions:
            filter_df = filter_df[filter_df["assembly_accession"].isin(accessions)]

        return filter_df["assembly_accession"].to_list() if accessions else []

    def filter_assembly_level(self, accessions=[], levels=[]):
        if not levels:
            levels = self._ASSEMBLY_LEVELS

        df = self.assembly_summary_df
        filter_df = df[df["assembly_level"].isin(levels)]

        if accessions:
            filter_df = filter_df[filter_df["assembly_accession"].isin(accessions)]

        return filter_df["assembly_accession"].to_list() if accessions else []

    def select_microbiome(self):
        accessions = self.filter_accession_by_group(groups=self.microbial_grops)
        return accessions

    @staticmethod
    def sum_md5(data):
        if isinstance(data, Path):
            with open(data, "rb") as file_to_check:
                content = file_to_check.read()
        else:
            content = data
        md5_returned = hashlib.md5(content).hexdigest()

        return md5_returned

    @classmethod
    def validate_md5(cls, data, md5_value):
        """
        input a file path or string
        """
        if data is None or md5_value is None:
            return False

        md5_returned = cls.sum_md5(data)
        if md5_returned == md5_value:
            return True
        else:

            return False

    @staticmethod
    def uniq_md5checksums(md5checksums):
        md5checksums = Path(md5checksums)
        tmp_checksums = md5checksums.parent / f"{md5checksums.name}.tmp"
        with open(tmp_checksums, "w") as f:
            f.write(
                "\n".join(set([i for i in md5checksums.read_text().split("\n") if i]))
                + "\n"
            )
        tmp_checksums.rename(md5checksums)

    def check_genomes(self, md5checksums):
        """
        find which genomes are not in md5checksum.txt, and remove them.
        """
        md5checksums = Path(md5checksums)
        self.uniq_md5checksums(md5checksums)
        checked_files = set()
        with open(md5checksums, "r") as f:
            for row in f.readlines():
                row = row.strip("\n")
                if not row:
                    continue
                md5_value, genome_file = row.split("\t")
                genome_file = md5checksums.parent / Path(genome_file).name
                checked_files.add(genome_file.name)
        existed_file = set(
            [
                file.name
                for file in md5checksums.parent.iterdir()
                if file.name.endswith(".fna.gz")
            ]
        )
        need_removed = existed_file - checked_files

        [(md5checksums.parent / file).unlink() for file in need_removed]
        return list(need_removed)

    def summary(
        self, accessions=[], taxids=[], groups=[], assembly_level=[], refseq_category=[]
    ):
        accessions = list(
            set(self.filter_accession_by_group(taxids))
            | set(self.filter_accession_by_taxid(groups))
            | set(accessions)
        )
        df = self.assembly_summary_df
        filter_df = df[df["assembly_accession"].isin(accessions)]
        if assembly_level:
            filter_df = df[df["assembly_level"].isin(assembly_level)]
        if refseq_category:
            filter_df = df[df["refseq_category"].isin(refseq_category)]
        return filter_df.to_dict("records")


def main():
    def get_argument():
        class ExplicitDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
            def _get_help_string(self, action):
                if action.default in (None, False):
                    return action.help
                return super()._get_help_string(action)

        asm_levels = ["chromosome", "complete", "contig", "scaffold"]
        refseq_categories = ["reference", "representative", "na"]
        groups = [
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "metagenomes",
            "other",
            "plant",
            "protozoa",
            "unknown",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        ]

        def check_file(path):
            if not path:
                raise TypeError("Please input path")
            else:
                path = Path(path)
                if not path.exists():
                    raise argparse.ArgumentTypeError("No such as a file or directory")
                else:
                    return path
            raise TypeError("Please input path")

        def parse_dir(path):
            if not path:
                raise TypeError("Please input path")

            path = Path(path)
            path.mkdir(exist_ok=True)
            return path

        def split_comma(str_list):
            if str_list:
                return [i.strip() for i in str_list.split(",")]

        def check_levels(str_lvls):
            lvl_ls = split_comma(str_lvls)
            if set(lvl_ls) - set(asm_levels):
                raise ImportError(f"assembly-level must be {','.join(asm_levels)}")
            else:
                lvl_map = dict(
                    zip(
                        ["chromosome", "complete", "contig", "scaffold"],
                        ["Chromosome", "Complete Genome", "Contig", "Scaffold"],
                    )
                )
                return [lvl_map[lvl] for lvl in lvl_ls]

        def check_refseq_category(str_category):
            cat_ls = split_comma(str_category)
            if set(cat_ls) - set(refseq_categories):
                raise ImportError(
                    f"refseq-category must be {','.join(refseq_categories)}"
                )
            else:
                return cat_ls

        def limit_parellel(n):
            n = int(n)
            max_task = 8
            if n > max_task:
                raise TypeError(f"maximum number of task is {max_task}")
            else:
                return n

        parser = argparse.ArgumentParser(
            formatter_class=ExplicitDefaultsHelpFormatter,
            description=f"{makura.__name__.capitalize()}: Download NCBI genomes",
        )
        parser.add_argument(
            "--version",
            "-v",
            action="version",
            version=f"{makura.__name__} {__version__}",
        )

        subcmd = parser.add_subparsers(
            dest="subcmd",
            description="subcommands",
            metavar="SUBCOMMAND",
        )

        update = subcmd.add_parser(
            "update",
            help="update assembly summary and taxonomy information",
            formatter_class=ExplicitDefaultsHelpFormatter,
        )

        update.add_argument(
            "--assembly-source",
            default="refseq",
            choices=["refseq", "genbank"],
            help="select RefSeq(GCF_) or Genbank(GCA_) genomes",
        )

        download = subcmd.add_parser(
            "download",
            help="Download genomes from RefSeq or GenBank database",
            formatter_class=ExplicitDefaultsHelpFormatter,
        )
        download.add_argument(
            "--out_dir",
            "-o",
            help="output directory",
            type=parse_dir,
            default=Path().cwd(),
        )
        download.add_argument(
            "--update", "-U", help="update the assembly summary", action="store_true"
        )
        input_group = download.add_mutually_exclusive_group(required=True)
        input_group.add_argument(
            "--accessions",
            "-a",
            default=[],
            help="Download by assembly accessions (comma-separated)",
        )
        input_group.add_argument(
            "--accession-list",
            help="Download by assembly accession list",
            type=check_file,
        )
        input_group.add_argument(
            "--taxids",
            "-t",
            default=[],
            type=split_comma,
            help="Download by taxid (comma-separated)",
        )
        input_group.add_argument(
            "--taxid-list", help="Download by taxid list", type=check_file
        )
        input_group.add_argument(
            "--groups",
            "-g",
            default=[],
            type=split_comma,
            help=f"""
                Limit to genomes at one or more groups (comma-separated)\n
                {','.join(groups)}
                """,
        )
        input_group.add_argument(
            "--microbiome", "-m", help="Download microbiome", action="store"
        )

        download.add_argument(
            "--assembly-level",
            default=[],
            type=check_levels,
            help=f"""Limit to genomes at one or more assembly levels (comma-separated)
            {','.join(asm_levels)}
            """,
        )

        download.add_argument(
            "--refseq-category",
            default=[],
            type=check_refseq_category,
            help=f"""
                Limit to genomes at one or more refseq category (comma-separated)\n
                {','.join(refseq_categories)}
                """,
        )

        download.add_argument(
            "--assembly-source",
            default="refseq",
            choices=["refseq", "genbank"],
            help="select RefSeq(GCF_) or Genbank(GCA_) genomes",
        )

        download.add_argument(
            "--download-method",
            default="rsync",
            choices=["rsync", "https"],
            help="select the method for downloading genomes",
        )
        download.add_argument(
            "-p",
            "--parallel",
            dest="parallel",
            type=limit_parellel,
            default=4,
            help="Download genome with multiple processes (max: 8)",
        )
        summary = subcmd.add_parser(
            "summary",
            formatter_class=ExplicitDefaultsHelpFormatter,
            help="Print a data report with genome metadata from RefSeq or GenBank database in JSON format",
        )
        summary.add_argument(
            "--update", "-U", help="update the assembly summary", action="store_true"
        )
        summary_input_group = summary.add_mutually_exclusive_group(required=True)
        summary_input_group.add_argument(
            "--accessions",
            "-a",
            default=[],
            help="print records by assembly accessions (comma-separated)",
        )
        summary_input_group.add_argument(
            "--accession-list",
            help="print records by assembly accession list",
            type=check_file,
        )
        summary_input_group.add_argument(
            "--taxids",
            "-t",
            default=[],
            type=split_comma,
            help="print records by taxid (comma-separated)",
        )
        summary_input_group.add_argument(
            "--taxid-list", help="print records by taxid list", type=check_file
        )
        summary_input_group.add_argument(
            "--microbiome",
            "-m",
            help="print records of microbiome",
            action="store_true",
        )

        summary.add_argument(
            "--assembly-level",
            default=[],
            type=check_levels,
            help=f"""Limit to genomes at one or more assembly levels (comma-separated)
            {','.join(asm_levels)}
            """,
        )

        summary.add_argument(
            "--refseq-category",
            default=[],
            type=check_refseq_category,
            help=f"""
                Limit to genomes at one or more refseq category (comma-separated)\n
                {','.join(refseq_categories)}
                """,
        )

        summary.add_argument(
            "--assembly-source",
            default="refseq",
            choices=["refseq", "genbank"],
            help="select RefSeq(GCF_) or Genbank(GCA_) genomes",
        )
        args = parser.parse_args()
        return args

    args = get_argument()

    assembly_source = args.assembly_source
    asmsum = AssemblySummary(db_type=assembly_source)

    if args.subcmd == "update":
        asmsum.update_assembly_summary()
        sys.exit(0)

    refseq_category = args.refseq_category
    assembly_level = args.assembly_level
    if args.accessions or args.accession_list:
        if args.accessions:
            accessions = args.accessions
        else:
            with open(args.accession_list, "r") as f:
                accessions = [
                    row.strip("\n ") for row in f.readlines() if row.strip("\n ")
                ]

    elif args.taxids or args.taxid_list:
        if args.taxids:
            taxids = args.taxids
        else:
            with open(args.taxid_list, "r") as f:
                taxids = [row.strip("\n ") for row in f.readlines() if row.strip("\n ")]

        accessions = asmsum.filter_accession_by_taxid(taxids=taxids)

    elif args.groups:
        accessions = asmsum.filter_accession_by_group(groups=args.groups)

    elif args.micrbiome:
        accessions = asmsum.select_microbiome()

    else:
        accessions = []

    if refseq_category:
        accessions = asmsum.filter_refseq_category(
            accessions=accessions, categories=refseq_category
        )

    if assembly_level:
        accessions = asmsum.filter_assembly_level(
            accessions=accessions, levels=assembly_level
        )

    if args.subcmd == "download":
        out_dir = args.out_dir
        parallel = args.parallel
        use_rsync = True if args.download_method == "rsync" else False

        md5checksums_f = out_dir / "md5checksums.txt"
        if md5checksums_f.is_file():
            print("uniq md5checksums.txt and remove genomes not in md5checksums.txt")
            removed = asmsum.check_genomes(md5checksums_f)
            print(f"Remove {len(removed)} genomes")
        asmsum.download_by_accession(
            out_dir, accessions=accessions, use_rsync=use_rsync, parallel=parallel
        )
    elif args.subcmd == "summary":
        records = asmsum.summary(accessions=accessions)
        print(json.dumps(records, indent=4))
