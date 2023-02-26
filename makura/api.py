#!/usr/bin/env python
from flask import Flask, jsonify, request


app = Flask(__name__)


@app.route("/summary")
def summary():
    from makura.assembly import AssemblySummary

    accessions = request.args.get("accessions")
    assembly_source = "refseq" if "GCF_" in accessions else "genbank"
    accessions_ls = accessions.split(",")
    asmsum = AssemblySummary(db_type=assembly_source)
    return asmsum.summary(accessions=accessions_ls)


if __name__ == "__main__":
    app.run(host="0.0.0.0")
