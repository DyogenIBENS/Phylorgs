#!/usr/bin/env python3


from sys import stdout
import argparse
import requests
import logging
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter("%(levelname)s:%(filename)s:%(message)s"))
logger.addHandler(ch)


URL = "http://www.ensembl.org/biomart/"
#ARCHIVE = "http://{}.archive.ensembl.org/biomart/"
ARCHIVES = {
            103: "http://feb2021.archive.ensembl.org/biomart/",
            102: "http://nov2020.archive.ensembl.org/biomart/",
            101: "http://aug2020.archive.ensembl.org/biomart/",
            100: "http://apr2020.archive.ensembl.org/biomart/",
            99: "http://jan2020.archive.ensembl.org/biomart/",
            98: "http://sep2019.archive.ensembl.org/biomart/",
            97: "http://jul2019.archive.ensembl.org/biomart/",
            96: "http://apr2019.archive.ensembl.org/biomart/",
            95: "http://jan2019.archive.ensembl.org/biomart/",
            94: "http://oct2018.archive.ensembl.org/biomart/",
            93: "http://jul2018.archive.ensembl.org/biomart/",
            87: "http://dec2016.archive.ensembl.org/biomart/",
            86: "http://oct2016.archive.ensembl.org/biomart/",
            85: "http://jul2016.archive.ensembl.org/biomart/",
            84: "http://mar2016.archive.ensembl.org/biomart/"}

FORM="martservice?query="

# example query

QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "{formatter}" header = "{header}" uniqueRows = "{uniqueRows}" count = "{count}" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        {filters}
        {attributes}
    </Dataset>
</Query>"""

FILTER = '<Filter name = "{}" value = "{}"/>'
ATTRIBUTE = '<Attribute name = "{}" />'


def query_fromfile(queryfile, outfile='-', ensembl_version=None):
    with open(queryfile) as qf:
        query = qf.read()
    do_query(query, outfile, ensembl_version)


def query_fromargs(outfile='-', ensembl_version=None, formatter='TSV', header=0,
                    uniqueRows=0, count='', dataset='hsapiens_gene_ensembl',
                    filters=None, attributes=None):
    if filters:
        filter_lines = '\n'.join(FILTER.format(*f.split('=')) for f in filters)
    else:
        filter_lines = ''

    if attributes:
        attr_lines = '\n'.join(ATTRIBUTE.format(a) for a in attributes)
    else:
        attr_lines = ''

    query = QUERY.format(formatter=formatter, header=header, uniqueRows=uniqueRows,
                         count=count, dataset=dataset, filters=filter_lines,
                         attributes=attr_lines)
    
    do_query(query, outfile, ensembl_version=None)


def do_query(query, outfile='-', ensembl_version=None):
    url = ARCHIVES.get(ensembl_version, URL) + FORM + query
    response = requests.get(url)
    response.raise_for_status()
    
    if response.text:
        out = stdout if outfile == '-' else open(outfile, 'w')
        out.write(response.text)
        if outfile != '-': out.close()
    else:
        logger.error("No content. Status code: %d", response.status_code)


def main(**kwargs):
    if kwargs.get('queryfile'):
        for kw in ('formatter', 'header', 'uniqueRows', 'count', 'dataset', 
                   'filters', 'attributes'):
            kwargs.pop(kw)

        query_fromfile(**kwargs)
    
    else:
        kwargs.pop('queryfile')
        query_fromargs(**kwargs)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('queryfile', nargs='?')
    parser.add_argument('-o', '--outfile', default='-')
    parser.add_argument('-e', '--ensembl-version')
    qgroup = parser.add_argument_group('Query arguments', description='When no queryfile is given.')
    qgroup.add_argument('-F', '--formatter', choices=['TSV', 'FASTA'], default='TSV')
    qgroup.add_argument('-H', '--header', type=int, choices=[0, 1], default=0)
    qgroup.add_argument('-u', '--uniqueRows', type=int, choices=[0, 1], default=0)
    qgroup.add_argument('-c', '--count', default='')
    qgroup.add_argument('-d', '--dataset', default='hsapiens_gene_ensembl')
    qgroup.add_argument('-f', '--filters', action='append',
                        help='syntax: "<name>=<value>"')
    qgroup.add_argument('-a', '--attributes', nargs='+')

    args = parser.parse_args()
    main(**vars(args))

