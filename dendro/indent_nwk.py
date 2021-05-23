#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import stdin, stdout, stderr
import argparse
import re
from io import StringIO


# The following regex:
# - matches structural tokens: parentheses and comma;
# - includes white spaces before/after (stripped later);
# - matches delimiters to be jumped over: quotes and comments brackets.
PROTECT = {'[&&NHX:': ']', 
           '[&': ']', #'{': '}', # Beast nexus format
           '"': '"',
           "'": "'"}
IGNORE = re.compile('\033\[\d\d?(;\d+)*m')  # Shell escape codes (colored terminal output).
RE_STRUCT = re.compile(r'\s*(?:' + IGNORE.pattern + r'|,|\(|\)|' + r'|'.join(re.escape(s) for s in PROTECT) + r'|;)\s*')
INDENT = '  '


def go_in(nindent, indent):
    nindent += 1
    # return new indent size, and string to print before the rest of the text
    newstruct = '(\n' + nindent*indent
    return nindent, newstruct

def go_out(nindent, indent):
    nindent -= 1
    newstruct = '\n' + nindent*indent + ')'
    return nindent, newstruct

def go_same(nindent, indent):
    return nindent, ',\n' + nindent*indent

def go_quit(nindent, indent):
    return 0, ';\n'

STRUCT_DO = {',': go_same,
             '(': go_in,
             ')': go_out,
             ';': go_quit}


def main(treefiles, inplace=False, indent=INDENT):
    if not treefiles:
        treefiles = [stdin]
    for treefile in treefiles:
        if treefile is stdin:
            treetxt = stdin.read().rstrip(' \t\n\r')
        else:
            with open(treefile) as inputtree:
                treetxt = inputtree.read().rstrip(' \t\n\r')

        if inplace:
            out = StringIO()  # Buffer the output in-memory, in case of error.
        else:
            out = stdout
        
        try:
            # Split text by semantic units: nodes, commas, parentheses.
            nindent = 0
            #prev_pos = 0
            #for m in RE_STRUCT.finditer(treetxt):
            structmatch = RE_STRUCT.search(treetxt)
            while structmatch:
                struct = structmatch.group().strip()
                start,end = structmatch.start(), structmatch.end()
                if struct in STRUCT_DO:
                    nindent, newstruct = STRUCT_DO[struct](nindent, indent)

                    out.write(treetxt[:start].rstrip().lstrip() + newstruct)

                    treetxt = treetxt[end:]
                elif struct in PROTECT:
                    #out.write(treetxt[:end])
                    #treetxt = treetxt[structmatch.end():]
                    closing_char = PROTECT[struct]
                    try:
                        closing_pos = treetxt.index(closing_char, end)
                    except ValueError as err:
                        err.args = ('Missing closing structure for %r' % struct,)
                        raise
                    out.write(treetxt[:(closing_pos+1)])
                    treetxt = treetxt[(closing_pos+1):]
                elif IGNORE.match(struct):
                    out.write(treetxt[:end])
                    treetxt = treetxt[end:]
                else:
                    raise ValueError("Invalid structure character %r" % struct)

                structmatch = RE_STRUCT.search(treetxt)
            if inplace:
                inputtree.close()
                with open(inputtree.name, 'w') as realout:
                    realout.write(out.getvalue())
        
        except BrokenPipeError:  # IOError in Python2.7
            if out is stdout:
                print('indent_nwk:Caught BrokenPipeError to stdout.', file=stderr)
            else:
                raise
        finally:
            if out is not stdout:
                out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treefiles', nargs='*', help='[stdin]')
    parser.add_argument('-i', '--inplace', action='store_true',
                        help='Edit input file in place')
    
    main(**vars(parser.parse_args()))

