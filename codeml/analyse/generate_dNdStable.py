#!/usr/bin/env python

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""

from __future__ import print_function

import sys
import os.path
import re
import argparse
import ete3
import LibsDyogen.utils.myPhylTree as PhylTree

from select_leaves_from_specieslist import convert_gene2species


ANCGENE2SP = re.compile(r'([A-Z][A-Za-z_.-]+)ENSGT')


def print_if_verbose(*args, **kwargs):
    pass

def showtree(fulltree, ages):
    pass

def def_showtree(show=False):
#    fulltree.show()
    if show:
        def showtree(fulltree, ages):
            ages_dict = {row[0]: row[1:] for row in ages}
            for n in fulltree.traverse():
                info = ages_dict.get(n.name, [0, 'leaf'])
                n.add_feature('age', info[0])
                n.add_feature('type', info[1])
            fulltree.show()
    else:
        def showtree(fulltree, ages):
            pass
    return showtree


def branch2nb(mlc, replace_nwk='.mlc'):
    """mlc: filehandle"""
    print_if_verbose()
    regex = re.compile(r'^(.*);$')
    regex_lnL = re.compile(r'^lnL\(')
    # get the line listing all branches
    line = mlc.next().rstrip()
    while not regex_lnL.match(line):
        line = mlc.next().rstrip()
    branches_line = mlc.next().rstrip()
    branches = [b.split('..') for b in branches_line.split()]
    # get translation by looking at the newick tree lines.
    # I could also get it by the lines #1: ENS...
    line = mlc.next().rstrip()
    while not regex.match(line):
        line = mlc.next().rstrip()
    tree_nbs = ete3.Tree(line)

    line = mlc.next().rstrip()
    while not regex.match(line):
        line = mlc.next().rstrip()
    tree_ids = ete3.Tree(line)

    id2nb = dict(zip(tree_ids.get_leaf_names(), tree_nbs.get_leaf_names()))
    nb2id = dict(zip(tree_nbs.get_leaf_names(), tree_ids.get_leaf_names()))

    # get internal nodes nb-to-id conversion (need the tree with internal node
    # annotations)
    nwkfile = mlc.name.replace(replace_nwk, '.nwk')
    try:
        fulltree = ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            print("\nNewickError: Malformed newick tree structure in %r" \
                    % nwkfile,
                    file=sys.stderr)
        else:
            print("\nNewickError: Unexisting tree file %r" % nwkfile,
                    file=sys.stderr)
        sys.exit(1)


    for leafnb, leafid in zip(tree_nbs.get_leaves(), fulltree.get_leaves()):
        leafid.add_feature('nb', leafnb.name)

    # branches follow the order of the newick string.
    while branches:
        base, tip = branches.pop()
        print_if_verbose("%-8s" % (base + '..' + tip), end=' ')
        try:
            tip_id = nb2id[tip]
            try:
                base_node = fulltree.search_nodes(name=tip_id)[0].up
                while len(base_node.children) == 1:
                    base_node = base_node.up
                base_id = base_node.name
                print_if_verbose("%s -> %s  " % (base_id, tip_id), end=" ")
                nb2id[base] = base_id
                id2nb[base_id] = base
                # Add number in the fulltree:
                base_node.add_feature('nb', base)
                # Also update the tree_nbs
                base_nb_node = tree_nbs.search_nodes(name=tip)[0].up
                base_nb_node.name = base
                print_if_verbose('ok')
            except AttributeError:
                print('root (%r: %r) cannot have ancestors' % (tip, tip_id))
                pass
        except KeyError as e:
            # Not found now, put it back in the queue for later
            branches.insert(0, (base, tip))
            print_if_verbose('KeyError')

    print_if_verbose(tree_nbs.get_ascii())
    return id2nb, nb2id, tree_nbs, fulltree


def get_dNdS(mlc):
    """Parse table of dN/dS from codeml output file.
    
    mlc: filehandle
    """
    #reg_dNdS = re.compile(r'dN & dS for each branch')
    reg_dNdS = re.compile(r'\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS')
    line = mlc.next().rstrip()
    while not reg_dNdS.match(line):
        line = mlc.next().rstrip()
    assert mlc.next().rstrip() == ''  # skip blank line
    dNdS = {}
    line = mlc.next().rstrip()
    while line != '':
        #print line
        fields = line.split()
        dNdS[fields[0]] = [float(x) for x in fields[1:]]
        line = mlc.next().rstrip()
    return dNdS
    

def set_dS_fulltree(fulltree, id2nb, dNdS):
    for node in fulltree.get_descendants('postorder'): # exclude root
        if len(node.children) != 1:
            parent = node.up
            intermediates = []
            while parent.up and len(parent.children) == 1: # under root, up is None
                intermediates.append(parent)
                parent = parent.up

            if parent.is_root():
                # then this is a root with a single children: cannot have dS
                dS_tot = 0
            else:
                dS_tot = dNdS[id2nb[parent.name] + '..' + id2nb[node.name]][5]
    
            if intermediates: # this if is not necessary, but maybe for efficiency
                dist_tot = sum(n.dist for n in intermediates + [node])

                for inter_node in intermediates + [node]:
                    inter_node.add_feature('dS', dS_tot * inter_node.dist / dist_tot)
            else:
                node.add_feature('dS', dS_tot)


def rm_erroneous_ancestors(fulltree, phyltree):
    """Some ancestors with only one child are present in the species phylogeny.
    They must be remomved when they have an age of zero"""
    for node in fulltree.iter_descendants():
        if not node.is_leaf():
            try:
                taxon = ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % node.name)
            if phyltree.ages[taxon] == 0 and hasattr(node, 'reinserted'):
                node.delete(prevent_nondicotomic=False, preserve_branch_length=True)



def sum_average_dNdS(dNdS, nb2id, tree_nbs):
    """Compute the average dS depth of each node (between this node and
    leaves).
    
    Probably the roughest method. Use rather `bound_average_dS`"""
    dS, dN = {}, {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        if n.is_leaf():
            dS[nb2id[node]] = 0
            dN[nb2id[node]] = 0
            dS[node] = 0
            dN[node] = 0
        else:

            dS_branch = dNdS[node + '..' + c.name][5]
            dN_branch = dNdS[node + '..' + c.name][4]
            dS_children = [dS[nb2id[c.name]] + dS_branch for c in n.children]
            dN_children = [dN[nb2id[c.name]] + dN_branch for c in n.children]
            #dS_children = [dS[c.name]] for c in n.children]
            #dN_children = [dN[c.name]] for c in n.children]
            average_dS_ch = sum(dS_children) / len(dS_children)
            average_dN_ch = sum(dN_children) / len(dN_children)

            dS[nb2id[node]] = average_dS_ch
            dN[nb2id[node]] = average_dN_ch
            dS[node] = average_dS_ch
            dN[node] = average_dN_ch

    return dN, dS


def bound_average_dS(dNdS, nb2id, tree_nbs, phyltree):
    """calibrate tree using speciation nodes."""
    ages = []
    subtree = {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        scname = nb2id[node] # scientific name
        print_if_verbose("%2s. %s:" % (node, scname), end=' ')
        if n.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[node] = {'taxon': taxon, 'dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[node] = {'taxon': taxon}
            # determine if speciation or not
            children_taxa = set(subtree[c.name]['taxon'] for c in n.children)
            # compute and store average dS to next spe, for each child
            children_dS = [subtree[c.name]['dS'] + \
                            dNdS[node + '..' + c.name][5] \
                            for c in n.children]
            if len(children_taxa) == 1:
                # It is a duplication: average children_dS
                subtree[node]['dS'] = sum(children_dS) / len(children_dS)
                print_if_verbose("Dupl; dS=%s" % subtree[node]['dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[node]['age'] = subtree[n.children[0].name]['age']
            else:
                print_if_verbose("Spe")
                # speciation: store the age of this taxon
                node_age = subtree[node]['age'] = phyltree.ages[taxon]
                ages.append([scname, str(node_age), "spe"])
                # compute average dS to each next speciation.
                # update age of each posterior dup
                if len(n.children) > 1:
                    pass
                for c in n.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[c.name]['age']
                    dS_genebranch = dNdS[node + '..' + c.name][5]
                    scaling_dS = subtree[c.name]['dS'] + dS_genebranch
                    print_if_verbose("    climb up to next speciation: " \
                                      "scaling_dS=%s; br_len=%s" % \
                                        (scaling_dS, branch_length))
                    nextnodes = [c]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['dS']
                        print_if_verbose("    - %2s. %s: dS=%s" % \
                                            (nextnode.name,
                                             nb2id[nextnode.name],
                                             nextnode_dS))
                        if nextnode_dS > 0:
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nb2id[nextnode.name], str(age), "dup"])
                            nextnodes.extend(nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[node]['dS'] = 0
    return ages


def bound_average_dS(dNdS, id2nb, fulltree, phyltree):
    ages = []
    subtree = {}
    for node in fulltree.traverse('postorder'):
        scname = node.name
        print_if_verbose("* %3s. %s:" % (id2nb.get(scname), scname), end=' ')
        if node.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[scname] = {'taxon': taxon, 'tmp_dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon}
            #for child in node.children:
            #    child_name = child.name
            #    child_taxon = subtree[child_name]
            children_taxa = set((subtree[ch.name]['taxon'] for ch in node.children))

            if len(children_taxa & set((taxon,))) == 1:
                # it is a duplication:
                children_dS = [ch.dS + subtree[ch.name]['tmp_dS'] for ch in node.children]
                subtree[scname]['tmp_dS'] = sum(children_dS) / len(children_dS)

                print_if_verbose("Dupl; dS=%s" % subtree[scname]['tmp_dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[scname]['age'] = subtree[node.children[0].name]['age']

            else:
                # it is a speciation.
                print_if_verbose("Spe")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = phyltree.ages[taxon]
                ages.append([scname, str(node_age), "spe"])
                # climb up tree and assign an age to each duplication
                for ch in node.children:
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    # NOTE: could actually have used the phyltree.branches
                    branch_length = node_age - subtree[ch.name]['age']
                    dS_genebranch = ch.dS
                    scaling_dS = subtree[ch.name]['tmp_dS'] + dS_genebranch
                    print_if_verbose("    climb up to next speciation: " \
                                      "scaling_dS=%s; br_len=%s" % \
                                        (scaling_dS, branch_length))
                    nextnodes = [ch]
                    while nextnodes:
                        nextnode = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['tmp_dS']
                        print_if_verbose("    - %2s. %s: dS=%s" % \
                                            (id2nb.get(nextnode.name),
                                             nextnode.name,
                                             nextnode_dS))
                        if nextnode_dS > 0:
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nextnode.name, str(age), "dup"])
                            nextnodes.extend(nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[scname]['tmp_dS'] = 0
    return ages


def bound_average_dS_2(dNdS, id2nb, fulltree, phyltree):
    """each dS of a duplication is normalized by the D' + dS, where:
        - dS: average of branch length leading to next speciation, in dS
        - D': branch length from the previous speciation to this duplication"""
    ages = []
    subtree = {}
    for node in fulltree.traverse('postorder'):
        scname = node.name
        print_if_verbose("* %3s. %s:" % (id2nb.get(scname), scname), end=' ')
        if node.is_leaf():
            taxon = convert_gene2species(scname)
            subtree[scname] = {'taxon': taxon, 'tmp_dS': 0, 'age': 0} #, 'next_age': 0}
            #ages[scname] = 0
            print_if_verbose("Leaf")
        else:
            try:
                taxon = ANCGENE2SP.match(scname).group(1).replace('.', ' ')
            except AttributeError:
                raise RuntimeError("Can not match species name in %r" % scname)

            subtree[scname] = {'taxon': taxon}
            #for child in node.children:
            #    child_name = child.name
            #    child_taxon = subtree[child_name]
            children_taxa = set((subtree[ch.name]['taxon'] for ch in node.children))

            if len(children_taxa & set((taxon,))) == 1:
                # it is a duplication:
                children_dS = [ch.dS + subtree[ch.name]['tmp_dS'] for ch in node.children]
                subtree[scname]['tmp_dS'] = sum(children_dS) / len(children_dS)

                print_if_verbose("Dupl; dS=%s" % subtree[scname]['tmp_dS'])
                # store the age of the next speciation event.
                # Since it is a duplication, this should be the same for
                # children[0] and children[1].
                # This 'age' will later be modified (rescaled according to dS)
                subtree[scname]['age'] = subtree[node.children[0].name]['age']

            else:
                # it is a speciation.
                print_if_verbose("Spe")
                # store the age of this taxon
                node_age = subtree[scname]['age'] = phyltree.ages[taxon]
                ages.append([scname, str(node_age), "spe"])
                # climb up tree and assign an age to each duplication
                for ch in node.children:
                    ### This is where version _2 is different
                    # walk descendants until speciation.
                    # need to get the age of next speciation and compute the time
                    # between the two speciation.
                    branch_length = node_age - subtree[ch.name]['age']
                    print_if_verbose("    climb up to next speciation: " \
                                      "br_len=%s" % branch_length)
                    nextnodes = [(ch, ch.dS)]
                    while nextnodes:
                        nextnode, next_path_dS = nextnodes.pop(0)
                        nextnode_dS = subtree[nextnode.name]['tmp_dS']
                        print_if_verbose("    - %2s. %s: dS(to speciation)=%s"\
                                         "; dS(from speciation)=%s" % \
                                            (id2nb.get(nextnode.name),
                                             nextnode.name,
                                             nextnode_dS, next_path_dS))
                        if nextnode_dS > 0:
                            scaling_dS = next_path_dS + nextnode_dS
                            age = node_age - (1 - nextnode_dS / scaling_dS) * branch_length
                            #ages[nextnode.name] = age
                            ages.append([nextnode.name, str(age), "dup"])
                            nextnodes.extend((nch, nch.dS + next_path_dS) \
                                                for nch in nextnode.children)
                        subtree.pop(nextnode.name)

                # then reset dS to zero
                subtree[scname]['tmp_dS'] = 0
    return ages


def save_ages(ages, opened_outfile):
    for row in ages:
        opened_outfile.write("\t".join(row) + "\n")


def process_mlc(mlcfile, phyltree, replace_nwk='.mlc'):
    with open(mlcfile) as mlc:
        id2nb, nb2id, tree_nbs, fulltree = branch2nb(mlc, replace_nwk)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    rm_erroneous_ancestors(fulltree, phyltree)
    set_dS_fulltree(fulltree, id2nb, dNdS)
    ages = bound_average_dS(dNdS, id2nb, fulltree, phyltree)
    showtree(fulltree, ages)
    return ages

def process_mlc2(mlcfile, phyltree, replace_nwk='.mlc'):
    with open(mlcfile) as mlc:
        id2nb, nb2id, tree_nbs, fulltree = branch2nb(mlc, replace_nwk)
        dNdS = get_dNdS(mlc)

    #dN, dS = sum_average_dNdS(dNdS, nb2id, tree_nbs)
    rm_erroneous_ancestors(fulltree, phyltree)
    set_dS_fulltree(fulltree, id2nb, dNdS)
    ages = bound_average_dS_2(dNdS, id2nb, fulltree, phyltree)
    showtree(fulltree, ages)
    return ages


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outfile')
    parser.add_argument('mlcfiles', nargs='+')
    parser.add_argument('--method2', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='print progression along tree')
    parser.add_argument('--show', action='store_true',
                        help='start the interactive ete3 tree browser')
    parser.add_argument('-r', '--replace-nwk', default='.mlc',
                        help='string to be replaced by .nwk to find the tree'\
                               ' file [%(default)s]')
    args = parser.parse_args()
    outfile = args.outfile
    mlcfiles = args.mlcfiles
    nb_mlc = len(mlcfiles)
    
    # "compile" some functions to avoid redondant tests ("if verbose: ...")
    if args.verbose:
        def print_if_verbose(*args, **kwargs):
            print(*args, **kwargs)
    if args.method2:
        process_mlc = process_mlc2

    showtree = def_showtree(args.show)

    phyltree = PhylTree.PhylogeneticTree("/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data85/PhylTree.Ensembl.85.conf")
    
    with open(outfile, 'a') as out:
        for i, mlcfile in enumerate(mlcfiles, start=1):
            percentage = float(i) / nb_mlc * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_mlc, percentage, mlcfile),
                    file=sys.stderr, end=' ')
            try:
                ages = process_mlc(mlcfile, phyltree, args.replace_nwk)
                save_ages(ages, out)
            except:
                print()
                raise
    print()
    
