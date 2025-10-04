from collections import defaultdict
from statistics import median
from eris.io import SeqFile
from eris.graph import Edge, Graph
from eris.seq import Record
from eris.scan import ScannerResult


def get_end_node(contigs: list[str], outdegree: dict, indegree: dict) -> str:
    """
    Search a node/contig that has no outging edge to IS identified contigs.

    A contig is considered an end if it has no outgoing edges to other IS contigs
    (checked up to 2 level), but still has  at least one incoming edge.
    If none is found, the graph is assumed circular and the first contig is returned.

    :params contigs: A list of contigs of an identified IS
    :param outdegree: Dictionary mapping each node to its outgoing edges
    :param indegree: Dictionary mapping each node to its incoming edges.
    :return: the end contig ID
    """
    # check whether a node has outgoing edge to the same IS node
    # False, if the node has no outgiong edges to IS contigs,
    # at least two level searching.

    to_is = {n: False for n in contigs}

    for node in contigs:
        edges1 = outdegree[node]
        for edge1 in edges1:
            to = edge1.to
            if to in contigs:
                to_is[node] = True
            else:
                # if IS contig is not found in the first neighbor,
                # check one level up, if that non IS identified contig
                # adjacent with contig that has identified IS
                edges2 = outdegree[to]
                for edge2 in edges2:
                    to = edge2.to
                    if to in contigs:
                        to_is[node] = True

    for k, v in to_is.items():
        # make sure that the end has indegree node
        if v is False:
            if indegree.get(k) is not None:
                return k

    # need to think if no above scenario are met
    print("Looks like the IS contigs graph is circular")
    return contigs[0] # just pick the first node


def get_start_node(contigs: list, indegree: dict) -> str:
    """
    Search a node/contig that has no outging edge to IS identified contigs.

    A contig is considered a start if it has no incoming edges from other IS contigs 
    (checked up to 2 levels). If no such node is found, the graph is assumed circular and the first 
    contig is returned.

    :params contigs: A list of contigs of an identified IS
    :param indegree: Dictionary mapping each node to its incoming edges.
    :return: the start contig ID
    """

    fr_is = {node: False for node in contigs}
    for node in contigs:

        # if no indegree, just start from here
        if indegree.get(node) is None:
            return node

        # look whether a node has an edge from IS identified node/contig
        edges1 = indegree.get(node)
        for edge1 in edges1:
            if edge1.fr in contigs:
                fr_is[node] = True
            else:
                # search one contig/node before it
                edges2 = indegree.get(edge1.fr)
                if edges2 is not None:
                    for edge2 in edges2:
                        if edge2.fr in contigs:
                            fr_is[node] = True

    for k, v in fr_is.items():
        if v is False:
            return k

    # need to think if no above scenario are met
    print("Looks like the IS contigs graph is circular")
    return contigs[0] # just pick the first node

def get_indegree(graph: Graph):
    """
    Create dictionary to stored indegree edges for each node

    :params graph: Graph object from eris
    """
    indegree = defaultdict(set)
    for _, edges in graph.adj.items():
        for edge in edges:
            indegree[edge.to].add(edge)
    return indegree

# find a shortest path from start to end contigs
def find_path(graph: Graph, start: str, end: str, depth_records: dict,
              is_with_contigs: dict,
              is_name: str, median_depth: float) -> tuple[list[str], int]:
    """
    Using BFS implementation, find a shortest path from start to end node,
    also quantify the copies of the elements.

    :params graph: Graph object from eris
    :params start: a start node of the IS contigs
    :params end: an end node of the IS contigs
    :params depth_record: record for the depth of each node in the genome graph
    :params is_with_contigs: a key of IS name with value of a list of the corresponding contigs in the graph
    :params is_name: the IS name
    :params median_depth: the median depth of overal nodes in the genome graph
    :return: the give IS path in the genome graph and its copy number prediction
    """

    path = []
    marked = {node: False for node in graph._node_ids}

    queue = [start]

    # if the path traver node that higher n times more than median,
    # the predicted copies then is approx n.
    copies = 1

    check_contigs = is_with_contigs[is_name]
    check_contigs_din = check_contigs.copy() # can shrink without affect the original

    while queue:
        v = queue.pop(0)
        if not marked[v]:
            path.append(v)

            # check if we traversing the node that higher than the median,
            # that likely is the collapse repeat.
            depth = depth_records[v]
            copies = max(round(depth/median_depth), copies)

            if v in check_contigs_din:
                check_contigs_din.remove(v)
            marked[v] = True
            if v == end: # if already find the end node, just return it
                return path, copies


        to_traverse = None
        deepest = 0
        for edge in graph.adj.get(v, []):

            # I want a node to be stored into queue
            # if this node is in the IS check contigs list
            neighbor = edge.to
            if neighbor in check_contigs:
                if depth_records[neighbor] > deepest:
                    to_traverse = neighbor

            # or this node has outgoing edegs to the IS node(s)
            else:
                for neighbor_edge in graph.adj.get(neighbor, []):
                    if neighbor_edge.to in check_contigs:
                        if depth_records[neighbor] > deepest:
                            to_traverse = neighbor

        # if it finds more than one path, select the one that has the highest depth
        if to_traverse is not None:
            queue.append(to_traverse)


    # if it finds no outgoing edges from last node, and no queue left but
    # not all IS nodes have been visited,
    # try to search again from the new start node
    if len(queue) == 0 and len(check_contigs_din) > 0:
        indegree = get_indegree(graph)
        start_node = get_start_node(check_contigs_din, indegree)
        queue = [start_node]

    return path, copies


def collapseis(eris_result: ScannerResult,
             genome_graph: SeqFile) -> dict[str: tuple[list, int]]:

    """
    Collapse and create a path from a fragmented IS in the given genome graph.
    Return the predicted path of the collapsed sequence and the copies of the sequence
    """

    result = defaultdict(tuple)

    # organise the identified IS with their corresponding contigs
    is_with_contigs = defaultdict(list)
    for entry in eris_result:
        if entry.kind == "mobile_element":
            is_name  = entry.qualifiers[0].value.split("_")[0]
            is_with_contigs[is_name].append(entry.location.parent_id)

    # Load all edges representation of the de brujin graph, assembly genome result
    # and build eris graph object
    gemome_edges = []
    node_depth = defaultdict(int)

    for entry in genome_graph:
        if isinstance(entry, Edge):
            gemome_edges.append(entry)
        elif isinstance(entry, Record):
            node_depth[entry.id] = entry.qualifiers[0].value # assuming depth info at 0 qualifiers or first field of the tag
    genome_graph.close()

    # create directed graph
    genome_graph = Graph(*gemome_edges, directed=True)
    # .adj represent only outgoing edges,
    outdegree_graph = genome_graph.adj

    # so i create indegree edges
    indegree_graph = get_indegree(genome_graph)

    # calculate the median
    median_depth = median(list(node_depth.values()))

    for is_name, is_contigs in is_with_contigs.items():

        # get the boundaries of the IS in the genome graph
        start = get_start_node(is_contigs, indegree_graph)
        end = get_end_node(is_contigs, outdegree_graph, indegree_graph)

        # find the path from start to end for an IS
        # also quantify the copies
        path_copies = find_path(genome_graph, start, end,
                                node_depth, is_with_contigs, is_name, median_depth)

        result[is_name] = path_copies

    return result
