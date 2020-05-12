from graphviz import render
from graphviz import Source
from Bio import SeqIO
from time import time


class Edge:
    def __init__(self, dna, coverage, node_from, node_to, weight=1):
        self.dna = dna
        self.node_from = node_from
        self.node_to = node_to
        self.coverage = coverage
        self.weight = weight
        self.deleted = False

    def __repr__(self):
        return 'Edge(dna=' + self.dna + ', coverage=' + str(self.coverage) + ')'


class Node:
    def __init__(self, dna, coverage=1):
        self.dna = dna
        self.coverage = coverage
        self.amount_in = 0
        self.edges_in = {}
        self.amount_out = 0
        self.edges_out = {}
        self.deleted = False

    def __repr__(self):
        return '\nNode(dna=' + self.dna + ', coverage=' + str(self.coverage) + \
               ',\n\tedges_in=' + str(self.edges_in) + ',\n\tedges_out=' + str(self.edges_out) + ')\n'


# функция, которая по строке возвращает ее же и комплиментарно-реверснутую ей, учитывая вырожденности.
def straight_and_reverse(dna):
    dna = dna.upper()
    straight = dna.replace('B', 'N'). \
        replace('D', 'N'). \
        replace('H', 'N'). \
        replace('V', 'N')
    reverse = ''
    reverse_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                    'R': 'Y', 'K': 'M', 'S': 'W', 'Y': 'R', 'M': 'K', 'W': 'S',
                    'B': 'N', 'D': 'N', 'H': 'N', 'V': 'N',
                    'N': 'N', '-': '-'}
    for s in dna:
        try:
            reverse += reverse_dict[s]
        except KeyError:
            reverse += 'N'
    return straight, reverse[::-1]


def visualize(name):
    render('dot', 'png', name)
    # To render an existing file in a notebook
    Source.from_file(name)


class Graph:
    def __init__(self):
        self.indexes = {}
        self.indexes_counter = 0
        self.nodes = {}

    def add_node(self, dna, coverage=1):
        straight, reverse = straight_and_reverse(dna)

        def add_seq(seq, seq_cov):
            if seq in self.nodes:
                self.nodes[seq].coverage += seq_cov
            else:
                self.indexes.update({seq: self.indexes_counter})
                self.indexes_counter += 1
                self.nodes.update({seq: Node(seq, seq_cov)})

        add_seq(straight, coverage)
        add_seq(reverse, coverage)

    def add_edge(self, dna, coverage=1):
        straight, reverse = straight_and_reverse(dna)

        def add_seq(seq, seq_cov):
            node_from = self.nodes[seq[:len(seq) - 1]]
            node_to = self.nodes[seq[1:]]
            # Добавляем ребро точкам, или просто увеличиваем coverage
            # (Если ребро уже выходит из одной точки, то обязательно входит в другую.
            if seq in node_from.edges_out:
                node_from.edges_out[seq].coverage += seq_cov
                node_to.edges_in[seq].coverage += seq_cov
            else:
                edge = Edge(seq, seq_cov, node_from, node_to)
                node_to.amount_in += 1
                node_from.amount_out += 1
                node_to.edges_in.update({seq: edge})
                node_from.edges_out.update({seq: edge})

        add_seq(straight, coverage)
        add_seq(reverse, coverage)

    def add_read(self, read, k=1):
        for index in range(len(read) - k + 1):
            self.add_node(read[index:index + k])
        for index in range(len(read) - k):
            self.add_edge(read[index:index + k + 1])

    def compression(self, kmer_len):
        for node in self.nodes:
            if not self.nodes[node].deleted and self.nodes[node].amount_in == 1 and self.nodes[node].amount_out == 1:
                # Находим те входящие и исходящие ребра, которые действующие (их по одной штуке)
                for k in self.nodes[node].edges_in:
                    if not self.nodes[node].edges_in[k].deleted:
                        for l in self.nodes[node].edges_out:
                            if not self.nodes[node].edges_out[l].deleted:
                                self.nodes[node].deleted = True
                                self.nodes[node].edges_in[k].deleted = True
                                self.nodes[node].edges_out[l].deleted = True
                                weight_in = self.nodes[node].edges_in[k].weight
                                weight_out = self.nodes[node].edges_in[k].weight
                                cov_in = self.nodes[node].edges_in[k].coverage
                                cov_out = self.nodes[node].edges_in[k].coverage
                                seq = k + l[kmer_len:]
                                edge = Edge(seq,
                                            (weight_out * cov_out + weight_in * cov_in) / (weight_out + weight_in),
                                            self.nodes[node].edges_in[k].node_from,
                                            self.nodes[node].edges_out[l].node_to,
                                            weight_in + weight_out)

                                self.nodes[node].edges_in[k].node_from.edges_out.update({seq: edge})
                                self.nodes[node].edges_out[l].node_to.edges_in.update({seq: edge})

    def graph_to_files(self, path, k):
        dot_path = path + '.dot'
        fasta_path = path + '.fasta'
        aver_len = 0
        aver_cov = 0
        amount_of_reads = 0
        with open(dot_path, 'w') as dot_out:
            with open(fasta_path, 'w') as fasta_out:
                dot_out.write('digraph eazy_g {\n')
                visited = {}
                for seq in self.nodes:
                    visited.update({seq: False})
                for node in self.nodes:
                    if not self.nodes[node].deleted and not visited[node]:
                        visited[node] = True
                        for edge in self.nodes[node].edges_out:
                            if not self.nodes[node].edges_out[edge].deleted:
                                dna_from = self.nodes[node].dna
                                dna_to = self.nodes[node].edges_out[edge].node_to.dna
                                # dot_out.write(f'{self.nodes[dna_from].dna}->{self.nodes[dna_to].dna} '
                                #           f'[label={self.nodes[node].edges_out[edge].coverage}]\n')
                                index_from = self.indexes[dna_from]
                                index_to = self.indexes[dna_to]
                                length = len(self.nodes[node].edges_out[edge].dna) - k
                                coverage = self.nodes[node].edges_out[edge].coverage
                                dna = self.nodes[node].edges_out[edge].dna
                                dot_out.write(f'{index_from}->{index_to} '
                                              f'[label=\"len = {length},'
                                              f'cov = {coverage}\"]\n')
                                fasta_out.write(f'>{index_from}-{index_to}___length={length}___cov={coverage}\n'
                                                f'{dna}\n')
                                aver_cov += coverage
                                aver_len += length+k
                                amount_of_reads += 1
                dot_out.write('}')
        return aver_len / amount_of_reads, aver_cov / amount_of_reads

    def pruning(self, length, coverage):
        visited = {}
        for seq in self.nodes:
            visited.update({seq: False})
        for node in self.nodes:
            if not self.nodes[node].deleted and not visited[node]:
                visited[node] = True
                if self.nodes[node].amount_in == 1 and self.nodes[node].amount_out == 0:
                    for edge in self.nodes[node].edges_in:
                        if len(self.nodes[node].edges_in[edge].dna) < length and \
                                self.nodes[node].edges_in[edge].coverage < coverage:
                            self.nodes[node].edges_in[edge].node_from.amount_out -= 1
                            self.nodes[node].edges_in[edge].deleted = True
                            self.nodes[node].amount_in = 0
                if self.nodes[node].amount_out == 1 and self.nodes[node].amount_in == 0:
                    for edge in self.nodes[node].edges_out:
                        if len(self.nodes[node].edges_out[edge].dna) < length and \
                                self.nodes[node].edges_out[edge].coverage < coverage:
                            self.nodes[node].edges_out[edge].node_to.amount_in -= 1
                            self.nodes[node].edges_out[edge].deleted = True
                            self.nodes[node].amount_out = 0


def iterate_reads(path, graph, k):
    records = SeqIO.parse(path, 'fastq')
    for read in records:
        graph.add_read(str(read.seq), k)
    graph.compression(k)


def in_float(s='Введите число', integer=False, check=[False, 0, 0]):
    flag = True
    while flag:
        flag = False
        try:
            if integer:
                val = int(input(s + ': '))
            else:
                val = float(input(s + ': '))
            if check[0] and (val < check[1] or val > check[2]):
                raise ValueError
        except ValueError:
            flag = True
            if check[0]:
                print(f'Попробуйте снова! Введенное число должно принадлежать интервалу [{check[1]}; {check[2]}]\n')
            else:
                print(f'Попробуйте снова!\n')
    return val


def assemble_and_save(name_of_file, k=55):
    start = time()
    print(f'File {name_of_file}')
    print('Assembling...')
    fastqc_file = f'data/{name_of_file}.fastq'
    graph = Graph()
    iterate_reads(fastqc_file, graph, k)
    print(f'Assemble done in {time() - start} seconds')
    print('Saving and drawing...')
    start2 = time()
    name_save = f'graphs/{name_of_file}'
    length, coverage = graph.graph_to_files(name_save, k)
    visualize(name_save + '.dot')
    print(f'Saving and vizializing done in {time() - start2} seconds')

    print(f'\nAverage length = {length}\nAverage coverage = {coverage}')
    print('Pruning!')
    treshold_len = in_float('Input treshold for length')
    treshold_cov = in_float('Input treshold for coverage')

    start2 = time()
    graph.pruning(treshold_len+k, treshold_cov)
    graph.compression(k)
    name_save = f'graphs/{name_of_file}+_pr'
    graph.graph_to_files(name_save, k)
    visualize(name_save + '.dot')
    print(f'Pruning, saving and vizializing done in {time() - start2} seconds')

    print(f'Overall time {time() - start} seconds')
    print()


if __name__ == '__main__':
    name_of_file = 's_6.first10000'
    k = 55
    assemble_and_save(name_of_file, k)
