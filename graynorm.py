#!/usr/bin/env python
# ------------------------------------------------------------------------
# GrayNorm is an algorithm that helps the researcher to identify
# suitable reference genes that are minimally influenced by the experiment.
# Gene expression levels are determined by quantitative reverse
# transcription PCR (RT-qPCR) and one is interested in up- or down
# regulation of various genes as a consequence of varying experimental
# conditions.
# Copyright (C) 2013 Tony Remans <tony.remans@uhasselt.be>,
#                    Geert Jan Bex <geertjan.bex@uhasselt.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------
#
from argparse import ArgumentParser, FileType
from math import sqrt
import csv, itertools, logging, re, sys

class Data(object):

    def __init__(self, headers):
        self._headers = headers
        self._header_idx = {y: x for x, y in enumerate(self._headers)}
        self._non_data = []
        self._data = []
        self._control_idx = []
        self._cond_idx = []
        self._cond_groups = {}
        self._cond_group_idx = []
        self._gene_idx = None
    
    def add(self, data):
        '''add a sample, represented as a list'''
        for idx, header in enumerate(self._headers):
            if not header in self._non_data:
                data[idx] = float(data[idx])
        cond_values = []
        for x in self._cond_idx:
            try:
                cond_values.append(float(data[x]))
            except ValueError:
                cond_values.append(data[x])
        if self._control_cond_vals == cond_values:
            self._control_idx.append(len(self._data))
        cond_label = str(cond_values)
        if cond_label not in self._cond_groups:
            self._cond_groups[cond_label] = []
            self._cond_group_idx.append(cond_values)
        self._cond_groups[cond_label].append(len(self._data))
        self._data.append(data)

    @property
    def control_values(self):
        return self._control_cond_vals

    def set_condition_col_names(self, conditions):
        '''set the names of the columns that hold the experimental
           conditions'''
        self._conditions = conditions
        self._cond_idx = []
        for idx, header in enumerate(self._headers):
            if header in self._conditions:
                self._cond_idx.append(idx)
        self._non_data.extend(conditions)

    def set_gene_col_names(self, genes):
        self._gene_names = genes

    def genes(self, idx=None):
        if self._gene_names:
            return self._gene_names
        if idx is None and self._gene_idx is None:
            sys.stderr.write('### error: no gene sequence specified\n')
            sys.exit(3)
        if idx:
            self._gene_idx = idx
        return [self._headers[i] for i in self._gene_idx]

    @property
    def condition_columns(self):
        '''return the indices of the columns that hold the experimental
           conditions'''
        return self._cond_idx

    @property
    def condition_values(self):
        '''return the distinct condition values as a list of lists'''
        return self._cond_group_idx

    @property
    def nr_conditions(self):
        '''return the number of distinct conditions in the experiment'''
        return len(self._cond_group_idx)

    def condition_group(self, cond_values):
        cond_label = str(cond_values)
        if cond_label in self._cond_groups:
            return self._cond_groups[cond_label]
        else:
            sys.stderr.write('### error: no such condition {0}\n'.format(cond_label))
            sys.exit(1)

    def set_control_vals(self, values):
        self._control_cond_vals = values

    @property
    def control_rows(self):
        return self._control_idx

    def set_sample_col_name(self, sample):
        self._sample = sample
        for idx, header in enumerate(self._headers):
            if header == self._sample:
                self._sample_idx = idx
        self._non_data.append(sample)

    @property
    def sample_column(self):
        return self._sample_idx

    def compute_nf(self, genes):
        col_sel = [header in genes for header in self._headers]
        nfs = []
        for row in self._data:
            nf = reduce(lambda x, y: x*y,
                        itertools.compress(row, col_sel),
                        1.0)**(1.0/len(genes))
            nfs.append(nf)
        return nfs

    def compute_inv_nf(self, genes):
        return [1.0/x for x in self.compute_nf(genes)]

    def compute_inv_nf_vs_control(self, genes):
        inv_nfs = self.compute_inv_nf(genes)
        avg_control = sum([inv_nfs[idx] for idx in self.control_rows])
        avg_control /= len(self.control_rows)
        return [x/avg_control for x in inv_nfs]

    def compute_condition_stats(self, genes):
        inv_nfs_vs_ctrl = self.compute_inv_nf_vs_control(genes)
        stats = []
        for cond_values in self.condition_values:
            avg, stddev, stderr = compute_stats([inv_nfs_vs_ctrl[i]
                                                     for i in self.condition_group(cond_values)])
            stats.append({
                'cond':   cond_values,
                'avg':    avg,
                'stddev': stddev,
                'stderr': stderr,
                'cv_intra': stddev/avg
                })
        return stats

    def compute_overall_stats(self, genes):
        stats = self.compute_condition_stats(genes)
        avgs = [stat['avg'] for stat in stats]
        avg, stddev, _ = compute_stats(avgs)
        cumm = sum([abs(1.0 - x) for x in avgs])
        return {
            'avg':    avg,
            'stddev': stddev,
            'cumm':   cumm,
            'cv_inter': stddev/avg
        }

    def compute_all(self, cand_genes):
        gene_combinations = []
        for n in xrange(1, len(cand_genes) + 1):
            for genes in itertools.combinations(cand_genes, n):
                cond_stats = self.compute_condition_stats(genes)
                overall_stats = self.compute_overall_stats(genes)
                gene_combinations.append({
                    'genes':   genes,
                    'conds':   cond_stats,
                    'overall': overall_stats
                })
        gene_combinations.sort(key=lambda x: x['overall']['cv_inter'])
        return gene_combinations
        
    def __str__(self):
        s = ''
        s += '# genes: ' + ','.join(self.genes()) + '\n'
        s += '# controls: ' + ','.join([str(idx) for idx in self.condition_columns]) + '\n'
        s += '# control rows: ' + ', '.join([str(idx) for idx in self.control_rows]) + '\n'
        s += '\t'.join(self._headers)
        for data in self._data:
            s += '\n' + '\t'.join(['{0}'.format(x) for x in data])
        return s

    @property
    def header_row(self):
        row = ['gene combination', 'CV inter']
        row.extend(['CV intra cond {0}'.format(i)
                       for i in xrange(1, self.nr_conditions + 1)])
        for i in xrange(1, self.nr_conditions + 1):
            row.extend(['{quant} cond {i}'.format(quant=q, i=i)
                            for q in ['avg', 'stddev', 'stderr']])
        row.extend(['avg 1/NF', 'stddev 1/NF', 'cummulative 1/NF'])
        return row

    def output_row(self, stats):
        row = [' + '.join(stats['genes']), stats['overall']['cv_inter']]
        row.extend([cond['cv_intra'] for cond in stats['conds']])
        for cond in stats['conds']:
            row.extend([cond['avg'], cond['stddev'], cond['stderr']])
        row.extend([stats['overall']['avg'], stats['overall']['stddev'],
                    stats['overall']['cumm']])
        return row

    def write_results(self, stats, output_file):
        with output_file:
            csv_writer = csv.writer(output_file)
            csv_writer.writerow(self.header_row)
            for stat in stats:
                csv_writer.writerow(self.output_row(stat))


def compute_stats(numbers):
    n = len(numbers)
    s = sum(numbers)
    s2 = sum([x**2 for x in numbers])
    stddev = sqrt((s2 - s**2/n)/(n - 1.0))
    return (s/n, stddev, stddev/sqrt(n))

def read_file(data_file_name):
    data = None
    meta_info_re = re.compile(r'\s*#\s*(\w+)\s*:\s*(.+?)\s*$')
    sample_col_name = None
    gene_col_names = None
    cond_col_names = None
    control_values = None
    with open(data_file_name, 'rb') as data_file:
        dialect = csv.Sniffer().sniff(data_file.read(2048))
        data_file.seek(0)
        data_reader = csv.reader(data_file, dialect=dialect)
        for row in data_reader:
            if row[0].strip().startswith('#'):
                match = meta_info_re.match(row[0])
                if match:
                    key = match.group(1)
                    value = match.group(2)
                    if key == 'sampleid':
                        sample_col_name = value
                    elif key == 'refgenes':
                        gene_col_names = re.split(r'\s*,\s*', value)
                    elif key == 'controls':
                        cond_col_names = []
                        control_values = []
                        controls = re.split(r'\s*,\s*', value)
                        for control in controls:
                            parts = re.split(r'\s*=\s*', control)
                            if len(parts) < 2:
                                sys.stderr.write('### error: no control ' \
                                                 'value for control ' \
                                                 'parameter \'{0}\', ' \
                                                 'check input ' \
                                                 'data format\n'.format(control))
                                sys.exit(1)
                            cond_col_names.append(parts[0])
                            try:
                                control_values.append(float(parts[1]))
                            except ValueError:
                                control_values.append(parts[1])
            elif row[0].isspace():
                pass
            else:
                headers = row
                status = check_input(sample_col_name, cond_col_names,
                                     control_values, gene_col_names,
                                     headers)
                if status:
                    sys.exit(status)
                data = Data(headers)
                data.set_sample_col_name(sample_col_name)
                data.set_condition_col_names(cond_col_names)
                data.set_control_vals(control_values)
                data.set_gene_col_names(gene_col_names)
                break
        for row in data_reader:
            if not row[0].strip().startswith('#') and not row[0].isspace():
                data.add(row)
    return data

def check_input(sample_col_name, cond_col_names, control_values,
                gene_col_names, headers):
    col_info = [
        ('sample ID', sample_col_name),
        ('control', cond_col_names),
        ('genes', gene_col_names),
    ]
    for info in col_info:
        status = check_column(info[0], info[1], headers)
        if status:
            return status
    return 0

def check_column(name, col_names, headers):
    if not col_names:
        sys.stderr.write('### error: no column name(s) for {0} specified, check input data format\n'.format(name))
        return 4
    unknowns = unknown_columns(col_names, headers)
    if unknowns:
        sys.stderr.write('### error: no column for {0}(s) {1} present, check input data format\n'.format(name, ','.join(["'{0}'".format(x) for x in unknowns])))
        return 5

def unknown_columns(cols, headers):
    if isinstance(cols, str):
        cols = [cols]
    c_set = set(cols)
    h_set = set(headers)
    if c_set < h_set:
        return set()
    else:
        return c_set - h_set

def compute_gene_idx(gene_str):
    parts = re.split(r'\s*,\s*', gene_str)
    idx = []
    for part in parts:
        match = re.match(r'^(\d+)$', part.strip())
        if match is not None:
            idx.append(int(match.group(1)) - 1)
            continue
        match = re.match(r'^(\d+)-(\d+)$', part.strip())
        if match is not None:
            idx.extend(range(int(match.group(1)) - 1, int(match.group(2))))
            continue
        sys.stderr.write('### error: invalid gene specs: {0}\n'.format(gene_str))
    return idx

def main():
    arg_parser = ArgumentParser(description='compute support')
    arg_parser.add_argument('-in', dest='data_file', required=True,
                            help='CSV file to read input from')
    arg_parser.add_argument('-refgenes', dest='cand_genes',
                            help='comma-separated list of candidate'
                                 ' normalization genes')
    arg_parser.add_argument('-out', dest='output', type=FileType('wb'),
                            required=True, help='CSV file to write the'
                                                ' results')
    arg_parser.add_argument('-verbose', dest='verbose', action='store_true',
                            help='print feedback during run')
    options = arg_parser.parse_args()
    if options.verbose:
        logging.basicConfig(level=logging.INFO)
    gene_idx = None
    if options.cand_genes:
        gene_idx = compute_gene_idx(options.cand_genes)
    try:
        data = read_file(options.data_file)
    except IOError as e:
        print dir(e)
        sys.stderr.write("### error: {0}: '{1}'\n".format(e.strerror,
                                                          e.filename))
        sys.exit(2)
    genes = data.genes(gene_idx)
    logging.info('sample column: ' + str(data.sample_column))
    logging.info('condition columns: ' + str(data.condition_columns))
    logging.info('control values: ' + str(data.control_values))
    logging.info('condition values: ' + str(data.condition_values))
    logging.info('control data rows: ' + str(data.control_rows))
    logging.info('genes: ' + str(genes))
    stats = data.compute_all(genes)
    data.write_results(stats, options.output)
    return 0

if __name__ == '__main__':
    status = main()
    sys.exit(status)

    
