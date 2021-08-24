from bz2 import open as bz2_open
from collections import defaultdict
from gzip import open as gzip_open
from os import mkdir
from os.path import basename, exists
import re
import csv

from Bio import SeqIO
from Bio.Seq import reverse_complement
from dict_trie import Trie
from fastools import guess_file_format, guess_header_format
from jit_open import Handle, Queue
from .match import multi_align


_get_barcode = {
    'normal': lambda record: record.id.split('#')[1].split('/')[0],
    'x': lambda record: record.description.split(':')[-1],
    'umi': lambda record: record.description.split(' ')[0].split(':')[-1],
    'unknown': lambda record: str(record.seq)}

_type_handler = defaultdict(lambda: open, {
    'bz2': bz2_open,
    'bzip2': bz2_open,
    'gz': gzip_open,
    'gzip': gzip_open})


def _name(handle):
    if hasattr(handle.buffer, '_fp'):
        return handle.buffer._fp.name
    return handle.name


class Extractor(object):
    def __init__(self, handle, in_read=False, fmt=None, start=None, end=None):
        """Configure a barcode extractor.

        :arg stream handle: Handle to an NGS data file.
        :arg bool in_read: Inspect the read instead of the header.
        :arg str fmt: Header format.
        :arg int start: Start of the barcode.
        :arg int end: End of the barcode.
        """
        self._start = start
        self._end = end

        if self._start:
            self._start -= 1

        if not fmt:
            if not in_read:
                self._get_barcode = _get_barcode[guess_header_format(handle)]
            else:
                self._get_barcode = _get_barcode['unknown']
        else:
            self._get_barcode = _get_barcode[fmt]

    def get(self, record):
        return self._get_barcode(record)[self._start:self._end]


def count(handle, extractor, sample_size, threshold, use_freq=False):
    """Get the most frequent barcodes from an NGS data file.

    :arg stream handle: Handle to an NGS data file.
    :arg Extractor extractor: A barcode extractor.
    :arg int sample_size: Number of records to probe.
    :arg int threshold: Threshold for the selection method.
    :arg bool use_freq: Select frequent barcodes instead of a fixed amount.

    :returns list: A list of barcodes.
    """
    barcodes = defaultdict(int)

    for i, record in enumerate(SeqIO.parse(handle, guess_file_format(handle))):
        if i > sample_size:
            break
        barcodes[extractor.get(record)] += 1

    if use_freq:
        return filter(lambda x: barcodes[x] >= threshold, barcodes)
    return sorted(barcodes, key=barcodes.get, reverse=True)[:threshold]


def _open_files(path, filenames, barcode, queue):
    """For a list of input files, open the corresponding output files.

    :arg str path: Output directory.
    :arg list filename: List of input filenames.
    :arg str barcode: Name of the barcode.
    :arg Queue queue: Queue for open files.

    :returns list: List of handles of output files.
    """
    if not exists(path):
        mkdir(path)
    handles = []
    for filename in filenames:
        base, ext = basename(filename).split('.', True)
        handles.append(
            Handle('{}/{}_{}.{}'.format(path, base, barcode, ext), queue,
                f_open=_type_handler[ext.split('.')[-1]]))

    return handles

def _open_diag_files(path, filenames, barcode, queue, case):
    """Open the corresponding output files for a diagnostic using input file and diag case.

    :arg str path: Output directory.
    :arg list filename: List of input filenames.
    :arg str barcode: Name of the barcode.
    :arg Queue queue: Queue for open files.
    :arg str case: case diagnostic.

    :returns list: List of handles of output files.
    """
    if not exists(path):
        mkdir(path)
    handles = []
    base, ext = basename(filenames[0]).split('.', True)
    base = 'result'
    barcode = barcode + '_' + case
    barcode1 = barcode + '_R1'
    barcode2 = barcode + '_R2'
    handles.append(
        Handle('{}/{}_{}.{}'.format(path, base, barcode1, ext), queue,
            f_open=_type_handler[ext.split('.')[-1]]))
    handles.append(
        Handle('{}/{}_{}.{}'.format(path, base, barcode2, ext), queue,
            f_open=_type_handler[ext.split('.')[-1]]))

    return handles

def _write(handles, records, file_format):
    for i, record in enumerate(records):
        SeqIO.write(record, handles[i], file_format)


def demultiplex(
        input_handles, barcodes_handle, extractor, mismatch, use_edit,
        path='.'):
    """Demultiplex a list of NGS data files.

    :arg list input_handles: List of handles to NGS data files.
    :arg stream barcodes_handle: Handle to a file containing barcodes.
    :arg Extractor extractor: A barcode extractor.
    :arg int mismatch: Number of allowed mismatches.
    :arg bool use_edit: Use Levenshtein distance instead of Hamming distance.
    :arg str path: Output directory.
    """
    filenames = list(map(lambda x: _name(x), input_handles))
    queue = Queue()
    default_handles = _open_files(path, filenames, 'UNKNOWN', queue)

    barcodes = {}
    for line in barcodes_handle.readlines():
        try:
            name, barcode = line.strip().split()
        except ValueError:
            raise ValueError('invalid barcodes file format')
        barcodes[barcode] = _open_files(path, filenames, name, queue)

    trie = Trie(barcodes.keys())
    distance_function = trie.best_hamming
    if use_edit:
        distance_function = trie.best_levenshtein

    file_format = guess_file_format(input_handles[0])

    readers = list(map(
        lambda x: SeqIO.parse(x, file_format), input_handles))

    while True:
        records = list(map(lambda x: next(x), readers))
        if not records:
            break

        barcode = distance_function(extractor.get(records[0]), mismatch)
        if barcode:
            _write(barcodes[barcode], records, file_format)
        else:
            _write(default_handles, records, file_format)

    queue.flush()


def compare_paired_assignation(records, mismatch, indel_score, path, filenames, queue, list_id_barcode, file_format, dict_id_count_R1, dict_id_count_R2, list_overall_stats):
    """
    """
    # sametag R1 and R2 have same tag
    # onetage only one of the two reads have a tag
    # incomptag R1 and R2 have different tag
    # notag none R1 or R2 have tag
    #case=["sametag", "onetag", "incomptag", "notag" ]
 
    default_handles = _open_files(path, filenames, 'UNKNOWN', queue)
    stat_handles = _open_files(path, filenames, 'STAT', queue)

    current_r1_id = records[0].id
    reference_r1 = str(records[0].seq.upper())
    current_r2_id = records[1].id
    reference_r2 = str(records[1].seq.upper())
    foundr1 = False
    foundr2 = False
    save_tag_used = ""
    save_tagr2_used = ""
    last_tag_used = ""
    for id, barcode in list_id_barcode:
        
        # R1
        if multi_align(reference_r1, barcode, mismatch, indel_score):
            save_tag_used = id
            last_tag_used = id
            foundr1 = True
            if save_tag_used in dict_id_count_R1:
                if dict_id_count_R1[save_tag_used][0] == 0:
                    dict_id_count_R1[save_tag_used] = [1, 0, 0, 0]
                else:
                    dict_id_count_R1[save_tag_used][0] = dict_id_count_R1[save_tag_used][0]+1
        # R2
        if multi_align(reference_r2, barcode, mismatch, indel_score):
            save_tagr2_used = id
            last_tag_used = id
            foundr2 = True
            if save_tagr2_used in dict_id_count_R2:
                if dict_id_count_R2[save_tagr2_used][0] == 0:
                    dict_id_count_R2[save_tagr2_used] = [1, 0, 0, 0]
                else:
                    dict_id_count_R2[save_tagr2_used][0] = dict_id_count_R2[save_tagr2_used][0]+1
 
    if not foundr1 and not foundr2:
        # notag none R1 or R2 have tag
        _write(_open_diag_files(path, filenames, "UNKNOWN", queue, "notag"), records, file_format)
        list_overall_stats[0] +=1

    elif (not foundr1 and foundr2) or (foundr1 and not foundr2):
        # onetag only one of the two reads have a tag
        if foundr1:
            dict_id_count_R1[last_tag_used][1] =  dict_id_count_R1[last_tag_used][1] + 1 
        if foundr2:
            dict_id_count_R2[last_tag_used][1] =  dict_id_count_R2[last_tag_used][1] + 1
        _write(_open_diag_files(path, filenames, str(last_tag_used), queue, "onetag"), records, file_format)
        list_overall_stats[1] +=1

    elif foundr1 and foundr2:
        if save_tag_used == save_tagr2_used:         
            # sametag R1 and R2 have same tag
            _write(_open_diag_files(path, filenames, str(save_tag_used) , queue, "sametag"), records, file_format)
            dict_id_count_R1[save_tag_used][2] =  dict_id_count_R1[save_tag_used][2] + 1 
            dict_id_count_R2[save_tag_used][2] =  dict_id_count_R2[save_tag_used][2] + 1
            list_overall_stats[2] +=1

        else: 
            # incomptag R1 and R2 have different tag
            _write(_open_diag_files(path, filenames, "conflict", queue, "incomptag"), records, file_format)
            dict_id_count_R1[save_tag_used][3] =  dict_id_count_R1[save_tag_used][3] + 1 
            dict_id_count_R2[save_tagr2_used][3] =  dict_id_count_R2[save_tagr2_used][3] + 1
            list_overall_stats[3] +=1
            
    return dict_id_count_R1, dict_id_count_R2, list_overall_stats

def match(input_handles, barcodes_handle, mismatch, use_edit, path='.'):
    """Demultiplex a list of NGS data files.

    :arg list input_handles: List of handles to NGS data files.
    :arg stream barcodes_handle: Handle to a file containing barcodes.
    :arg int mismatch: Number of allowed mismatches.
    :arg bool use_edit: Use Levenshtein distance instead of Hamming distance.
    :arg str path: Output directory.
    """
    filenames = list(map(lambda x: _name(x), input_handles))
    queue = Queue()
    default_handles = _open_files(path, filenames, 'UNKNOWN', queue)

    indel_score = 1
    if not use_edit:
        indel_score = 1000
    barcodes = []
    list_id_barcode = []
    dict_id={}
    for line in map(lambda x: x.strip().upper().split(), barcodes_handle.readlines()):
        try:
            name = line.pop(0)
        except ValueError:
            raise ValueError('invalid barcodes file format')
        barcodes.append((_open_files(path, filenames, name, queue), line))
        list_id_barcode.append([name, line])
        dict_id.update({name: [0, 0, 0, 0]})

    dict_id_count_R1 = dict_id.copy()
    dict_id_count_R2 = dict_id.copy()
    
    file_format = guess_file_format(input_handles[0])
    readers = list(map(
        lambda x: SeqIO.parse(x, file_format), input_handles))
    # print(input_handles)
    # print(input_handles)
    # if len(input_handles) == 2:
    #     readers = list(map(
    #     lambda x: SeqIO.parse(x, file_format), input_handles[0]))
    # print(readers)
    overall_count_notag = 0
    overall_count_onetag = 0
    overall_count_sametag = 0
    overall_count_incomptag = 0
    list_overall_stats = [overall_count_notag, overall_count_onetag, overall_count_sametag, overall_count_incomptag]

    while True:
        records = list(map(lambda x: next(x), readers))
        if not records:
            break
        if len(records)==2:

            dict_id_count_R1, dict_id_count_R2, list_overall_stats= compare_paired_assignation(records, mismatch, indel_score, path, filenames, queue, list_id_barcode, file_format, dict_id_count_R1, dict_id_count_R2, list_overall_stats)
        else:
            reference = str(records[0].seq.upper())
            # TODO make reverse search optionnal 
            # reference_rc = reverse_complement(reference)
            found = False
            for handles, barcode in barcodes:
                save_tag_used = ""
                if multi_align(reference, barcode, mismatch, indel_score):
                    _write(handles, records, file_format)
                    if save_tag_used == "":
                        save_tag_used = barcode[0]
                    else:
                        print("WARNING, several TAG in same read")
                        print(records) 
                    found = True
                    continue
                #TODO chekc primer in reverse later
                # TODO make reverse search optionnal 
                # elif multi_align(reference_rc, barcode, mismatch, indel_score):
                #     _write(handles, records, file_format_rc)
                #     found = True
                #     continue

            if not found:
                _write(default_handles, records, file_format)

    # print(dict_id_count_R1, dict_id_count_R2, list_overall_stats)
    result_path1 = path + '/result_stat_R1.csv'
    result_path2 = path + '/result_stat_R2.csv'
    with open(result_path1, 'w') as f:
        f.write(';nb reads with no tag in the pair;nb reads with one tag on the two reads;nb reads with same tag for both read;nb reads with different tag for both read\n')
        f.write("%s;%s;%s;%s;%s\n"%('',list_overall_stats[0], list_overall_stats[1], list_overall_stats[2], list_overall_stats[3]))
        f.write('\n')
        f.write('ID;nb reads with the current tag;nb reads with one tag on the two reads;nb reads with same tag for both read;nb reads with different tag for both read\n')
        for key in dict_id_count_R1.keys():
            tag = key
            notag = dict_id_count_R1[key][0]
            onetag = dict_id_count_R1[key][1]
            sametag = dict_id_count_R1[key][2]
            incomptag = dict_id_count_R1[key][3]
            f.write("%s;%s;%s;%s;%s\n"%(tag,notag,onetag,sametag,incomptag))

    with open(result_path2, 'w') as f:
        f.write(';nb reads with no tag in the pair;nb reads with one tag on the two reads;nb reads with same tag for both read;nb reads with different tag for both read\n')
        f.write("%s;%s;%s;%s;%s\n"%('',list_overall_stats[0], list_overall_stats[1], list_overall_stats[2], list_overall_stats[3]))
        f.write('\n')
        f.write('ID;nb reads with the current tag;nb reads with one tag on the two reads;nb reads with same tag for both read;nb reads with different tag for both read\n')
        for key in dict_id_count_R2.keys():
            tag = key
            notag = dict_id_count_R2[key][0]
            onetag = dict_id_count_R2[key][1]
            sametag = dict_id_count_R2[key][2]
            incomptag = dict_id_count_R2[key][3]
            f.write("%s;%s;%s;%s;%s\n"%(tag,notag,onetag,sametag,incomptag))            

    queue.flush()

