# from ..utils.utils import reader
from tqdm import tqdm
from src.utils import utils
from bx.intervals.intersection import Interval, Intersecter
from rich.progress import track


def read_chain_file(chain_file, print_table=False):
    '''
    Read chain file.

    Parameters
    ----------
    chain_file : file
            Chain format file. Input chain_file could be either plain text, compressed file
            (".gz",".Z", ".z", ".bz", ".bz2", ".bzip2"), or a URL pointing to the chain file
            ("http://","https://", "ftp://"). If url was used, chain file must be plain text.

    print_table : bool, optional
            Print mappings in human readable table.

    Returns
    -------
    maps : dict
            Dictionary with source chrom name as key, IntervalTree object as value. An
            IntervalTree contains many intervals. An interval is a start and end position
            and a value. eg. Interval(11, 12, strand="-", value = "abc")

    target_chromSize : dict
            Chromosome sizes of target genome

    source_chromSize : dict
            Chromosome sizes of source genome
    '''

    maps = {}
    target_chromSize = {}
    source_chromSize = {}
    if print_table:
        blocks = []
    total_size = len(list(utils.reader(chain_file)))
    for line in track(utils.reader(chain_file), description="[bold green]Reading", total=total_size,):
        # Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        if not line.strip():
            continue
        line = line.strip()
        if line.startswith(('#', ' ')):
            continue
        fields = line.split()

        if fields[0] == 'chain' and len(fields) in [12, 13]:
            # score = int(fields[1])		  # Alignment score
            source_name = fields[2]		  # E.g. chrY
            source_size = int(fields[3])  # Full length of the chromosome
            source_strand = fields[4]	  # Must be +
            if source_strand != '+':
                raise Exception(
                    "Source strand in a chain file must be +. (%s)" % line)
            source_start = int(fields[5])  # Start of source region
            # source_end = int(fields[6])	  # End of source region

            target_name = fields[7]		  # E.g. chr5
            target_size = int(fields[8])  # Full length of the chromosome
            target_strand = fields[9]	  # + or -
            target_start = int(fields[10])
            # target_end = int(fields[11])
            target_chromSize[target_name] = target_size
            source_chromSize[source_name] = source_size

            if target_strand not in ['+', '-']:
                raise Exception("Target strand must be - or +. (%s)" % line)
            # chain_id = None if len(fields) == 12 else fields[12]
            if source_name not in maps:
                maps[source_name] = Intersecter()

            sfrom, tfrom = source_start, target_start

        # Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to)
        elif fields[0] != 'chain' and len(fields) == 3:
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
            if print_table:
                if target_strand == '+':
                    blocks.append((source_name, sfrom, sfrom + size, source_strand,
                                  target_name, tfrom, tfrom + size, target_strand))
                elif target_strand == '-':
                    blocks.append((source_name, sfrom, sfrom + size, source_strand, target_name,
                                  target_size - (tfrom + size), target_size - tfrom, target_strand))

            if target_strand == '+':
                maps[source_name].add_interval(
                    Interval(sfrom, sfrom + size, (target_name, tfrom, tfrom + size, target_strand)))
            elif target_strand == '-':
                maps[source_name].add_interval(Interval(
                    sfrom, sfrom + size, (target_name, target_size - (tfrom + size), target_size - tfrom, target_strand)))

            sfrom += size + sgap
            tfrom += size + tgap

        elif fields[0] != 'chain' and len(fields) == 1:
            size = int(fields[0])
            if print_table:
                if target_strand == '+':
                    blocks.append((source_name, sfrom, sfrom + size, source_strand,
                                  target_name, tfrom, tfrom + size, target_strand))
                elif target_strand == '-':
                    blocks.append((source_name, sfrom, sfrom + size, source_strand, target_name,
                                  target_size - (tfrom + size), target_size - tfrom, target_strand))

            if target_strand == '+':
                maps[source_name].add_interval(
                    Interval(sfrom, sfrom + size, (target_name, tfrom, tfrom + size, target_strand)))
            elif target_strand == '-':
                maps[source_name].add_interval(Interval(
                    sfrom, sfrom + size, (target_name, target_size - (tfrom + size), target_size - tfrom, target_strand)))
        else:
            raise Exception("Invalid chain format. (%s)" % line)
    # if (sfrom + size) != source_end  or (tfrom + size) != target_end:
    #	 raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)

    if print_table:
        for i in blocks:
            print('\t'.join([str(n) for n in i]))

    return (maps, target_chromSize, source_chromSize)


def intersectBed(lst1, lst2):
    '''
    Return intersection of two bed regions.

    Parameters
    ----------
    lst1 : list
            The 1st genomic region. List of chrom, start, end.
            Example: ['chr1',10, 100]

    lst2 : list
             The 2nd genomic region. List of chrom, start, end.
             Example: ['chr1',50, 120]

    Examples
    --------
    >>> intersectBed(['chr1',10, 100],['chr1',50, 120])
    ('chr1', 50, 100)
    >>> intersectBed(['chr1',10, 100],['chr1',20, 30])
    ('chr1', 20, 30)

    '''
    (chr1, st1, end1) = lst1
    (chr2, st2, end2) = lst2
    if int(st1) > int(end1) or int(st2) > int(end2):
        raise Exception("Start cannot be larger than end")
    if chr1 != chr2:
        return None
    if int(st1) > int(end2) or int(end1) < int(st2):
        return None
    return (chr1, max(st1, st2), min(end1, end2))


def map_coordinates(mapping, q_chr, q_start, q_end, q_strand='+'):
    '''
    Map coordinates from source (i.e. original) assembly to target (i.e. new) assembly.

    Parameters
    ----------
    mapping : dict
            Dictionary with source chrom name as key, IntervalTree object as value.

    q_chr : str
            Chromosome ID of query interval

    q_start : int
            Start position of query interval.

    q_end : int
            End position of query interval.

    q_strand : str
            Strand of query interval.

    '''

    matches = []
    complement = {'+': '-', '-': '+'}

    if q_chr in mapping:
        targets = mapping[q_chr].find(q_start, q_end)
    elif q_chr.replace('chr', '') in mapping:
        targets = mapping[q_chr.replace('chr', '')].find(q_start, q_end)
    elif ('chr' + q_chr) in mapping:
        targets = mapping['chr' + q_chr].find(q_start, q_end)
    else:
        return None
    if len(targets) == 0:
        return None
    elif len(targets) == 1:
        s_start = targets[0].start
        s_end = targets[0].end
        t_chrom = targets[0].value[0]
        # t_chrom = update_chromID(q_chr, t_chrom, chr_style = chrom_style)
        t_start = targets[0].value[1]
        t_end = targets[0].value[2]
        t_strand = targets[0].value[3]

        (chr, real_start, real_end) = intersectBed(
            (q_chr, q_start, q_end), (q_chr, s_start, s_end))
        l_offset = abs(real_start - s_start)
        # r_offset = real_end - s_end
        size = abs(real_end - real_start)

        matches.append((chr, real_start, real_end, q_strand))
        if t_strand == '+':
            i_start = t_start + l_offset
            if q_strand == '+':
                matches.append((t_chrom, i_start, i_start + size, t_strand))
            else:
                matches.append((t_chrom, i_start, i_start +
                               size, complement[t_strand]))
        elif t_strand == '-':
            i_start = t_end - l_offset - size
            if q_strand == '+':
                matches.append((t_chrom, i_start, i_start + size, t_strand))
            else:
                matches.append((t_chrom, i_start, i_start +
                               size, complement[t_strand]))
        else:
            raise Exception(
                "Unknown strand: %s. Can only be '+' or '-'." % q_strand)

    elif len(targets) > 1:
        for t in targets:
            s_start = t.start
            s_end = t.end
            t_chrom = t.value[0]
            # t_chrom = update_chromID(q_chr, t_chrom, chr_style=chrom_style)
            t_start = t.value[1]
            t_end = t.value[2]
            t_strand = t.value[3]

            (chr, real_start, real_end) = intersectBed(
                (q_chr, q_start, q_end), (q_chr, s_start, s_end))

            l_offset = abs(real_start - s_start)
            # r_offset = abs(real_end - s_end)
            size = abs(real_end - real_start)
            matches.append((chr, real_start, real_end, q_strand))
            if t_strand == '+':
                i_start = t_start + l_offset
                if q_strand == '+':
                    matches.append(
                        (t_chrom, i_start, i_start + size, t_strand))
                else:
                    matches.append(
                        (t_chrom, i_start, i_start + size, complement[t_strand]))
            elif t_strand == '-':
                i_start = t_end - l_offset - size
                if q_strand == '+':
                    matches.append(
                        (t_chrom, i_start, i_start + size, t_strand))
                else:
                    matches.append(
                        (t_chrom, i_start, i_start + size, complement[t_strand]))
            else:
                raise Exception(
                    "Unknown strand: %s. Can only be '+' or '-'." % q_strand)

    # if print_match:
    #     print(matches)
        # input: 'chr1',246974830,247024835
        # output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
        # [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]

    return matches
