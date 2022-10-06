# TALON: Techonology-Agnostic Long Read Analysis Pipeline
# Author: Dana Wyman
# -----------------------------------------------------------------------------
# This program reads in SAM-formatted long read alignments and adds a custom 
# tag to reflect the fraction of As in the sequence immediately following the
# alignment. This can help indicate the likelihood of an internal priming 
# artifact. 

import pyfaidx
import pysam
import multiprocessing as mp
from datetime import datetime, timedelta
import time
import os
import glob
from optparse import OptionParser
import sys

def get_options():
    """ Read input args """

    parser = OptionParser(description=("This program reads in SAM-formatted "
                          "long read alignments and adds a custom tag to "
                          "reflect the fraction of As in the sequence "
                          "immediately following the alignment. This can help "
                          "indicate the likelihood of an internal priming "
                          "artifact."))
    parser.add_option("--f", dest = "sam_file",
                      help = "SAM file of transcripts")
    parser.add_option("--g", dest = "genome_file",
                      help = "Reference genome fasta file")
    parser.add_option("--t", dest = "threads", type = int,
                      help = "Number of threads to run", default = 1)
    parser.add_option("--ar", dest = "fracA_range_size", type = int,
                      help = ("Size of post-transcript interval to compute "
                              "fraction As on. Default = 20"), default = 20)
    parser.add_option("--tmpDir", dest = "tmp_dir",
                      help = ("Path to directory for tmp files. "
                              "Default = tmp_label_reads"), 
                      default = "tmp_label_reads")
    parser.add_option("--deleteTmp", dest = "delete_tmp",
                      action='store_true',
                      help = ("If this option is set, the temporary directory "
                              "generated by the program will be "
                              "removed at the end of the run."))
    parser.add_option("--o", dest = "outprefix", default = "talon_prelabels",
                      help = "Prefix for outfiles")

    (opts, args) = parser.parse_args()
    return opts

def fetch_seq(chrom: str, start: int, stop: int, strand: str, genome: pyfaidx.Fasta,
              indexing=0):
    """ Given a genomic interval, return the sequence with respect to the
        strand supplied.
        If 1-based indexing is specified, then 1 will be subtracted from the
        position to convert to the Python indexing. """

    if start > stop:
        raise ValueError("Start must be less than or equal to stop")

    if indexing != 0:
        if indexing == 1:
            start -= 1
        else:
            raise ValueError("Valid indexing modes include: 1 or 0")

    seq = genome[chrom][start:stop]

    if strand == "-":
        seq = seq.reverse.complement

    return str(seq)

def compute_frac_As(seq: str):
    """ Compute fraction of sequence made up of As """

    a = seq.count('A')
    n = len(seq)
    if n == 0:
        return 0
    else:
        return float(a)/n

def fetch_range_after_transcript(transcript_end: int, strand: str, length: int):
    """ Given the 1-based stop position of a transcript and its strand,
        return a 1-based genomic range of the specified length that starts with
        the base just after the end position. The smaller position is always
        reported first.
        Example:
              fetch_range_after_transcript(4, '+', 2) would yield (5, 6)
              fetch_range_after_transcript(4, '-', 2) would yield (2, 3)
    """
    if length < 1:
        raise ValueError("Length must be greater than or equal to 1")

    if strand == '+':
        range_start = transcript_end + 1
        range_end = range_start + length - 1
    elif strand == '-':
        range_start = transcript_end - 1
        range_end = range_start - length + 1
    else:
        raise ValueError("Strand must be + or -")

    return (min(range_start, range_end), max(range_start, range_end))

def compute_transcript_end(transcript: pysam.AlignedSegment):
    """ Compute the position of the final transcript base relative to the genome,
        taking strand into account. Position is 1-based. """

    strand = "-" if transcript.is_reverse else "+"
    if strand == '+':
        return transcript.reference_end
    if strand == '-':
        return transcript.reference_start + 1 # (make 1-based)

def compute_frac_as_after_transcript(chrom: str, transcript_end: int, strand: str,
                                     range_size: int, genome: pyfaidx.Fasta):
    """ Given a transcript end, strand, range size, and genome object,
        compute the fraction of sequence in the range immediately after
        the transcript end that is made up of As."""

    # Get sequence of range immediately after transcript
    range_start, range_end = fetch_range_after_transcript(transcript_end,
                                                          strand, range_size)
    range_seq = fetch_seq(chrom, range_start, range_end, strand, genome,
                          indexing = 1)

    # Get fraction As in sequence
    return compute_frac_As(range_seq)

  
def split_reads_by_chrom(sam_file, tmp_dir = "tmp_label_reads", n_threads = 1):
    """ Reads a SAM/BAM file and splits the reads into one file per chromosome.
        Returns a list of the resulting filenames."""

    ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    print("[ %s ] Splitting SAM by chromosome..." % (ts))
    sys.stdout.flush()

    tmp_dir = tmp_dir + "/raw"
    os.system("mkdir -p %s" %(tmp_dir))

    if sam_file.endswith(".sam"):
        # Convert to bam
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print("[ %s ] -----Converting to bam...." % (ts))
        sys.stdout.flush()
        bam_file = tmp_dir + "/all_reads.bam"
        pysam.view("-b", "-S", "-@", str(n_threads), "-o", bam_file, sam_file, 
                   catch_stdout=False)
    elif sam_file.endswith(".bam"):
        bam_file = sam_file
    else:
        raise ValueError("Please provide a .sam or .bam file")

    # Index the file if no index exists
    if not os.path.isfile(bam_file + ".bai"):
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print("[ %s ] -----Sorting and indexing..." % (ts))
        sys.stdout.flush()
        sorted_bam = tmp_dir + "/all_reads.sorted.bam"
        pysam.sort("-@", str(n_threads), "-o", sorted_bam, bam_file)
        bam_file = sorted_bam
        pysam.index(bam_file)
        
    # Open bam file
    tmp_dir += "/chroms"
    os.system("mkdir -p %s" %(tmp_dir))
    read_files = []
    ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    print("[ %s ] -----Writing chrom files..." % (ts))
    sys.stdout.flush()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over chromosomes and write a reads file for each
        chromosomes = [ x.contig for x in bam.get_index_statistics() \
                        if x.mapped > 0 ]
        for chrom in chromosomes:
           records = bam.fetch(chrom)
           fname = tmp_dir + "/" + chrom + ".sam"
           with pysam.AlignmentFile(fname, "w", template = bam) as o: 
               for record in records:
                   o.write(record)
           read_files.append(fname)

    return read_files


def run_chrom_thread(sam_file, options):
    """ """
    outname = sam_file.split("/")[-1].split(".sam")[0]
    genome = pyfaidx.Fasta(options.genome_file, sequence_always_upper=True,
                           one_based_attributes=False)

    os.system("mkdir -p %s" % (options.tmp_dir + "/labeled")) 
    out_log_fname = options.tmp_dir + "/labeled/" + outname + "_read_labels.tsv"
    out_sam_fname = options.tmp_dir + "/labeled/" + outname + ".sam"

    # Iterate over reads
    out_log = open(out_log_fname, 'w')
    pos_seen_fracA = {} # Store fraction As for previously seen positions
    with pysam.AlignmentFile(sam_file) as sam:
        out_sam = pysam.AlignmentFile(out_sam_fname, "w", template=sam)

        for record in sam:  # type: pysam.AlignedSegment
            if record.is_secondary == True or record.is_unmapped == True:
                continue
            read_id = record.query_name
            chrom = record.reference_name
            strand = "-" if record.is_reverse else "+"
            transcript_end = compute_transcript_end(record)

            # Add custom fraction A tag to the read
            location_str = "_".join([chrom, str(transcript_end), strand])
            if location_str in pos_seen_fracA:
                frac_As = pos_seen_fracA[location_str]
            else:
                frac_As = compute_frac_as_after_transcript(chrom, transcript_end, 
                                                           strand,
                                                           options.fracA_range_size,
                                                           genome)
                pos_seen_fracA[location_str] = frac_As

            record.tags += [('fA', round(frac_As,3))]

            # TODO: Add other labels to the read, i.e. CAGE, canonical polyA

            # Write to output files
            out_sam.write(record)
            out_log.write("\t".join([read_id, str(frac_As)]) + '\n')

        out_sam.close()
    out_log.close()
    return

def pool_outputs(indir, outprefix):
    """ Given an input directory containing SAM files and log files,
        concatenate them to form the final output. """

    sam_fname = outprefix + "_labeled.sam"
    log_fname = outprefix + "_read_labels.tsv" 
    
    # Get list of files to combine
    sam_files = glob.glob(indir + "/*.sam")
    log_files = glob.glob(indir + "/*_read_labels.tsv")

    # Add headers
    with open(log_fname, 'w') as f:
        f.write("\t".join(["read_name", "fraction_As"]) + '\n')

    os.system('cp %s %s' % (sam_files[0], sam_fname))

    # Concatenate 
    for sam in sam_files[1:]:
        os.system('grep -v "^@" %s >> %s' % (sam, sam_fname))

    for logfile in log_files:
        os.system('cat %s >> %s' % (logfile, log_fname))

    return

def main(options=None):
    if options == None:
        options = get_options()

    # Initialize worker pool
    with mp.Pool(processes=options.threads) as pool:
        # Print start message
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print("[ %s ] Started talon_label_reads run." % (ts))
        sys.stdout.flush()

        # Remove tmp dir if it exists
        if os.path.exists(options.tmp_dir):
            os.system("rm -r %s" % (options.tmp_dir))

        # Partition reads by chromosome
        read_files = split_reads_by_chrom(options.sam_file, tmp_dir = options.tmp_dir, 
                                          n_threads = options.threads)

        # Now launch the parallel TALON read label jobs
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        print("[ %s ] Launching parallel jobs..." % (ts))
        sys.stdout.flush()
        jobs = [(sam, options) for sam in read_files]
        pool.starmap(run_chrom_thread, jobs)

        pool.close()
        pool.join()

    # Pool together output files
    ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    print("[ %s ] Pooling output files..." % (ts))
    sys.stdout.flush()
    pool_outputs(options.tmp_dir + "/labeled", options.outprefix)

    # Delete tmp_dir if desired
    if options.delete_tmp:
        os.system("rm -r %s" % (options.tmp_dir))
 
    ts = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    print("[ %s ] Run complete" % (ts))
    sys.stdout.flush()

if __name__ == '__main__':
    options = get_options()
    main(options) 
