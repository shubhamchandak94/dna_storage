
Flexbar installation:
=====================

Flexbar binaries are available for Linux 64, Mac OSX and Windows 32 and 64 on the web site:
https://sourceforge.net/projects/theflexibleadap

To run binaries, make sure that the TBB library (Intel threading building blocks) is available to the system:

### Linux
One possibility is to simply put the file libtbb.so.2 in your working directory. To use it permanently, copy libtbb.so.2 to the shared library directory with the following command or ask the administrator to install it:
cp FLEXBAR_DIR/libtbb.so.2 /usr/local/lib

Or adjust the lib search path to include the directory of the lib file for the current terminal session:
export LD_LIBRARY_PATH=FLEXBAR_DIR

### Mac OSX
Use the file libtbb.dylib for this system. It applies the same as for Linux. Copy the file to the lib directory:
cp FLEXBAR_DIR/libtbb.dylib /usr/local/lib

Or use this command to include the lib file in the lib search path:
export DYLD_LIBRARY_PATH=FLEXBAR_DIR

### Windows
Leave the file tbb.dll in same directory together with the Flexbar binary. Additionally, Visual Studio 10 sp1 has to be installed. Alternatively, those who have not Visual Studio can download the Visual Studio 10 sp1 redistributable from Microsoft:

Windows 32 (x86):
https://www.microsoft.com/en-us/download/details.aspx?id=8328

Windows 64 (x64):
https://www.microsoft.com/en-us/download/details.aspx?id=13523


Flexbar usage:
==============

Please refer to the help screen or documentation on:
https://sourceforge.net/apps/mediawiki/theflexibleadap

Flexbar needs at least the following input:
1) file in fasta format containing barcode or adapter sequences (b or a)
2) file in fasta/q or csfasta/q format with sequencing reads (s)

Additionally the target name and format of reads have to be set.
See the help screen for more information: flexbar -h

SYNOPSIS
    flexbar -t target -f format -s reads { -b barcodes | -a adapters } [options]

DESCRIPTION
    -h, --help
          Display this help message
    -v, --version
          Display version information
    -n, --threads NUM
          Number of threads, default: 1
    -t, --target STR
          Prefix for output file names
    -s, --source FILE
          Input file with reads, that may contain barcodes
    -p, --source2 FILE
          Second input file for paired read scenario
    -f, --format STR
          Input format of reads: csfasta, csfastq, fasta, fastq-sanger, fastq-solexa, fastq-i1.3,
          fastq-i1.5, fastq-i1.8 (illumina)

  Barcode detection:
    -b, --barcodes FILE
          Fasta file with barcodes, specify (br) to use seperate barcode reads
    -br, --barcode-reads FILE
           Fasta or fastq file with barcode reads, if barcodes not within reads
    -be, --barcode-trim-end STR
           Type of barcoding within source reads, see section trim-end types, default: ANY
    -bt, --barcode-threshold NUM
           Allowed mismatches and indels per 10 bases for barcode, default: 1.0
    -bo, --barcode-min-overlap NUM
           Minimum overlap for barcodes (default is length of first barcode)
    -bm, --barcode-match NUM
           Match score, default: 1
    -bi, --barcode-mismatch NUM
           Mismatch score, default: -1
    -bg, --barcode-gap-cost NUM
           Gap score, default: -7
    -bv, --barcode-remove
           Remove barcodes within reads based on barcoding parameters

  Adapter removal:
    -a, --adapters STR
          Fasta file with adapters, or barcodes to remove within reads
    -as, --adapter-seq STR
           Single adapter sequence as alternative to adapters option
    -ae, --adapter-trim-end STR
           Type of alignment for removal, see section trim-end types, default: RIGHT
    -at, --adapter-threshold NUM
           Allowed mismatches and indels per 10 bases for adapter, default: 3.0
    -ao, --adapter-min-overlap NUM
           Minimum overlap of adapter and read in base pairs, default: 8
    -am, --adapter-match NUM
           Match score, default: 1
    -ai, --adapter-mismatch NUM
           Mismatch score, default: -1
    -ag, --adapter-gap-cost NUM
           Gap score, default: -7
    -aa, --adapter-no-adapt
           Do not treat parameter min-overlap as adaptive measure, see doc

  Filtering:
    -m, --min-readlength NUM
          Minimum read length to remain after removal, default: 18
    -u, --max-uncalled NUM
          Allowed uncalled bases (N or .) in reads, default: 0
    -q, --phred-pre-trim NUM
          Trim reads from 3' end until quality larger than specified
    -x, --fixed-pre-trim NUM
          Trim reads to fixed number of bases before removal
    -y, --fixed-post-trim NUM
          Trim reads to fixed number of bases after removal
    -d, --no-length-dist
          Prevent length distribution for each read output file
    -r, --removal-tag
          Tag reads for which adapter or barcode is removed
    -l, --log-level STR
          Print alignments for all or modified reads. One of ALL, MOD, and TAB

  Trim-end types:
    ANY: longer part of read remains
    LEFT: adapters that align <= read end position
    RIGHT: adapters that align >= read start position
    LEFT_TAIL: considers first n (adapter length) bases
    RIGHT_TAIL: considers last n bases of reads

EXAMPLES
    flexbar -t target -f csfastq -s reads.csfastq -a adapters.fasta
    flexbar -t target -f fastq-i1.3 -s reads.fastq -b bar.fasta -a adap.fasta


In the first example, adapters in fasta format are removed from color-space reads that have quality scores (csfastq format). The resulting reads are written to the file target.csfastq in same format. In the second example, barcoded reads in illumina version 1.3 fastq format are demultiplexed by specifying a file with barcodes in fasta format. After the reads are seperated with respect to barcodes, adapters given in fasta format are removed from the sequencing data.

Although default parameters of Flexbar are optimized to deliver good results in a large number of scenarios, it might improve results, if you change parameters, e.g. --adapter-min-overlap, --adapter-threshold, and --min-readlength.


Building Flexbar from SVN:
==========================

It should be pssible to compile Flexbar for any platform. What you need is:

1) Check out the SVN repository to a local directory FLEXBAR_DIR via the following command or use a graphical SVN client like tortoisesvn.
svn co svn://svn.code.sf.net/p/theflexibleadap/code/trunk FLEXBAR_DIR

2) Download the threading library TBB if you dont have a Linux, Windows or Mac OSX running. For these systems the lib files are supplied with binaries. Download the latest stable release, it should work with version >= 3.0:
http://www.threadingbuildingblocks.org/file.php?fid=77

unpack the archive:
tar xzvf tbb30_20100406oss_lin.tgz

and in the unpacked folder run: gmake

3) Make the TBB library available in your library searchpath:
- If you have Linux 64, Mac OSX or Windows follow the steps described for binaries at the beginning.
- If you compiled TBB yourself copy the compiled lib from your $TBB_INSTALL_DIR/build/release folder to FLEXBAR_DIR/lib and to your library searchpath.

4) Download cmake from http://www.cmake.org and install it.

5) Start cmake via gui and open FLEXBAR_DIR, or use the command line to change to FLEXBAR_DIR and type (do not forget the dot):
cmake .

CMake also allows you to create makefiles or project files for development evironments via the -G switch, e.g. you can build Eclipse project files:
cmake -G "Eclipse CDT4 - Unix Makefiles" .

6) Compile the source code via the make build system by simply typing make in FLEXBAR_DIR. In general, the seqan and tbb library (in FLEXBAR_DIR/lib) need to be available to the compiler and linker. In case of eclipse, this means you should start eclipse, import the project from "FLEXBAR_DIR/build" and compile Flexbar after you set the lib path in your eclipse project settings.


Project files and folders:
==========================

lib: contains shared tbb libs for different platforms (Linux 64, Mac OSX, Windows).
include: contains an adapted version of the seqan and tbb framework.
test: folder contains small test datasets for testing Flexbar after modifications.

Developers test dataset:
Flexbar comes with a test dataset, designed to test most of the adapter removal scenarios. If you change Flexbar sourcecode and want to make sure it does the same as before, just make sure the compiled flexbar binary is found via the path variable and run flexbar_validate.sh from the test folder.
