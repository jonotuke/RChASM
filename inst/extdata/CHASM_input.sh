#!/bin/bash
# KP 2026-01-27

# Set the following variables to fully qualified paths if samtools or bedtools are not in your PATH:
SAMTOOLS=samtools
BEDTOOLS=bedtools
PERL=perl

help()
{
	echo "CHASM_input.sh: Produce input for CHASM from bam files."
	echo "Usage: $0 -i [in-file] -o [out-file] -p [protocol-name] -b [bed-file] -w [+/-bp for bed] -q [min-MQ] -l [min-length]"
	echo "Parameters:"
	echo "     -i [file]     --  text file listing bam files, one per line; run with -i <(ls -1 *.bam) if you like to use bam files from the current directory." 
	echo "     -o [file]     --  output text table for CHASM."
	echo "   [ -p [name]     --  optional protocol name; default: \"default-protocol\" ]"
	echo "   [ -b [bed-file] --  optional bed file with coordinates for filtering ]"
	echo "   [ -w [number]   --  optional number of basepairs to extend the size of regions in the bed file. ]"
	echo "   [ -q [number]   --  optional minimum MQ; default: 25 ]"
	echo "   [ -l [number]   --  optional minimum length in bp; default: 35 ]"
	echo "   [ -N            --  do not print header line for output; default: print headerline ]"
	echo "   [ -h            --  prints this help ]"  
}

# Check that executables are present
if ! `$SAMTOOLS --version >/dev/null 2>&1` ; then
	echo "Error: Can't find samtools. Please install samtools or adjust the path to the samtools executable at the beginning of shell script $0."
	help ;
	exit 1
fi

if ! `$BEDTOOLS --version >/dev/null 2>&1` ; then
	echo "Error: Can't find bedtools. Please install bedtools or adjust the path to the bedtools executable at the beginning of shell script $0."
	help
	exit 1
fi

if ! `$PERL --version >/dev/null 2>&1` ; then
	echo "Error: Can't find perl. Please install perl." # o.O
	help
	exit 1
fi

# Parameters
MQ=25
MINL=35
BED_FILE=
BED_EXTEND=0
BAM_FILES_FILE=
OUTFILE=
PROTOCOL=default-protocol
NOHEAD=no

while getopts "i:o:p:b:w:q:l:Nh" opt; do
	case $opt in
		i)
			if [ -f $OPTARG -a -r $OPTARG ] ; then
				BAM_FILES_FILE=$OPTARG
			else 
				echo "Error: Input file $OPTARG doesn't exist or is not readable."
				exit 2
			fi
			;;
		o)
			if echo -n "" > $OPTARG ; then
				OUTFILE=$OPTARG
			else 
				echo "Error: Output file $OPTARG can't be created."
				exit 2
			fi
			;;
		p)
			PROTOCOL=$OPTARG
			;;
		b)
			if [ -f $OPTARG -a -r $OPTARG ] ; then
				BED_FILE=$OPTARG
			else 
				echo "Error: Input bed-file $OPTARG doesn't exist or is not readable."
				exit 2
			fi
			;;
		w)
			if [[ ! "$OPTARG" =~ ^[0-9]+$ ]] ; then
				echo "Error: Option -w requires a positive integer number."
				exit 2
			fi
			BED_EXTEND=$OPTARG
			;;
		q)
			if [[ ! "$OPTARG" =~ ^[0-9]+$ ]] ; then
				echo "Error: Option -q requires a positive integer number."
				exit 2
			fi
			MQ=$OPTARG
			;;
		l)
			if [[ ! "$OPTARG" =~ ^[0-9]+$ ]] ; then
				echo "Error: Option -l requires a positive integer number."
				exit 2
			fi
			MINL=$OPTARG
			;;
		N)
			NOHEAD=yes
			;;
		h)
			help
			exit 0 ;
			;;
	esac
done

# check basic params
if [ -z "$BAM_FILES_FILE" -o -z "$OUTFILE" ] ; then
	echo "Error: Need at least in-file and out-file as parameters."
	help
	exit 1
fi

if [ "$NOHEAD" == "no" ] ; then
	(echo -en "sample\tprotocol\t" ;
	for i in `seq 1 22` ; do echo -en "chr$i\t" ; done ;
	echo -e "X\tY")>$OUTFILE
fi

# Go through all BAM files and get counts
cat $BAM_FILES_FILE | while read bam ; do
	if ! samtools quickcheck $bam ; then
		echo "Error: bam file $bam does not appear to be a bam file or is broken."
		exit 3
	fi

	# we filter with a bed file
	if [ -s "$BED_FILE" ] ; then
		$BEDTOOLS intersect -a $bam -b <(cat $BED_FILE|perl -F"\t" -lanse 'if (/^#/) { print ; } else { $F[1]-=$L ; $F[2]+=$L ; print join "\t", @F ; }' -- -L=$BED_EXTEND) | $SAMTOOLS view -q $MQ /dev/stdin | $PERL -s -- <(cat <<'ENDL'
			use strict ;
			our ( $LEN, $NAME, $PROT ) ;

			my %chr ;
			for ( my $i = 1 ; $i <= 22 ; $i++ ) { $chr{$i}=$i-1 ; }
			$chr{X} = 22 ;
			$chr{Y} = 23 ;
			#$chr{MT} = 24 ;

			my @res = (0) x keys %chr ;
			while (<STDIN>) {
				chomp ; 
				my @r = split /\t/ ;
				next if ( length($r[9]) < $LEN ) ;
				$res[$chr{$r[2]}]++ if ( exists $chr{$r[2]} ) ;
			}

			print join "\t", ($NAME, $PROT, @res) ; print "\n" ;
						
ENDL
) -LEN=$MINL -NAME=$bam -PROT=$PROTOCOL
		
	else
	# No bed file, just all reads
		$SAMTOOLS view -q $MQ $bam | $PERL -s -- <(cat <<'ENDL' 
			use strict ;
			our ( $LEN, $NAME, $PROT ) ;
			my %chr ;
			for ( my $i = 1 ; $i <= 22 ; $i++ ) { $chr{$i}=$i-1 ; }
			$chr{X} = 22 ;
			$chr{Y} = 23 ;
			#$chr{MT} = 24 ;

			my @res = (0) x keys %chr ;
			while (<STDIN>) {
				chomp ; 
				my @r = split /\t/ ;
				next if ( length($r[9]) < $LEN ) ;
				$res[$chr{$r[2]}]++ if ( exists $chr{$r[2]} ) ;
			}

			print join "\t", ($NAME, $PROT, @res) ; print "\n" ;
						
ENDL
) -LEN=$MINL -NAME=$bam -PROT=$PROTOCOL
	fi 
done >>$OUTFILE

