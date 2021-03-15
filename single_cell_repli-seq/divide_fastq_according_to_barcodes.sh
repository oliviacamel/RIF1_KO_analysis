while getopts hb:o:f: option;
do
case "${option}" in
b) barcodefile=${OPTARG};;
o) outputprefix=${OPTARG};;
f) fastq=${OPTARG};;
h|*) echo "Usage: $0 [-b barcode] [-o outputfile] [-f fastq]" 1>&2; exit 1 ;;
esac
done
if [ $OPTIND -eq 1 ]; then echo "Usage: $0 [-b barcode] [-o outputfile] [-f fastq]" 1>&2; exit 1; fi
shift $((OPTIND -1))

for linenumber in $(seq 1 $(wc -l $barcodefile| awk '{print $1}')); do barcode=$( cat $barcodefile | awk -v field="$linenumber" 'NR==field {print $1}' ) ; cellcycle=$( cat $barcodefile | awk -v field="$linenumber" 'NR==field {print $2}' ) ; cat $fastq  | awk -v barcodeinfile="$barcode" '$1~barcodeinfile {{print f}{print $0}{for(i=1; i<3; i++){getline; print}}}{f=$0} ' > $outputprefix\_$barcode\_$cellcycle\.fastq
done
exit 0
