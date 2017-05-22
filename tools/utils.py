import sys, os
import textwrap
import configparser
from configparser import ConfigParser, ExtendedInterpolation

# Read config File
config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','settings.conf'))

# set -e: If a command fails, set -e will make the whole script exit, instead of just resuming on the next line.
# set -o pipefail: set -o pipefail causes a pipeline (for example, curl -s http://sipb.mit.edu/ | grep foo) to produce a failure return code if any command errors.
# The mktemp utility takes the given filename template and overwrites a portion of it to create a unique filename. 
# Tee command is used to store and view (both at the same time) the output of any other command. 
# Crea una cartella temporanea e si posiziona in questa

def cmd_init():
    cmd_init = textwrap.dedent(r"""
            currentDir=$(pwd) && export CURRENT_DIR=$currentDir;
            tmpDir=$(mktemp -d --tmpdir={config[opt][scratch]}) && export TMPDIR=$tmpDir;
            printf "%s %s at %s\n" "$(date "+%T %D")" "$(hostname)" "$tmpDir" | tee /dev/stderr && cd $tmpDir;
            """.format(config=config))
    return cmd_init

def cmd_out():
    cmd_out  = textwrap.dedent(r"""
            echo ""$(date "+%T %D")" Moving files to Storage"; 
            [[ -a $tmpDir/out.bam ]] && mv -f $tmpDir/out.bam $OUT_FILE;
            OUT_FILE_BAI=$(echo $OUT_FILE | sed -e "s/.bam/.bai/g");
            [[ -a $tmpDir/out.bai ]] && mv -f $tmpDir/out.bai $OUT_FILE_BAI;
            cd && /bin/rm -rf $tmpDir;
            """)
    return cmd_out

def cmd_out_vcf():
    cmd_out_vcf  = textwrap.dedent(r"""
            echo ""$(date "+%T %D")" Moving files to Storage"; 
            [[ -a $tmpDir/out.vcf ]] && mv -f $tmpDir/out.vcf $OUT_FILE;
            OUT_FILE_IDX=$(echo $OUT_FILE".idx");
            [[ -a $tmpDir/out.vcf.idx ]] && mv -f $tmpDir/out.vcf.idx $OUT_FILE_IDX;
            echo ""$(date "+%T %D")" Moving done";
            cd && /bin/rm -rf $tmpDir;
            """)
    return cmd_out_vcf

def awk_cmd():
    return """xargs awk -F'\t' '$1 == "Filename" {{sub("_[0-9]*.fastq",".fastq",$2); sample=$2}}; $1 == "Total Sequences" {{test=$2"\t"$1"\t"sample}}; END {{print test}}' """

def _list2input(l, opt):
    return opt + (" "+opt).join(map(lambda x: str(x), l))

def Load_Fastq(chunk,sample_name,rgid,pair,path,platform,platform_unit,procedure,library,adapter1,adapter2,manifest,assembly,result_dir):
    return

