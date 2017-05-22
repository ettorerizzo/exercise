import sys, os
import textwrap
import configparser
from configparser import ConfigParser, ExtendedInterpolation
from utils import cmd_init, cmd_out, _list2input

# Read config File
config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','settings.conf'))

def MarkDuplicates(core_req, mem_req, in_bams, rgid , out_bam, **tags):
        
    cmd_main = textwrap.dedent(r"""
        outFile={out_bam} && export OUT_FILE=$outFile; 
        
        {config[tools][java_MD]} -Xmx{mem_req}M -jar {config[tools][picard_dir]}/MarkDuplicates.jar \
        TMP_DIR=$tmpDir \
        OUTPUT=out.bam \
        METRICS_FILE=out.metrics \
        ASSUME_SORTED=True \
        CREATE_INDEX=True \
        MAX_RECORDS_IN_RAM=1000000 \
        MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=1000 \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        VALIDATION_STRINGENCY=SILENT \
        VERBOSITY=INFO \
        {inputs}
        """.format(config = config,
                   inputs = _list2input(in_bams,"INPUT="),
                   **locals()))
    
    return cmd_init()+cmd_main+cmd_out()