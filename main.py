"""
Program: DNA-Seq Germline Pipeline 
Contact: <info@engenome.com>

Input file is a json of the following format:
[
  {
    "chunk"          : "001",           
    "sample_name"    : "F1-Sample",         # Single biological sample
    "rgid"           : "F1-RG",             # Read Group assigned to all samples belonging to the same patient.   
    "pair"           : "1",                 # 1 or 2
    "path"           : "/path/to/fastq",
    "platform"       : "Illumina",          # Illumina, PacBio, ThermoFischer
    "platform_unit"  : "FlowCellId",    
    "procedure"      : "WGS",               # WGS, Capture, Amplicon
    "library"        : "/path/to/bed",      # Target regions BED file 
    "adapter1"       : "CTGTCTCTTGATCACA",  # Trimming
    "adapter2"       : "CTGTCTCTTGATCACA",  # Trimming
    "manifest"       : "/path/to/manifest", # Trimming: Illumina Manifest file
    "assembly"       : "b37",               # b37, hg38
    "result_dir"     : "/path/to/result_dir"# Final BAM and VCF folder 
  },
  {..}
]
"""
import sys, os
import json
import configparser
import subprocess as sp
from itertools import groupby
from functools import partial
from argparse import RawTextHelpFormatter
from configparser import ConfigParser, ExtendedInterpolation
from cosmos.api import Cosmos, Dependency, draw_stage_graph, draw_task_graph, pygraphviz_available, default_get_submit_args

def settings(svr,output_dir):
    # Set settings.conf file 
    config = configparser.SafeConfigParser()
    config_file=os.path.join(sys.path[0],'settings.conf')
    config.read(config_file)
    if svr=="aws":
        config.set('opt', 'scratch', '${aws:scratch}')
        config.set('opt', 'share_ref', '${aws:share_ref}')
        config.set('opt', 'share_tools', '${aws:share_tools}')
    elif svr=="development":
        config.set('opt', 'scratch', '${development:scratch}')
        config.set('opt', 'share_ref', '${development:share_ref}')
        config.set('opt', 'share_tools', '${development:share_tools}')
    elif svr=="production":
        config.set('opt', 'scratch', '${production:scratch}')
        config.set('opt', 'share_ref', '${production:share_ref}')
        config.set('opt', 'share_tools', '${production:share_tools}')
    config.write(open(config_file,'w'))

def get_drmaa_native_specification(task, default_queue='all.q', parallel_env='orte'):
    """
    Default method for determining the extra arguments to pass to the DRM.
    For example, returning `"-n 3" if` `task.drm == "lsf"` would cause all jobs
    to be submitted with `bsub -n 3`.
    :param cosmos.api.Task task: The Task being submitted.
    :param str default_queue: The default_queue to submit to
    :rtype: str
    """
    drm = task.drm
    default_job_priority = None
    use_mem_req = True
    use_time_req = False
    core_req = task.core_req
    mem_req = task.mem_req if use_mem_req else None
    time_req = task.time_req if use_time_req else None
    jobname = '%s[%s]' % (task.stage.name, task.uid.replace('/', '_'))
    default_queue = ' -q %s' % default_queue if default_queue else ''
    priority = ' -p %s' % default_job_priority if default_job_priority else ''
    
    if drm in ['lsf', 'drmaa:lsf']:
        rusage = '-R "rusage[mem={mem}] ' if mem_req and use_mem_req else ''
        time = ' -W 0:{0}'.format(task.time_req) if task.time_req else ''
        return '-R "{rusage}span[hosts=1]" -n {task.core_req}{time}{default_queue} -J "{jobname}"'.format(**locals())
    elif drm in ['ge', 'drmaa:ge']:
        '''        
        h_vmem = int(math.ceil(mem_req / float(core_req))) if mem_req else None
    
        def g():
            resource_reqs = dict(h_vmem=h_vmem, slots=core_req, time_req=time_req)
            for k, v in resource_reqs.items():
                if v is not None:
                    yield '%s=%s' % (k, v)

        resource_str = ','.join(g())
        '''
        return '{queue} -l spock_mem={mem_req}M,num_proc={cpu_req} {priority} -N "{jobname}"'.format(queue=default_queue, mem_req=mem_req, cpu_req=core_req,jobname=jobname,priority=priority)
        #return '-cwd -pe {parallel_env} {core_req} {priority} -N "{jobname}"{queue}'.format(resource_str=resource_str,priority=priority,queue=default_queue,jobname=jobname, core_req=core_req,parallel_env=parallel_env)

    elif drm == 'local':
        return None
    else:
        raise Exception('DRM not supported: %s' % drm)

def recipe(workflow, input_file, output_dir):
    # Read config File
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(os.path.join(sys.path[0],'settings.conf'))
      
    # Load Fatsq files
    input_json = json.load(open(input_file,'r'))
    load_fastq=[]
    for input in input_json:
        load_fastq.append(workflow.add_task(
                                    func=utils.Load_Fastq,
                                    params=dict(
                                            chunk=input['chunk'],
                                            sample_name=input['sample_name'], 
                                            rgid=input['rgid'],
                                            pair=input['pair'], 
                                            path=input['path'],
                                            platform=input['platform'],
                                            platform_unit=input['platform_unit'],
                                            procedure=input['procedure'],
                                            library=input['library'],
                                            adapter1=input['adapter1'],
                                            adapter2=input['adapter2'],
                                            manifest=input['manifest'],
                                            assembly=input['assembly'],
                                            result_dir=input['result_dir']),                                        
                                    uid="_".join(["LoadFastq",input['sample_name'], input['chunk'],input['pair']])))
    
    # Alignment
    get_chunk = lambda task: {'chunk':task.params['chunk'], 
                              'sample_name':task.params['sample_name'], 
                              'rgid':task.params['rgid'], 
                              'platform':task.params['platform'], 
                              'platform_unit':task.params['platform_unit'], 
                              'procedure':task.params['procedure'], 
                              'library':task.params['library'], 
                              'assembly':task.params['assembly'], 
                              'result_dir':task.params['result_dir']}
    
    fastq_task_group={}; align_fastq=[]; check_fastq=[]; count=1
    for tags, load_tasks_group in groupby(load_fastq, get_chunk):
        fastq_task_group[count]=load_tasks_group 
        # fastq_task_group[count] = itertools._grouper object
        # tags = Tags dict
        # parent = Task list
        parent=[]                                   
        for element in fastq_task_group[count]:
            parent.append(element)
        """
        stage_name="FastQC"
        check_fastq.append(workflow.add_task(
                                    func=fastqc.Fastqc,
                                    stage_name=stage_name,
                                    params=dict(core_req=config['requirements']['aligment_core_req'], 
                                                mem_req=int(config['requirements']['aligment_mem_req'])*1024, 
                                                in_txts=[Dependency(t, 'path') for t in parent],
                                                out_dir='%s/%s/' %  (output_dir, stage_name ),
                                                out_subdir='%s/%s/%s/%s/' %  (output_dir, stage_name,tags['sample_name'], tags['chunk']),
                                                out_summary='%s/%s/%s_%s.summary.txt' %  (output_dir, stage_name,tags['sample_name'], tags['chunk'])),
                                    parents=parent,
                                    uid="_".join([stage_name,tags['sample_name'], tags['chunk']])))
        """
        stage_name="Alignment"
        align_fastq.append(workflow.add_task(
                                    func=alignment.Alignment,
                                    stage_name=stage_name,
                                    params=dict(core_req=config['requirements']['aligment_core_req'], 
                                                mem_req=int(config['requirements']['aligment_mem_req'])*1024, 
                                                in_txts=[Dependency(t, 'path') for t in parent],
                                                out_bam='%s/%s/%s/%s.bam' %  (output_dir, stage_name, tags['sample_name'], tags['chunk']),
                                                **tags),
                                    parents=parent,
                                    uid="_".join([stage_name,tags['sample_name'], tags['chunk']])))
    """
    stage_name="FastQC_Statistic"
    workflow.add_task(
                    func=fastqc.FastqcStatistic,
                    stage_name=stage_name,
                    params=dict(core_req=1,
                                mem_req=int(config['requirements']['aligment_mem_req'])*1024, 
                                in_summary=[Dependency(t, 'out_summary') for t in check_fastq],
                                in_dir=[Dependency(t, 'out_dir') for t in check_fastq],
                                out_dir='%s/%s/' %  (output_dir, stage_name),
                                out_html='%s/%s/summary.html' %  (output_dir, stage_name)
                                ),
                                uid="_".join([stage_name]))
    """
    # Group By rgid, sorted by rgid and sample_type
    get_sample = lambda task: {'sample_name':task.params['sample_name'], 
                              'rgid':task.params['rgid'], 
                              'procedure':task.params['procedure'], 
                              'library':task.params['library'], 
                              'assembly':task.params['assembly'], 
                              'result_dir':task.params['result_dir']}
    bam_aligned={}
    aligned=[]
    coverage=[]
    for tags, aligned_samples in groupby(align_fastq, get_sample):
        bam_aligned[count]=aligned_samples
        parent=[]                                   
        for element in bam_aligned[count]:
            parent.append(element)
        stage_name="Mark_Duplicates"
        duplicate_task = workflow.add_task(
                                func=picard.MarkDuplicates,
                                stage_name=stage_name,
                                params=dict(core_req=config['requirements']['picard_markduplicates_core_req'], 
                                            mem_req=int(config['requirements']['picard_mem_req'])*1024, 
                                            in_bams=[Dependency(t, 'out_bam') for t in parent],
                                            out_bam = '%s/%s/%s.bam' % (output_dir, stage_name, tags['sample_name']),
                                            **tags),
                                parents=parent,
                                uid="_".join(tags['sample_name']))
        aligned.append(duplicate_task)
        """ 
        stage_name="Coverage"
        coverage_task = workflow.add_task(
                                func=pipes.Coverage,
                                stage_name=stage_name,
                                params=dict(core_req=1,
                                            mem_req=int(config['requirements']['bedtools_mem_req'])*1024, 
                                            in_bam=[Dependency(duplicate_task, 'out_bam')],
                                            out_count='%s/%s/%s/%s.count' %  (output_dir, stage_name, tags['rgid'], tags['rgid']),
                                            rgid=tags['rgid']
                                            ),
                                uid="_".join([tags['rgid']]))
        coverage.append(coverage_task)
        """
        
    # Indel Realignment and Base Quality Score Recalibration
    intervals = (range(1,23) + ['X', 'Y'])
    #intervals = [20, 21]
    
    bam_aligned_md={}
    gatk_preprocessed=[]
    for n in intervals:
        for tags, aligned_samples in groupby(aligned, get_sample):
            bam_aligned_md[count]=aligned_samples
            parent=[]                                   
            for element in bam_aligned_md[count]:
                parent.append(element)
            stage_name="Indel_Realignment"
            indel_task = workflow.add_task(
                                    func=gatk.IndelRealigner,
                                    stage_name=stage_name,
                                    params=dict(core_req=config['requirements']['indelRealign_core_req'], 
                                                mem_req=int(config['requirements']['indelRealign_mem_req'])*1024, 
                                                in_bam=[Dependency(t, 'out_bam') for t in parent],
                                                out_bam='%s/%s/%s/%s.bam' %  (output_dir, stage_name, tags['sample_name'], n),
                                                chr=n, 
                                                **tags),
                                    parents=parent,
                                    uid="_".join([stage_name,tags['sample_name'], str(n)]))
            stage_name="Haplotype_Caller"
            gatk_preprocessed.append(workflow.add_task(
                                    func=gatk.HaplotypeCaller,
                                    stage_name=stage_name,
                                    params=dict(core_req=config['requirements']['haploCaller_core_req'], 
                                                mem_req=int(config['requirements']['haploCaller_mem_req'])*1024, 
                                                in_bam=[Dependency(indel_task, 'out_bam')],
                                                out_vcf='%s/%s/%s/%s.vcf' %  (output_dir, stage_name, tags['sample_name'], n),
                                                chr=n,
                                                **tags),
                                    uid="_".join([stage_name,tags['sample_name'],str(n)])))

            


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description=__doc__,formatter_class=RawTextHelpFormatter)
    p.add_argument('-i', '--input_json', type=str, required=True, help='Input Json')
    p.add_argument('-o', '--output_dir', type=str, required=True, help='Directory for Intermediate files')
    p.add_argument('-n', '--analysis', default='DNA-Seq_Germline', type=str, help="Analysis ID")
    p.add_argument('-drm', default='drmaa:ge', choices=('local', 'drmaa:ge', 'ge'), help='Distributed resource management system')
    p.add_argument('-svr', default='development', choices=('aws', 'development', 'production'), help='Server settings')
    # Display help message with python argparse when script is called without any arguments
    if len(sys.argv)==1:
        p.print_help()
        sys.exit(1)
        
    args = p.parse_args()
    print args
    
    sp.check_call('mkdir -p '+os.path.join(args.output_dir,args.analysis), shell=True)
    cosmos = Cosmos('sqlite:///%s/%s.db' % (os.path.dirname(os.path.join(args.output_dir,'')),args.analysis),
                    # example of how to change arguments if you're NOT using default_drm='local'
                    #get_submit_args=partial(default_get_submit_args, parallel_env='orte'),
                    get_submit_args=get_drmaa_native_specification,
                    default_drm=args.drm)
    cosmos.initdb()
    os.chdir(os.path.join(args.output_dir,args.analysis))
    
    settings(args.svr, args.output_dir)
    from tools import fastqc, alignment, picard, utils, gatk
    
    workflow = cosmos.start(args.analysis,
                        restart=False,
                        skip_confirm=True)
    recipe(workflow, args.input_json, os.path.join(args.output_dir,args.analysis))
    
    workflow.make_output_dirs()
    workflow.run(max_attempts=4)

    if pygraphviz_available:
        # These images can also be seen on the fly in the web-interface
        draw_stage_graph(workflow.stage_graph(), '/tmp/ex1_task_graph.png', format='png')
        draw_task_graph(workflow.task_graph(), '/tmp/ex1_stage_graph.png', format='png')
    else:
        print 'Pygraphviz is not available :('
