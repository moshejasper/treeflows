"""Functions for submitting jobs to spartan"""

### IMPORTS ###

from pathlib import Path
import inspect
#import logging
#from time import localtime, strftime
import subprocess
import random
import uuid
import textwrap
import argparse
import inspect
#from datetime import datetime

def caller():
    """Return the caller module `__name__` from the stack.

    Returns:
        The `__name__` string of the caller's globals (or None if unavailable).
    """
    # get context
    frame = inspect.stack()[2].frame # apparently the level we need
    caller_globals = frame.f_globals
    return caller_globals.get("__name__")

def unique_name():
    """Generate a short unique identifier suitable for job names."""
    return uuid.uuid4().hex[:8] # possible use name main loop here... 

def quickparse(*args):
    """NOTE: this function, despite appearances, does not actually return the datatype... 
    Fix in future via parsing & examining args internally. """
    parser = argparse.ArgumentParser()
    for arg in args:
        atype = int
        try:
            int(arg)
        except ValueError:
            atype = float
            try:
                float(arg)
            except ValueError:
                atype = str
        parser.add_argument(arg, type=atype)
    return parser.parse_args()

class SlurmJobRun:
    """Holds individual Slurm job Run """
    def __init__(self, slurmstring):
        """Create a runnable SLURM submission payload (raw sbatch script text)."""
        self.slurmstring = slurmstring

    def submit(self):
        """Submit this job via `sbatch` and return the job id string."""
        result = subprocess.run(['sbatch', '--parsable'], input=self.slurmstring, capture_output=True, text=True)
        if result.stderr is not None and len(result.stderr.strip()) is not None and len(result.stderr.strip()) > 0:
            print("STDERR:")
            print(result.stderr.strip())
        return result.stdout.strip() 
    
    def submit_conditional(self, jobids, afterany=True):
        """Submit this job with dependencies on `jobids` (afterany/afterok)."""
        d_mode = "afterany" if afterany else "afterok"
        if not type(jobids) == list:
            jobids = [jobids,]
        print(jobids)
        if len(jobids) == 0:
            print("No dependencies: submitting as normal job")
            return self.submit()
        jobids = ":".join([str(x).strip() for x in jobids])
        jobdepends = f"{d_mode}:" + jobids
        print(jobdepends)
        result = subprocess.run(['sbatch', '--parsable', '-d', jobdepends], input=self.slurmstring, capture_output=True, text=True)
        if result.stderr is not None and len(result.stderr.strip()) is not None and len(result.stderr.strip()) > 0:
            print("STDERR:")
            print(result.stderr.strip())
        return result.stdout.strip()
    
    def __str__(self):
        """Return the underlying SLURM script text."""
        return self.slurmstring



class SlurmJob:
    """Class for running and launching slurm jobs. 
    PARAMETERS: 
        pyfile - Path to python script to execute. If `None` assigns currents script and (by default) sets runblock to `__slurm__`. 
                    if a file, sets runblock by default to `__main__`
        threads, threadmem(GB), jobtime('D-HH:MM:SS'), jobname - parameters for SLURM script
        use_tmp, tmp_mem - whether to ask for space on the /tmp/ directory for the job, and how much (GB) to ask for
        modules - list of spartan modules to load (if needed by job)
        conda_env - conda environment to load (appended to default conda path)
        partitions, partitiondict - partitions and part-to-qos mapping (if it needs to change)
        runblock - the __name__ that script is run under (primarily to dodge __main__ or override default)
        slurmdir - an alternative slurm directory to store log files (either relative or absolute)
        
    FUNCTIONS: 
        render_slurm - updates internal str represenation of slurm file and returns it to user (params affect render but not SlurmJob object)
        submit_slurm - renders, updates, & submits slurm file to slurm (params supplied affect render/submit but leave SlurmJob otherwise unaltered)
        """

    def __init__(self, pyfile=None, threads=1, threadmem=10, jobtime="0-12:00:00", modules = ['GSL/2.7'], conda_env = 'gendb', 
                 partitiondict = None, tmp_mem = None, use_tmp = False, partitions = None, jobname="default", runblock=None, 
                 slurmdir = None):
        """Configure a SLURM job wrapper for running a Python script on Spartan."""
        if pyfile is None:
            self.pyfile = Path(inspect.stack()[1].filename).resolve()
            self.runblock = "__slurm__" if runblock is None else runblock
        else:
            self.pyfile = Path(pyfile)
            self.runblock = "__main__" if runblock is None else runblock
        self.threads = threads
        self.threadmem = threadmem # Gb / thread
        self.jobtime = jobtime
        self.modules = list(modules) # enable to be just string if one entry
        self.conda_env = conda_env
        if partitiondict is None:
            self.partitiondict = {'fos': 'fos', 'cascade': 'normal', 'sapphire': 'normal', 'bigmem': 'normal', 'cascade,sapphire': 'normal'}
        else: self.partitiondict = partitiondict
        self.tmp_mem = 10 if tmp_mem is None else tmp_mem
        self.use_tmp = use_tmp
        self.partitions = ['fos', 'cascade,sapphire', 'cascade,sapphire'] if partitions is None else [partitions] if type(partitions) == str else partitions
        self.jobname = jobname
        self.slurmdir = slurmdir
        self.slurmstring = self.update_slurmstring(self.threads, self.threadmem, self.jobtime, self.jobname, self.conda_env, None, self.runblock, self.slurmdir)

    def _prewrapper(self, cpus, mem, jtime, jobname, conda_env, slurmdir):
        """Render the `#!/bin/bash` + `#SBATCH` preamble for a job submission."""

        jname = jobname
        xx = len(self.partitions) - 1 # 2 if all here... 
        idx = random.randint(0, xx) # going to submit across three servers, essentially... 
        jpart = self.partitions[idx]
        qqos = self.partitiondict[jpart]
        if len(self.modules) > 0:
            modscript = 'module load ' + ' '.join(self.modules)
        else:
            modscript = ''
        tmpscript = f"#SBATCH --tmp={self.tmp_mem}GB" if self.use_tmp else ""
        if slurmdir is not None:
            slurmdir = Path(slurmdir)
            if not slurmdir.is_absolute(): # assume is relative
                slurmdir = f"./{slurmdir}"
        sl_dir = f"#SBATCH -o {slurmdir}/slurm-%j.out" if slurmdir is not None else ""
        twrap = "#!/bin/bash"
        twrap += textwrap.dedent(f"""
            #SBATCH --partition={jpart}
            #SBATCH --qos={qqos}

            # Multithreaded (SMP) job: must run on one node 
            #SBATCH --nodes=1
            #SBATCH --job-name="{jname}"
            #SBATCH --account="punim1778"

            #SBATCH --ntasks=1
            #SBATCH --cpus-per-task={cpus}
            #SBATCH --mem-per-cpu={mem}G
            {tmpscript}
            {sl_dir}

            # The maximum running time of the job in days-hours:mins:sec
            #SBATCH --time={jtime}

            # check that the script is launched with sbatch
            if [ "x$SLURM_JOB_ID" == "x" ]; then
            echo "You need to submit your job to the queuing system with sbatch"
            exit 1
            fi

            module load unimelb-mf-clients
            module load Miniconda3
            {modscript}
            eval "$(conda shell.bash hook)"
            conda activate /data/gpfs/projects/punim1778/conda_envs/{conda_env}
            """)
        return(twrap)
    
    def _postwrapper(self):
        """Return a small footer appended to submitted jobs (stats command)."""
        return("\n\nmy-job-stats -a -n -s")
    
    def update_slurmstring(self, threads=None, threadmem=None, jobtime=None, jobname=None, conda_env=None, argstring=None, runblock=None, slurmdir=None):
        """Render and store the slurm script inside this object."""
        slurmstring = self.render_slurm(threads, threadmem, jobtime, jobname, conda_env, argstring, runblock, slurmdir)
        self.slurmstring = slurmstring
    
    def render_slurm(self, threads=None, threadmem=None, jobtime=None, jobname=None, conda_env=None, argstring=None, runblock=None, slurmdir=None):
        """Render a `SlurmJobRun` with current/default parameters."""
        threads = self.threads if threads is None else threads
        threadmem = self.threadmem if threadmem is None else threadmem
        jobtime = self.jobtime if jobtime is None else jobtime
        jobname = self.jobname if jobname is None else jobname
        conda_env = self.conda_env if conda_env is None else conda_env
        argstring = "" if argstring is None else argstring
        runblock = self.runblock if runblock is None else runblock
        slurmdir = self.slurmdir if slurmdir is None else slurmdir
        # py_temp_loc = py_temp + subfile + ".py"
        # py_temp_add = py_temp_loc.replace(':', r'\:')
        slurmstring = self._prewrapper(threads, threadmem, jobtime, jobname, conda_env, slurmdir)
        slurmstring += textwrap.dedent(f"""
            python -c "
            import runpy
            runpy.run_path('{str(self.pyfile)}', run_name='{runblock}')
            " {argstring}
            """)
        slurmstring += self._postwrapper()
        return SlurmJobRun(slurmstring)

    def submit_slurm(self, threads=None, threadmem=None, jobtime=None, jobname=None, conda_env=None, argstring=None, runblock=None, slurmdir=None):
        """Submit Slurm Job. Params override defaults without updating them. 
        Params: 
            threads     - number of threads to use
            threadmem   - memory per thread (GB)
            jobtime     - job running time (str) - format 'D-HH:MM:SS'
            jobname     - job name (str)
            conda_env   - conda environment to load (name, not path),  (str)
            argstring   - commandline arguments to pass to submitted python script "--arg1 1 --arg2 2" etc. 
            runblock    - __name__ block to run within file
            slurmdir    - subfolder or absolute path where log files should be stored
            """
        slurmstring = self.render_slurm(threads, threadmem, jobtime, jobname, conda_env, argstring, runblock, slurmdir)
        result = slurmstring.submit()
        return result
    
    def submit_slurm_conditional(self, jobids, threads=None, threadmem=None, jobtime=None, jobname=None, 
                                 conda_env=None, argstring=None, runblock=None, afterany=True, slurmdir=None):
        """Submit Slurm Job. Params override defaults without updating them. 
        Params: 
            threads     - number of threads to use
            jobids      - list of slurm jobs ids (4353636) to make new job conditional on completing
            threadmem   - memory per thread (GB)
            jobtime     - job running time (str) - format 'D-HH:MM:SS'
            jobname     - job name (str)
            conda_env   - conda environment to load (name, not path),  (str)
            argstring   - commandline arguments to pass to submitted python script "--arg1 1 --arg2 2" etc. 
            runblock    - __name__ block to run within file
            slurmdir    - subfolder or absolute path where log files should be stored
            """
        slurmstring = self.render_slurm(threads, threadmem, jobtime, jobname, conda_env, argstring, runblock, slurmdir)
        result = slurmstring.submit_conditional(jobids, afterany)
        return result
    
    def get_jobids(self, jobname=None, user=None, partition=None):
        """Query `squeue` and return job ids (optionally filtered by jobname/user/partition)."""
        usein = "--me" if user is None else f"-u {user}"
        jraw = subprocess.run(f"squeue {usein}", shell=True, text=True, stdout=subprocess.PIPE).stdout
        titles = ['jobid', 'partition', 'jobname', 'user', 'jobstatus', 'jobtime', 'nnodes', 'nodelist']
        idlist = []
        for jline in jraw.strip().split('\n')[1:]:
            jdict = dict(zip(titles, jline.strip().split()))
            if jobname is not None:
                if jdict['jobname'] != jobname:
                    continue
            if partition is not None:
                if type(partition) == str:
                    partition = [partition,]
                if jdict['partition'] not in partition:
                    continue
            else:
                if jdict['partition'] == "interactive":
                    continue
            idlist.append(jdict['jobid'])
        return idlist
    
    def quicksubmit(self, trigger="__main__", term=True):
        """Submit the job when called from a particular `__name__` block."""
        result = None
        if caller() == trigger:
            result = self.submit_slurm()
            if term:
                exit(f"Submitted slurmjob {result}")
        return result

    def __str__(self):
        """Return a human-readable description of this job configuration."""
        ostring = textwrap.dedent(f"""
            SLURM SUBMIT OBJECT:
            "-----------------"
            Py Script:  {self.pyfile}
            JOB NAME:   {self.jobname}
            JOB TIME:   {self.jobtime}
            THREADS:    {self.threads}
            MEM/THREAD: {self.threadmem}
            TEMPDIR:    {self.use_tmp}
            TEMP_MEM:   {self.tmp_mem}
            """)
        return ostring



if __name__ == "__main__":
    my_slurm = SlurmJob(pyfile=__file__, use_tmp=True)
    print(my_slurm.render_slurm())

if __name__ == "__slurm__":
    print("Hello World")
