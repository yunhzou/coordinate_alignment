"""Process-local AAM profiling without changing search decisions.

CPU accounting is additive across workers; wall durations are NOT subtracted
across parallel processes. Inclusive phase records must not be summed together.
"""
from collections import defaultdict
from functools import wraps
import json
import os
from pathlib import Path
import resource
import time


def cpu_total():
    return sum(r.ru_utime+r.ru_stime for r in
               (resource.getrusage(resource.RUSAGE_SELF),resource.getrusage(resource.RUSAGE_CHILDREN)))


class SearchProfiler:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=False)
        self.parent_pid = os.getpid()
        self.restore = []

    def wrapper(self, function, phase):
        @wraps(function)
        def call(*args, **kwargs):
            wall, cpu = time.perf_counter(), time.process_time()
            try:
                return function(*args, **kwargs)
            finally:
                event = dict(phase=phase, pid=os.getpid(), parent=os.getpid()==self.parent_pid,
                    wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu)
                # One writer per file, inherited safely by forked workers.
                with (self.directory/f'{os.getpid()}.jsonl').open('a') as stream:
                    stream.write(json.dumps(event, separators=(',',':'))+'\n')
        return call

    def patch(self, owner, name, value):
        self.restore.append((owner,name,owner.__dict__[name]))
        setattr(owner,name,value)

    def __enter__(self):
        from rxn_core import aam, artifacts
        from rxn_core.search_graph import AAMSearchGraph
        for name in ('write_raw_cut','write_graph_checkpoint','write_aam_checkpoint'):
            self.patch(artifacts,name,self.wrapper(getattr(artifacts,name),'persistence'))
        for name in ('read_raw_cut','read_graph_checkpoint','read_aam_checkpoint'):
            self.patch(artifacts,name,self.wrapper(getattr(artifacts,name),'loading'))
        self.patch(aam,'_search_cut',self.wrapper(aam._search_cut,'cut_search_inclusive'))
        self.patch(aam,'finalize_graph_symmetry',self.wrapper(aam.finalize_graph_symmetry,'symmetry'))
        combine = AAMSearchGraph.combine.__func__
        self.patch(AAMSearchGraph,'combine',classmethod(self.wrapper(combine,'graph_merge')))
        self.wall_start, self.cpu_start = time.perf_counter(), cpu_total()
        return self

    def __exit__(self, *exception):
        self.wall_seconds = time.perf_counter()-self.wall_start
        self.cpu_seconds = cpu_total()-self.cpu_start
        for owner,name,value in reversed(self.restore):
            setattr(owner,name,value)

    def summary(self):
        phases = defaultdict(lambda: dict(calls=0, summed_wall_seconds=0., cpu_seconds=0.))
        for path in self.directory.glob('*.jsonl'):
            for line in path.read_text().splitlines():
                event = json.loads(line)
                key = event['phase']+('/parent' if event['parent'] else '/workers')
                row = phases[key]
                row['calls'] += 1
                row['summed_wall_seconds'] += event['wall_seconds']
                row['cpu_seconds'] += event['cpu_seconds']
        io_cpu = sum(v['cpu_seconds'] for k,v in phases.items()
                     if k.startswith(('persistence/','loading/')))
        assert self.cpu_seconds+1e-5 >= io_cpu
        return dict(elapsed_wall_including_io_seconds=self.wall_seconds,
            total_cpu_including_io_seconds=self.cpu_seconds,
            compute_cpu_excluding_persistence_and_loading_seconds=self.cpu_seconds-io_cpu,
            phases=dict(phases),
            semantics='CPU excludes measured checkpoint encoding/compression/writing and decoding/loading. '
                      'Includes IPC, process setup, and profiling overhead. Phase wall values are summed worker durations, '
                      'not parallel latency. cut_search includes within-cut graph merges; do not sum inclusive phases. '
                      'No invalid subtraction of summed worker IO from elapsed wall time.')
