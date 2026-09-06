"""CPU admission control for one retro run; Slurm supplies available capacity."""
from dataclasses import dataclass
from math import ceil


@dataclass(frozen=True)
class SlurmBudget:
    cpu_limit: int
    worker_cpus: int = 1
    coordinator_cpus: int = 1
    allocation_quantum: int = 1

    def __post_init__(self):
        if self.worker_cpus < 1 or self.coordinator_cpus < 1 or self.allocation_quantum < 1:
            raise ValueError('CPU reservations must be positive')
        object.__setattr__(self, 'worker_cpus',
            ceil(self.worker_cpus / self.allocation_quantum) * self.allocation_quantum)
        if self.cpu_limit < self.coordinator_cpus + self.worker_cpus:
            raise ValueError('CPU budget must fit the coordinator and one worker')

    @property
    def concurrency(self):
        return (self.cpu_limit - self.coordinator_cpus) // self.worker_cpus

    def arguments(self, shards):
        return [f'--array={shards}%{self.concurrency}',
                f'--cpus-per-task={self.worker_cpus}']
