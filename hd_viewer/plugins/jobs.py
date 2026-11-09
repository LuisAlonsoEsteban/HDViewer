"""Process-wide plugin jobs: identical requests share one run, sessions own them."""

from collections import Counter
from contextlib import contextmanager
import threading


class Cancelled(RuntimeError):
    """Every session awaiting a job closed or cancelled it."""

SESSION_CLOSED = "the browser session closed"
USER_CANCELLED = "cancelled by the user"


class Job:
    """One cache key's run. Cancelled when every session awaiting it has closed
    or withdrawn (the user cancelled it)."""

    def __init__(self, key):
        self.key = key
        self.lock = threading.Lock()  # held while checking the cache and running
        self.owners = Counter()  # owner -> threads of that owner awaiting the job
        self.closed = set()  # owners whose session ended or who withdrew
        self.cancelled = threading.Event()
        self.reason = SESSION_CLOSED
        self._kill = None
        self._kill_lock = threading.Lock()

    def live_owners(self):
        return {owner for owner, count in self.owners.items() if count and owner not in self.closed}

    def on_cancel(self, kill):
        """Register how to stop the running process (None to clear)."""
        with self._kill_lock:
            self._kill = kill
        if kill is not None and self.cancelled.is_set():
            kill()

    def check(self):
        if self.cancelled.is_set():
            raise Cancelled(f"cancelled: {self.reason}")

    @contextmanager
    def hold(self, semaphore, poll=0.5):
        """Acquire a run slot, giving up the wait if the job is cancelled."""
        while not semaphore.acquire(timeout=poll):
            self.check()
        try:
            self.check()
            yield
        finally:
            semaphore.release()

    def cancel(self, reason=SESSION_CLOSED):
        if not self.cancelled.is_set():
            self.reason = reason
        self.cancelled.set()
        with self._kill_lock:
            kill = self._kill
        if kill is not None:
            kill()


class JobTable:
    def __init__(self):
        self._lock = threading.Lock()
        self._jobs = {}
        self._closed = set()  # owners whose session ended (a few bytes each)

    @contextmanager
    def join(self, key, owner=None):
        """The job for `key`, shared with concurrent identical requests."""
        with self._lock:
            if owner is not None and owner in self._closed:
                raise Cancelled("cancelled: the browser session closed")
            job = self._jobs.get(key)
            if job is None or job.cancelled.is_set():
                # A cancelled job is still dying: later requests start afresh
                job = self._jobs[key] = Job(key)
            job.owners[owner] += 1
            job.closed.discard(owner)  # rejoining after withdrawing from it
        token = _current.__dict__.get("job")
        _current.job = job
        try:
            yield job
        finally:
            _current.job = token
            with self._lock:
                job.owners[owner] -= 1
                if not +job.owners and self._jobs.get(key) is job:
                    del self._jobs[key]

    def cancel_owner(self, owner):
        """A session closed: cancel the jobs no other live session awaits."""
        with self._lock:
            self._closed.add(owner)
            orphans = []
            for job in set(self._jobs.values()):
                if job.owners.get(owner):
                    job.closed.add(owner)
                    if not job.live_owners():
                        orphans.append(job)
        for job in orphans:
            job.cancel()
        return len(orphans)

    def withdraw(self, job, owner, reason=USER_CANCELLED):
        """`owner` no longer awaits `job` (the user cancelled it): the run is
        stopped unless another live session still awaits it. True if stopped."""
        with self._lock:
            if not job.owners.get(owner):
                return False
            job.closed.add(owner)
            orphan = not job.live_owners()
        if orphan:
            job.cancel(reason)
        return orphan

    def __len__(self):
        with self._lock:
            return len(self._jobs)


_current = threading.local()


def current_job():
    """The job the calling thread runs, if any (runners use it to register kills)."""
    return getattr(_current, "job", None)


# Shared by segmentation and deconvolution (keys are namespaced by service)
JOBS = JobTable()
