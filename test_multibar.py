#!/usr/bin/env python
"""
Test script for multi-progress bar implementation.
This will be the template for the actual implementation.
"""

from multiprocessing import Manager, Process
from rich.progress import Progress
import time
import random

def worker(chunk_id, progress_dict, total_work):
    """Simulate chunk processing with progress updates."""
    for i in range(total_work):
        time.sleep(random.uniform(0.01, 0.05))  # Simulate work
        progress_dict[chunk_id] = i + 1

def main():
    num_chunks = 4
    work_per_chunk = 100
    
    with Manager() as manager:
        progress_dict = manager.dict()
        
        # Initialize progress for all chunks
        for i in range(num_chunks):
            progress_dict[i] = 0
        
        # Start worker processes
        processes = []
        for i in range(num_chunks):
            p = Process(target=worker, args=(i, progress_dict, work_per_chunk))
            p.start()
            processes.append(p)
        
        # Monitor progress with individual bars
        with Progress() as progress:
            tasks = {}
            for i in range(num_chunks):
                tasks[i] = progress.add_task(f"Chunk {i}", total=work_per_chunk)
            
            # Update progress bars until all workers are done
            while any(p.is_alive() for p in processes):
                for i in range(num_chunks):
                    current = progress_dict.get(i, 0)
                    progress.update(tasks[i], completed=current)
                time.sleep(0.1)
            
            # Final update
            for i in range(num_chunks):
                progress.update(tasks[i], completed=work_per_chunk)
        
        # Wait for all processes to finish
        for p in processes:
            p.join()

if __name__ == "__main__":
    main()

